#!/N/soft/mason/python/2.7.3/bin/python
# Simo Zhang

import re
import os
import sys

from io_simo import check_files_if_exist
from io_simo import make_dir_if_necessary
from ftsensor import FtSensor

class FixPE(object):
	def __init__(self, args):
		self.fq1 = args.fq1
		self.fq2 = args.fq2
		check_files_if_exist(self.fq1, self.fq2)
		self._check_ft(self.fq1, self.fq2)
		self.outp = os.path.realpath(args.outp)
		make_dir_if_necessary(os.path.dirname(self.outp))
		self._print_cmd()

	def _print_cmd(self):
		sys.stdout.write("[CMD PARSER] FQ1: %s\n" %(self.fq1))
		sys.stdout.write("[CMD PARSER] FQ2: %s\n" %(self.fq2))
		sys.stdout.write("[CMD PARSER] Output Prefix: %s\n" %(self.outp))

	def _check_ft(self, *files):
		''' check input file format '''
		for file in files:
			if FtSensor(file).getfiletype() != "fq":
				sys.stderr.write("[FMT CHECKER] %s is not FastQ format\n" %(file))
				sys.exit(1)

	def _get_splitby(self, file):
		hfmt, _ = FtSensor(file).getencoding()
		if hfmt == "illumina14":
			return "/"
		elif hfmt == "casava18":
			return " "

	def start(self):
		self._fix_pairing()

	def _fix_pairing(self):
		id, str_except_id = "", ""
		fq1_dict = {}
		split_by = self._get_splitby(self.fq1)
		with open(self.fq1, 'r') as fFQ1:
			for i, line in enumerate(fFQ1):
				if i % 4 == 0 :								# every four line represents one read #
					tmp_line = line.split(split_by)
					id = tmp_line[0]
					str_except_id = tmp_line[1]
#					if i % 10000 == 0 and i != 0:
#						sys.stdout.write("[FIXPE] %s: %d reads processed\n" %(self.fq1, i/4))
#						sys.stdout.flush()
				else:
					str_except_id += line
					if i % 4 == 3:
						fq1_dict[id] = str_except_id
						str_except_id = ""
		sys.stdout.write("\nBefore:\n")
		sys.stdout.write("%s: %d reads\n" %(self.fq1, (i+1)/4))
		sys.stdout.flush()

		fOFQ1 = open(self.outp + "_r1.preqc.fastq", 'w')
		fOFQ2 = open(self.outp + "_r2.preqc.fastq", 'w')
		fORPHAN = open(self.outp + "_orphan.preqc.fastq", 'w')

		split_by = self._get_splitby(self.fq1)
		num_paired, num_orphan = 0, 0
		with open(self.fq2, 'r') as fFQ2:
			for i, line in enumerate(fFQ2):
				if i % 4 == 0:
					tmp_line = line.split(split_by)
					id = tmp_line[0]
					str_except_id = tmp_line[1]
				else:
					str_except_id += line
					if i % 4 == 3:
						if id in fq1_dict:
							fOFQ1.write(id + " %s" %(fq1_dict[id]))
							fOFQ2.write(id + " %s" %(str_except_id))
							num_paired += 1
							fq1_dict.pop(id)
						else:
							fORPHAN.write(id + " %s" %(str_except_id))
							num_orphan += 1
						str_except_id = ""
		sys.stdout.write("%s: %d reads\n" %(self.fq2, (i+1)/4))
		for id, str_except_id in fq1_dict.iteritems():
			fORPHAN.write(id + " %s" %(str_except_id))
			num_orphan += 1
		fOFQ1.close()
		fOFQ2.close()
		fORPHAN.close()

		sys.stdout.write("After:\n")
		sys.stdout.write("%s: %d\n" %(self.outp + "_r1.preqc.fastq", num_paired))
		sys.stdout.write("%s: %d\n" %(self.outp + "_r2.preqc.fastq", num_paired))
		sys.stdout.write("%s: %d\n\n" %(self.outp + "_orphan.preqc.fastq", num_orphan))

class FixPE_Save(FixPE):
	def __init__(self, args):
		super(FixPE_Save, self).__init(args)

	def get_header_pattern(fq):
		cmd = "head -n1 %s" %(fq1)
		p = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
		return p.communicate()[0].strip().split(':')[0]

	def get_ids(basename, fq, header_pattern):
		id_file = basename + ".ids"
		grep_cmd = "grep ^%s %s" %(header_pattern, fq)
		awk_cmd = "awk '{print $1}' "
		sort_cmd = "sort "
		pgrep = subprocess.Popen(shlex.split(grep_cmd), stdout = subprocess.PIPE)
		pawk = subprocess.Popen(shlex.split(awk_cmd), stdin = pgrep.stdout, stdout = subprocess.PIPE)
		psort = subprocess.Popen(shlex.split(sort_cmd), stdin = pawk.stdout, stdout = open(id_file, 'w'))
		psort.wait()
		return id_file

	def get_union_ids(id_file_1, id_file_2):
		all_id_file = "all.ids"
		cat_cmd = "cat %s %s" %(id_file_1, id_file_2)
		sort_cmd = "sort "
		uniq_cmd = "uniq "
		pcat = subprocess.Popen(shlex.split(cat_cmd), stdout = subprocess.PIPE)
		psort = subprocess.Popen(shlex.split(sort_cmd), stdin = pcat.stdout, stdout = subprocess.PIPE)
		puniq = subprocess.Popen(shlex.split(uniq_cmd), stdin = psort.stdout, stdout = open(all_id_file, 'w'))
		puniq.wait()
		return all_id_file

	def get_indices(basename, fq, header_pattern):
		index_file = basename + ".index"
		if os.path.exists(index_file):
			print "skip indexing %s" %(index_file)
			return index_file
		print "get indices: %s" %(fq)
		grep_cmd = "grep -n ^%s %s" %(header_pattern, fq)
		awk_cmd = "awk '{split($0, a, \":\"); b = substr($0, index($0, \":\")+1); print b,a[1]}' "
		sort_cmd = "sort -n -k3 "
		pgrep = subprocess.Popen(shlex.split(grep_cmd), stdout = subprocess.PIPE)
		pawk = subprocess.Popen(shlex.split(awk_cmd), stdin = pgrep.stdout, stdout = subprocess.PIPE)
		psort = subprocess.Popen(shlex.split(sort_cmd), stdin = pawk.stdout, stdout = open(index_file, 'w'))
		psort.wait()
		return index_file

	def get_seq(pdict, fq, fP, fSE):
		missing = -1
		for line_num, line in enumerate(open(fq, 'r')):
			if line_num % 4 == 0:
				if line_num+1 in pdict:
					if pdict[line_num+1] == 1:
						fP.write(line)
						missing = 0
					else:
						fSE.write(line)
						missing = 1
			else:
				if not missing:
					fP.write(line)
				else:
					fSE.write(line)
			fP.flush()
			fSE.flush()

	def fix_pairing(index_file_1, index_file_2, fq1, fq2, paired_file_1, paired_file_2, se_file):
		print "fixing pairing"
		awk_cmd = """
				awk 'NR==FNR {a[$1] = $3; next;}
				{if($1 in a)
					{print $1, a[$1], $3; delete a[$1];}
				 else
					{print $1, \"-2\", $3;}}
				END {for (i in a) {print i, \"-1\", a[i];}}' %s %s
				""" %(index_file_1, index_file_2)
		pawk = subprocess.Popen(shlex.split(awk_cmd.strip()), stdout = subprocess.PIPE)
		fPAIR1 = open(paired_file_1, 'w')
		fPAIR2 = open(paired_file_2, 'w')
		fSE = open(se_file, 'w')

		pdict1, pdict2 = {}, {}
		sedict = {}
		for line in pawk.communicate()[0].splitlines():
			tmp_line = line.split(' ')
			if int(tmp_line[1]) > 0:
				pdict1[int(tmp_line[1])] = 1
				pdict2[int(tmp_line[2])] = 1
			else:
				if int(tmp_line[1]) == -1:
					pdict1[int(tmp_line[2])] = int(tmp_line[1])
				elif int(tmp_line[1]) == -2:
					pdict2[int(tmp_line[2])] = int(tmp_line[1])

		get_seq(pdict1, fq1, fPAIR1, fSE)
		pdict1.clear()
		get_seq(pdict2, fq2, fPAIR2, fSE)
		pdict2.clear()

		fPAIR1.close()
		fPAIR2.close()
		fSE.close()

#if __name__ == "__main__":
#	parser = ArgumentParser(description="Checking pair-end reads, only work for Casava 1.8+ style header")
#	parser.add_argument("-fq1", metavar="FILE", dest="fq1", required=True, help="First-end fastq file. Required")
#	parser.add_argument("-fq2", metavar="FILE", dest="fq2", required=True, help="Second-end fastq file. Required")
#	parser.add_argument("-p", metavar="STR", dest = "outprefix", required=True, help="output prefix")
#	parser.add_argument("-interleaved", action="store_true", dest="interleaved", help="toggle on to output interleaved fastq file. Optional.")
#
#	args = parser.parse_args()
#
#	check_inputs(args.fq1, args.fq2)
#
#	fix_pairing(args.fq1, args.fq2, args.outprefix, args.interleaved)
