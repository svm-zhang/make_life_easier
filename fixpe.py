#!/N/soft/mason/python/2.7.3/bin/python
# Simo Zhang

import re
import os
import sys
import subprocess
import shlex

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
		hfmt = FtSensor(file).gethfmt()
		if hfmt == "illumina14":
			return "/"
		elif hfmt == "casava18":
			return " "
		elif hfmt == "mix":
			sys.stderr.write("[FMT CHECKER] %s includes at least two different header formats\n")
			sys.exit(1)

	def start(self):
		self._fixpe()

	def _fixpe(self):
		''' fix reads pairing '''
		id, str_except_id = "", ""
		fq1_dict = {}
		split_by1 = self._get_splitby(self.fq1)
		with open(self.fq1, 'r') as fFQ1:
			for i, line in enumerate(fFQ1):
				if i % 4 == 0 :								# every four line represents one read #
					tmp_line = line.split(split_by1)
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
		nreads1 = (i+1)/4

		fOFQ1 = open(self.outp + "_r1.preqc.fastq", 'w')
		fOFQ2 = open(self.outp + "_r2.preqc.fastq", 'w')
		fORPHAN = open(self.outp + "_orphan.preqc.fastq", 'w')

		split_by2 = self._get_splitby(self.fq2)
		num_paired, num_orphan = 0, 0
		with open(self.fq2, 'r') as fFQ2:
			for i, line in enumerate(fFQ2):
				if i % 4 == 0:
					tmp_line = line.split(split_by2)
					id = tmp_line[0]
					str_except_id = tmp_line[1]
				else:
					str_except_id += line
					if i % 4 == 3:
						if id in fq1_dict:
							fOFQ1.write(id + "%s%s" %(split_by1, fq1_dict[id]))
							fOFQ2.write(id + "%s%s" %(split_by2, str_except_id))
							num_paired += 1
							fq1_dict.pop(id)
						else:
							fORPHAN.write(id + "%s%s" %(split_by2, str_except_id))
							num_orphan += 1
						str_except_id = ""
		for id, str_except_id in fq1_dict.iteritems():
			fORPHAN.write(id + "%s%s" %(split_by1, str_except_id))
			num_orphan += 1
		fOFQ1.close()
		fOFQ2.close()
		fORPHAN.close()

		sys.stdout.write("\nBefore:\n")
		sys.stdout.write("%s: %d reads\n" %(self.fq1, nreads1))
		sys.stdout.write("%s: %d reads\n" %(self.fq2, (i+1)/4))
		sys.stdout.write("After:\n")
		sys.stdout.write("%s: %d\n" %(self.outp + "_r1.preqc.fastq", num_paired))
		sys.stdout.write("%s: %d\n" %(self.outp + "_r2.preqc.fastq", num_paired))
		sys.stdout.write("%s: %d\n\n" %(self.outp + "_orphan.preqc.fastq", num_orphan))
		sys.stdout.flush()

class FixPE_Save(FixPE):
	def __init__(self, args):
		sys.stdout.write("[FIXPE] Using the memory saving mode\n")
		super(FixPE_Save, self).__init__(args)

	def _get_indices(self, file, index):
		if os.path.exists(index):
			sys.stdout.write("skip indexing %s\n" %(index))
			return
		sys.stdout.write("get indices: %s\n" %(file))
		hfmt = FtSensor(file).gethfmt()
#		grep_cmd = "grep -n ^%s %s" %(header_pattern, fq)
#		awk_cmd = "awk '{split($0, a, \":\"); b = substr($0, index($0, \":\")+1); print b,a[1]}' "
#		pgrep = subprocess.Popen(shlex.split(grep_cmd), stdout = subprocess.PIPE)
#		sort_cmd = "sort -n -k2 "
		awk_cmd = "awk \'{if(NR %% 4 == 1) {split($0, a, \"/\"); print NR,a[1];}}\' %s" %(file)
		pawk = subprocess.Popen(shlex.split(awk_cmd), stdout = subprocess.PIPE)
#		psort = subprocess.Popen(shlex.split(sort_cmd), stdin = pawk.stdout, stdout = subprocess.PIPE)
#		psort_out = psort.communicate()[0]
		pawk_out = pawk.communicate()[0]
		fINDEX = open(index, 'w')
		fINDEX.write(pawk_out)
		fINDEX.close()
		return

	def _get_seq(self, pdict, fq, fP, fORPHAN):
		missing = -1
		for line_num, line in enumerate(open(fq, 'r')):
			if line_num % 4 == 0:
				if line_num+1 in pdict:
					if pdict[line_num+1] == 1:
						fP.write(line)
						missing = 0
						del pdict[line_num+1]
					else:
						fORPHAN.write(line)
						missing = 1
			else:
				if not missing:
					fP.write(line)
				else:
					fORPHAN.write(line)
			fP.flush()
			fORPHAN.flush()
#		if os.getpid() > 0:
#			os._exit(0)

	def start(self):
		self._fixpe_save()

	def _fixpe_save(self):
		''' fix reads pairng in memmory saving mode '''
		index1 = self.outp + "_r1.index"
		index2 = self.outp + "_r2.index"
		self._get_indices(self.fq1, index1)
		self._get_indices(self.fq2, index2)

		sys.stdout.write("[FIXPE] Matching Reads ID ...\n")
		awk_cmd = """
				awk 'NR==FNR {a[$2] = $1; next;}
				{if($2 in a)
					{print $2, a[$2], $1; delete a[$2];}
				 else
					{print $2, \"-2\", $1;}}
				END {for (i in a) {print i, \"-1\", a[i];}}' %s %s
				""" %(index1, index2)
		pawk = subprocess.Popen(shlex.split(awk_cmd.strip()), stdout = subprocess.PIPE)

		sys.stdout.write("[FIXPE] Making pairs and rescuing orphans ...\n")
		pdict1, pdict2 = {}, {}
		sedict = {}
		for line in pawk.communicate()[0].splitlines():
			tmp_line = line.split(' ')
			if int(tmp_line[1]) > 0:
				pdict1[int(tmp_line[1])] = 1
				pdict2[int(tmp_line[2])] = 1
			else:
				if int(tmp_line[1]) == -1:
					pdict1[int(tmp_line[2])] = -1
				elif int(tmp_line[1]) == -2:
					pdict2[int(tmp_line[2])] = -2

		sys.stdout.write("[FIXPE] Pushing to files ...\n")
		fOFQ1 = open(self.outp + "_r1.preqc.fastq", 'w')
		fOFQ2 = open(self.outp + "_r2.preqc.fastq", 'w')
		fORPHAN1 = open(self.outp + "_orphan.preqc.1.fastq", 'w')
		fORPHAN2 = open(self.outp + "_orphan.preqc.2.fastq", 'w')
		pid = os.fork()
		if pid > 0:
			sys.stdout.write("[FIXPE] Pushing to %s ...\n" %(self.outp+"_r2.preqc.fastq"))
			self._get_seq(pdict2, self.fq2, fOFQ2, fORPHAN2)
			pdict2.clear()
		else:
			child = pid
			sys.stdout.write("[FIXPE] Pushing to %s ...\n" %(self.outp+"_r1.preqc.fastq"))
			self._get_seq(pdict1, self.fq1, fOFQ1, fORPHAN1)
			pdict1.clear()
			os._exit(0)
		fOFQ1.close()
		fOFQ2.close()
		fORPHAN1.close()
		fORPHAN2.close()

		sys.stdout.write("[FIXPE] Concatenating orphan preqc fastq files ...\n")
		cat_cmd = "cat %s %s " %(self.outp+"_orphan.preqc.1.fastq", self.outp+"_orphan.preqc.2.fastq")
		pcat = subprocess.Popen(shlex.split(cat_cmd), stdout = open(self.outp + "_orphan.preqc.fastq", 'w'))
		pcat_out, pcat_err = pcat.communicate()
		if pcat.returncode != 0:
			sys.stderr.write("[FIXPE] Fail to concatenate 2 orphan fastq files\n")
			sys.exit(1)

		sys.stdout.write("[FIXPE] Cleaning ...\n")
		os.remove(self.outp+"_orphan.preqc.1.fastq")
		os.remove(self.outp+"_orphan.preqc.2.fastq")
