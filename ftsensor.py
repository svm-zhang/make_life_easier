import os
import sys
import subprocess
import shlex

class FtSensor(object):
	def __init__(self, infile):
		self.infile = infile

	def _getfiletype(self):
		''' detect format of input file '''
		for i, line in enumerate(self._getlines(4).splitlines()):
			if i == 0:
				if line.startswith('>'):
					return "fa"
				elif line.startswith('@'):
					continue
			elif i == 2:
				if line.startswith('+'):
					return "fq"
		return None

	def _getlines(self, num_line):
		''' get given number of lines in a file '''
		cmd = "head -n %d %s" %(num_line, self.infile)
		p = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
		pout = p.communicate()[0]
		return pout

	def _getencoding(self):
		''' sense the header format and quality encodings '''
		ft, hfmt, encoding = None, None, None
		quals = []
		for i, line in enumerate(self._getlines(100000).splitlines()):
			if i % 4 == 0:
				if hfmt == "mix":
					continue
				if len(line.strip().split(' ')) == 2:
					if hfmt == "":
						hfmt = "illumina14"
					else:
						if hfmt != "illumina14":
							hfmt = "mix"		# different header fmt detected
				elif len(line.strip().split('/')) == 2:
					if hfmt == "":
						hfmt = "casava18"
					else:
						if hfmt != "casava18":
							hfmt = "mix"		# different header fmt detected
			elif i % 4 == 3:
				for qual in line.strip():
					quals.append(ord(qual))

		# get the range of quality encodings
		if min(quals) >= 33 and max(quals) <= 73:
			encoding = "sanger"
		elif min(quals) >= 64 and max(quals) <= 104:
			encoding = "illumina13+"
		elif min(quals) >= 66 and max(quals) <= 105:
			encoding = "illumina15+"
		elif min(quals) >= 33 and max(quals) <= 74:
			encoding = "illumina18+"

		return hfmt, encoding

	def bloodhound(self):
		'''
			begin sensing process
			sensing like a bloodhound
		'''
		ft = self._getfiletype()
		if ft == "fa":
			return ft, None, None
		elif ft == "fq":
			hfmt, encoding = self._getencoding()
			return ft, hfmt, encoding

#fq = sys.argv[1]
#
#print FtSensor(fq).bloodhound()
