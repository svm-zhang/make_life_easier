import os
import sys
import subprocess
import shlex

class FtSensor(object):
	def __init__(self, infile):
		self.infile = infile

	def getfiletype(self):
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

	def _get_total_lines(self):
		''' get total number of lines in a file'''

	def _getlines(self, num_line):
		''' get given number of lines in a file '''
		cmd = "head -n %d %s" %(num_line, self.infile)
		try:
			p = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
			pout, perr = p.communicate()
		except (OSError, ValueError) as e:
			sys.stderr.write("%s\n" %(perr))
			sys.stderr.write("%s\n" %(e))
		return pout

	def gethfmt(self):
		''' get the header line format '''
		hfmt = None
		for i, line in enumerate(self._getlines(100000).splitlines()):
			if i % 4 == 0:
				if hfmt == "mix":
					break
				if len(line.strip().split(' ')) == 2:
					if hfmt == None:
						hfmt = "casava18"
					else:
						if hfmt != "casava18":
							hfmt = "mix"		# different header fmt detected
				elif len(line.strip().split('/')) == 2:
					if hfmt == None:
						hfmt = "illumina14"
					else:
						if hfmt != "illumina14":
							hfmt = "mix"		# different header fmt detected
		return hfmt

	def getencoding(self):
		''' sense the header format and quality encodings '''
		encoding = None
		quals = []
		for i, line in enumerate(self._getlines(100000).splitlines()):
			if i % 4 == 3:
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
		ft = self.getfiletype()
		if ft == "fa":
			return ft, None, None
		elif ft == "fq":
			hfmt, encoding = self.getencoding()
			return ft, hfmt, encoding
