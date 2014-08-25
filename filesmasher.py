import sys
import subprocess
import shlex
import os
import multiprocessing
import re

from io_simo import check_file_if_exists
from io_simo import make_dir_if_necessary
from ftsensor import FtSensor

class FileSmasher(object):
	def __init__(self, infile, num_chunk, outp, task_q):
		self.infile = os.path.realpath(infile)
		self.num_chunk = num_chunk
		self.outp = os.path.realpath(outp)
		self.task_q = task_q
		self._check_cmd()

	def _check_cmd(self):
		check_file_if_exists(self.infile)
		make_dir_if_necessary(os.path.dirname(self.outp))
		if self.num_chunk < 1:
			sys.stderr.write("[FileSmasher] Error: number of chunks must be larger than zero\n")
			sys.exit(1)
			sys.stderr.write("[FileSmasher] Error: only one chunk will be split to input file\n")
			sys.stderr.write("[FileSmasher] Warning: only one chunk will be split to input file\n")

	def _show(self):
		sys.stdout.write("[FileSmasher] input: %s\n" %(self.infile))
		sys.stdout.write("[FileSmasher] %d Bytes file smashes into %d chunks\n" %(self._get_fsize(), self.num_chunk))
		sys.stdout.write("[FileSmasher] ouput prefix: %s\n" %(self.outp))

	def _init_processes(self, nproc, out_suffix):
		for _ in range(nproc):
			p = multiprocessing.Process(target=self._smash, args=(out_suffix, ))
			p.daemon = True
			p.start()

	def _get_fsize(self):
		''' get input file size '''
		try:
			fsize = os.path.getsize(self.infile)
		except os.error as e:
			sys.stderr.write("[FileSmasher] Error: %s\n" %(e))
			sys.exit(1)
		else:
			return fsize

	def _get_chunk_size(self):
		return self._get_fsize()/self.num_chunk

	def _smash(self, out_suffix):
		while True:
			try:
				chunk_start, chunk_end, chunk_num = self.task_q.get()
				sys.stdout.write("[FileSmasher] %s reading %s from %s to %s\n" %(multiprocessing.current_process().name, self.infile, chunk_start, chunk_end))
				with open(self.infile, 'r') as fIN:
					fIN.seek(chunk_start)
					if chunk_end:
						seq = fIN.read(chunk_end)
					else:
						seq = fIN.read()
				outfile = self.outp + ".%d.%s" %(chunk_num, out_suffix)
				sys.stdout.write("[FileSmasher] %s pushing to %s ...\n" %(multiprocessing.current_process().name, outfile))
				fOUT = open(outfile, 'w')
				fOUT.write(seq)
				fOUT.flush()
				fOUT.close()
			finally:
				self.task_q.task_done()

class FqSmasher(FileSmasher):
	def __init__(self, infile, num_chunk, outp, task_q):
		super(FqSmasher, self).__init__(infile, num_chunk, outp, task_q)
		self._check_fmt(self.infile)
		self.delim = self._get_delim()
		self.out_suffix = "fastq"

	def _check_fmt(self, file):
		if FtSensor(file).getfiletype() != "fq":
			sys.stderr.write("[FqSmasher] Error: input file format error\n")
			sys.exit(1)

	def _get_delim(self):
		''' get split boundary symbol '''
		cmd = "head -n 1 %s" %(self.infile)
		p = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
		header = p.communicate()[0]
		delim = ""
		if len(header.split(' ')) == 2:		#@... 1:N:0:...\n
			delim = "[1|2]:N:0:"
		elif len(header.split('/')) == 2:		#@...#.../1[2]\n
			delim = "[\/1|\/2]"
		return delim

	def getchunk(self):
		''' get chunk boundaries '''
		fIN = open(self.infile, 'r')
		chunk_size = super(FqSmasher, self)._get_chunk_size()
		chunk_num = 0
		while True:
			start = fIN.tell()
			fIN.seek(chunk_size, 1)
			line = fIN.readline()
			if not line:
				break
			while not re.search(self.delim, line.strip()):
				line = fIN.readline()
			else:
				fIN.seek(-len(line), 1)
				chunk_num += 1
				self.task_q.put((start, fIN.tell()-start, chunk_num))
		self.task_q.put((start, None, chunk_num+1))
		fIN.close()

	def start(self, nproc):
		super(FqSmasher, self)._show()
		if nproc < 1:
			sys.stderr.write("[FqSmasher] Error: number of processors must be larger than zero\n")
			sys.exit(1)
		super(FqSmasher, self)._init_processes(nproc, self.out_suffix)
		self.getchunk()
		try:
			self.task_q.join()
		except KeyboardInterrupt:
			sys.stderr.write("[FqSmasher] Terminated unexpectedly by keyboard\n")
			sys.exit(0)

class FaSmasher(FileSmasher):
	def __init__(self, infile, num_chunk, outp, task_q):
		super(FaSmasher, self).__init__(infile, num_chunk, outp, task_q)
		self._check_fmt()
		self.delim = ">"
		self.out_suffix = "fasta"

	def _check_fmt(self):
		if FtSensor(self.infile).getfiletype() != "fa":
			sys.stderr.write("[FaSmasher] Error: input file format error\n")
			sys.exit()

	def getchunk(self):
		''' get chunk boundaries '''
		super(FaSmasher, self)._show()
		fIN = open(self.infile, 'r')
		chunk_size = super(FaSmasher, self)._get_chunk_size()
		chunk_num = 0
		while True:
			start = fIN.tell()
			fIN.seek(chunk_size, 1)
			line = fIN.readline()
			if not line:
				break
			while not line.startswith(self.delim):
				line = fIN.readline()
			else:
				fIN.seek(-len(line), 1)
				chunk_num += 1
				self.task_q.put((start, fIN.tell()-start, chunk_num))
		self.task_q.put((start, None, chunk_num+1))
		fIN.close()
		return

	def start(self, nproc):
		self._show()
		if nproc < 1:
			sys.stderr.write("[FaSmasher] Error: number of processors must be larger than zero\n")
			sys.exit(1)
		super(FaSmasher, self)._init_processes(nproc, self.out_suffix)
		self.getchunk()
		try:
			self.task_q.join()
		except KeyboardInterrupt:
			sys.stderr.write("[FaSmasher] Terminated unexpectedly by keyboard\n")
			sys.exit()
