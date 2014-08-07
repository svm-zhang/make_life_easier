import sys
import subprocess
import shlex
import os

class FileSmasher(object):
	def __init__(self, infile, num_chunk, task_q):
		self.infile = infile
		self.num_chunk = num_chunk
		self.task_q = task_q

	def show(self):
		sys.stdout.write("input: %s\n" %(self.infile))
		sys.stdout.write("%d Bytes file smashes into %d chunks\n" %(self._get_fsize(), self.num_chunk))

	def _get_fsize(self):
		try:
			fsize = os.path.getsize(self.infile)
		except os.error as e:
			sys.stderr.write("%s\n" %(e))
			sys.exit(1)
		else:
			return fsize

	def get_chunk_size(self):
		return self._get_fsize()/self.num_chunk

class FqSmasher(FileSmasher):
	def __init__(self, infile, num_chunk, task_q):
		super(FqSmasher, self).__init__(infile, num_chunk, task_q)
		self.delim = self._get_delim()

	def _get_delim(self):
		cmd = "head -n 1 %s" %(self.infile)
		p = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
		header = p.communicate()[0]
		if len(header.split(' ')) == 2:		#@... 1:N:0:...\n
			delim = header.strip().split(' ')[1]
		elif len(header.split('/')) == 2:		#@...#.../1[2]\n
			delim = ("/1", "/2")
		return delim

	def show(self):
		super(FqSmasher, self).show()
		sys.stdout.write("delim: %s\n" %(self.delim, ))

	def start(self):
		self._fqsmasher()
		return

	def _fqsmasher(self):
		''' split FASTQ into chunks/blocks '''
		self.show()
		fIN = open(self.infile, 'r')
		chunk_size = super(FqSmasher, self).get_chunk_size()
		chunk_num = 0
		while True:
			start = fIN.tell()
			fIN.seek(chunk_size, 1)
			line = fIN.readline().strip()
			if not line:
				break
			while not line.endswith(self.delim):
				line = fIN.readline().strip()
			else:
				fIN.seek(-len(line), 1)
				chunk_num += 1
				self.task_q.put((start, fIN.tell()-start, chunk_num))
		self.task_q.put((start, None, chunk_num+1))
		fIN.close()
		return

class FaSmasher(FileSmasher):
	def __init__(self, infile, num_chunk, task_q):
		super(FaSmasher, self).__init__(infile, num_chunk, task_q)
		self.delim = ">"
		self.show()

	def show(self):
		super(FaSmasher, self).show()
		sys.stdout.write("delim: %s\n" %(self.delim, ))

	def start(self):
		self._fasmasher()
		return

	def _fasmasher(self):
		''' split FASTQ into chunks/blocks '''
		self.show()
		fIN = open(self.infile, 'r')
		chunk_size = super(FaSmasher, self).get_chunk_size()
		chunk_num = 0
		while True:
			start = fIN.tell()
			fIN.seek(chunk_size, 1)
			line = fIN.readline().strip()
			if not line:
				break
			while not line.startswith(self.delim):
				line = fIN.readline().strip()
			else:
				fIN.seek(-len(line), 1)
				chunk_num += 1
				self.task_q.put((start, fIN.tell()-start, chunk_num))
		self.task_q.put((start, None, chunk_num+1))
		fIN.close()
		return

#fq = sys.argv[1]
#fa = sys.argv[2]
#
#task_q = multiprocessing.JoinableQueue()
#FqSmasher(fq, 4, task_q).start()
#FaSmasher(fa, 4, task_q).start()
