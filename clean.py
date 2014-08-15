import os
import sys
import multiprocessing
import subprocess
import shlex

from ftsensor import FtSensor
from io_simo import check_files_if_exist
from io_simo import make_dir_if_necessary

class Clean(object):
	def __init__(self, inputs, nproc, task_q):
		self.inputs = inputs
		self.nproc = nproc
		self.task_q = task_q

	def _initProcesses(self):
		for _ in range(self.nproc):
			p = multiprocessing.Process(target=self._clean)
			p.daemon = True
			p.start()

	def _allocate(self):
		for input in self.inputs:
			fq, blast, out = input.split(',')
			check_files_if_exist(fq, blast)
			make_dir_if_necessary(os.path.dirname(out))
			self.task_q.put((fq, blast, out))

	def _readBLAST(self, blast):
		sys.stdout.write("%s\tgetting BLAST hits from %s\n" %(multiprocessing.current_process().name, blast))

		awk_cmd = "awk \'{print $1}\' %s" %(blast)
		uniq_cmd = "uniq "
		pawk = subprocess.Popen(shlex.split(awk_cmd), stdout = subprocess.PIPE)
		puniq = subprocess.Popen(shlex.split(uniq_cmd), stdin = pawk.stdout, stdout = subprocess.PIPE)

		hits = {}
		for hit in puniq.communicate()[0].splitlines():
			hits[hit] = 1

		return hits


	def _clean(self):
		while True:
			try:
				fq, blast, out = self.task_q.get()

				blast_hits = self._readBLAST(blast)
				sys.stdout.write("%s\t%d\n" %(multiprocessing.current_process().name, len(blast_hits)))
				hfmt, _ = FtSensor(fq).getencoding()
				fOUT = open(out, 'w')
				flag = 0
				i = 0
				sys.stdout.write("%s\tfiltering reads from %s\n" %(multiprocessing.current_process().name, fq))
				with open(fq, 'r') as fFQ:
					for i, line in enumerate(fFQ):
						if i % 4 == 0:
							if hfmt == "illumina14":
								tmp_line = line.strip().split()
							elif hfmt == "casava18":
								tmp_line = line.strip().split(" ")
							if not tmp_line[0].lstrip('@') in blast_hits:
								flag = 1
								fOUT.write(line)
						else:
							if flag == 1:
								fOUT.write(line)
								if i % 4 == 3:
									fOUT.flush()
									flag = 0
				fOUT.close()
			finally:
				self.task_q.task_done()

	def start(self):
		self._run()

	def _run(self):
		self._initProcesses()
		self._allocate()
	
		try:
			self.task_q.join()
		except KeyboardInterrupt:
			sys.stderr.write("Terminated unexpectedly by keyboard\n")
			sys.exit()
