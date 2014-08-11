import sys
import os
import multiprocessing

import io_simo
from ftsensor import FtSensor

class Fx2Fy(object):
	'''
		from FileA to FileB
	'''
	def __init__(self, file_a, file_b):
		self.file_a = file_a
		self.file_b = os.path.realpath(file_b)
		io_simo.check_file_if_exists(file_a)
		io_simo.make_dir_if_necessary(os.path.dirname(self.file_b))

class Fq2Fa(Fx2Fy):
	''' from FastQ to FastA '''
	def __init__(self, file_a, file_b):
		super(Fq2Fa, self).__init__(file_a, file_b)

	def start(self):
		# make sure the input file is FastQ
		if FtSensor(self.file_a).getfiletype() != "fq":
			sys.stderr.write("[FOMAT Error] %s is not in FastQ format\n" %(self.file_a))
			sys.exit()
		else:
			sys.stdout.write("[Fq2Fa] converting %s to %s\n" %(self.file_a, self.file_b))
			return self._run()

	def _run(self):
		num_seq = 0
		fFA = open(self.file_b, 'w')
		with open(self.file_a, 'r') as fFQ:
			for i, line in enumerate(fFQ):
				if i % 4 == 0:
					fFA.write(">" + line[1:])
				elif i % 4 == 1:
					fFA.write(line)
				else:
					num_seq += 1
				fFA.flush()
		return num_seq

class Fqs2Fas(Fq2Fa):
	''' from FastQ to FastA mp mode '''
	def __init__(self, files, nproc):
		self.files = files
		self.nproc = nproc

	def start(self):
		task_q = multiprocessing.JoinableQueue()
		for _ in range(self.nproc):
			p = multiprocessing.Process(target=self._run, args=(task_q, ))
			p.daemon = True
			p.start()
		i = 0
		while i < len(self.files):
			fq = os.path.realpath(self.files[i])
			fa = os.path.realpath(self.files[i+1])
			io_simo.check_file_if_exists(fq)
			io_simo.make_dir_if_necessary(os.path.dirname(fa))
			if FtSensor(fq).getfiletype() != "fq":
				sys.stderr.write("[FOMAT Error] %s is not in FastQ format\n" %(fq))
				sys.exit()
			task_q.put((fq, fa))
			i += 2
		try:
			task_q.join()
		except KeyboardInterrupt:
			sys.stderr.write("Terminated unexpectedly by keyboard\n")
			sys.exit()

	def _run(self, task_q):
		while True:
			try:
				fq, fa = task_q.get()
				sys.stdout.write("%s\tconverting %s\t%s\n" %(multiprocessing.current_process().name, fq, fa))
				num_seq = 0
				fFA = open(fa, 'w')
				with open(fq, 'r') as fFQ:
					for i, line in enumerate(fFQ):
						if i % 4 == 0:
							fFA.write(">" + line[1:])
						elif i % 4 == 1:
							fFA.write(line)
						else:
							num_seq += 1
						fFA.flush()
			finally:
				task_q.task_done()

class Fa2Fq(Fx2Fy):
	'''
		from FastA to FastQ
	'''
	def __init__(self, file_a, file_b, qual_file = ""):
		super(Fa2Fq, self).__init__(file_a, file_b)
		self.qual_file = qual_file
