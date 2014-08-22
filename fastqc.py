#!/N/soft/mason/python/2.7.3/bin/python

import re
import os
import sys
import multiprocessing
import argparse
import glob
import subprocess
import shlex

from io_simo import make_dir_if_necessary
from io_simo import check_file_if_exists

class FastQC:
	def __init__(self, args):
		self.target = os.path.realpath(args.target)
		self.source = os.path.realpath(args.source)
		self.nproc = args.nproc
		self.nthreads = args.nthreads
		self.noextract = args.noextract
		self.kmer = args.kmer
		self.fastqc_prg = args.fastqc_prg
		self._check_args()
		self._print_cmd()

	def _check_args(self):
		check_file_if_exists(self.source)
		make_dir_if_necessary(self.target)
		if self.nproc < 1 or self.nthreads < 1:
			sys.stderr.write("[mleFQ fastqc] Error: integers larger than 1 are required for nproc and nthreads\n")
			sys.exit(1)

		# check FastQC executables
		fastqc_cmd = "%s -h" %(self.fastqc_prg)
		try:
			fastqc_attempt = subprocess.Popen(shlex.split(fastqc_cmd), stdout = subprocess.PIPE)
		except OSError, e:
			sys.stderr.write("[mleFQ fastqc] Error: cannot find the executable of FastQC as provided: %s\n" %(self.fastqc_prg))
			sys.exit(1)

	def _print_cmd(self):
		sys.stdout.write("[mleFQ fastqc] target: %s\n" %(self.target))
		sys.stdout.write("[mleFQ fastqc] source: %s\n" %(self.source))
		sys.stdout.write("[mleFQ fastqc] nproc: %d\n" %(self.nproc))
		sys.stdout.write("[mleFQ fastqc] nthreads: %d\n" %(self.nthreads))
		sys.stdout.write("[mleFQ fastqc] noextract: %s\n" %(self.noextract))
		sys.stdout.write("[mleFQ fastqc] kmer: %d\n\n" %(self.kmer))

	def _getFastqs(self):
		fastqs = []
		if os.path.isdir(self.source):
			for file in glob.glob(os.path.join(self.source, "*.f*q")):
				fastqs.append(file)
		elif os.path.isfile(self.source):
			with open(self.source, 'r') as fSOURCE:
				for line in fSOURCE:
					if not line.startswith('#'):
						if not os.path.exists(line.strip()):
							sys.stderr.write("[Error] cannot find %s\n" %(line.strip()))
							sys.exit()
						fastqs.append(line.strip())
		sys.stdout.write("[mleFQ fastqc] %d fastq files parsed\n" %(len(fastqs)))
		return fastqs

	def _initProcesses(self, task_q):
		''' initiate processes '''
		sys.stdout.write("[mleFQ fastqc] Initiate Processes ...\n")
		for _ in range(self.nproc):
			process = multiprocessing.Process(target=self._runFastQC, args=(task_q, ))
			process.daemon = True
			process.start()
		return

	def start(self):
		self._run()

	def _run(self):
		# get Fastq files
		fastqs = self._getFastqs()

		# initiate processes
		task_q = multiprocessing.JoinableQueue()
		self._initProcesses(task_q)

		# allocate tasks for each initiated process
		sys.stdout.write("[mleFQ fastqc] Allocating tasks ...\n")
		for fastq in fastqs:
			task_q.put(fastq)

		# wait for the task_q empty, unless keyboard intervention
		try:
			task_q.join()
		except KeyboardInterrupt:
			sys.stderr.write("[mleFQ fastqc] Terminated unexpectedly by keyboard\n")
			sys.exit()
		else:
			sys.stderr.write("[mleFQ fastqc] Program finishes normally\n")

	def _runFastQC(self, task_q):
		''' run FastQC on each given fastq file '''
		while True:
			try:
				fastq = task_q.get()
				if self.noextract:
					run_cmd = "%s -o %s --noextract -t %d -k %d -q %s" %(self.fastqc_prg, self.target, self.nthreads, self.kmer, fastq)
				else:
					run_cmd = "%s -o %s -t %d -k %d -q %s" %(self.fastqc_prg, self.target, self.nthreads, self.kmer, fastq)
				sys.stdout.write("[mleFQ fastqc] %s\trunning FastQC on %s\n" %(multiprocessing.current_process().name, fastq))
				try:
					p = subprocess.Popen(shlex.split(run_cmd), stdout = subprocess.PIPE)
					pout, perr = p.communicate()
					if p.returncode != 0:
						raise RuntimeError("[Error] %r failed, status code %s\n[Erorr] STDOUT: %r\n[Error] STDOUT: %r\n" %(run_cmd, p.returncode, pout, perr))
				except (OSError, ValueError), e:
					sys.stderr.write("[Erorr] %s\n" %(e))
					sys.exit()
				else:
					sys.stdout.write("[mleFQ fastqc] %s\tFastQC successfully finishes on %s\n" %(multiprocessing.current_process().name, fastq))
			finally:
				task_q.task_done()
