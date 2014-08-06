#!/N/soft/mason/python/2.7.3/bin/python

import re
import os
import sys
import multiprocessing
import argparse
import glob
import subprocess
import shlex

def allocate_tasks(fastqs, task_queue):
	''' allocate tasks for each initiated process '''
	sys.stdout.write("[Progress] Allocating tasks ...\n")
	for fastq in fastqs:
		task_queue.put(fastq)
	return

def create_processes(nproc, task_queue, target, nthreads, noextract, kmer, fastqc_prg):
	''' initiate processes '''
	sys.stdout.write("[Progress] Initiate Processes ...\n")
	for _ in range(nproc/nthreads):
		process = multiprocessing.Process(target=run_FastQC, args=(task_queue, target, nthreads, noextract, kmer, fastqc_prg))
		process.daemon = True
		process.start()
	return

def run_FastQC(task_queue, target, nthreads, noextract, kmer, fastqc_prg):
	''' run FastQC on each given fastq file '''
	while True:
		try:
			fastq = task_queue.get()
			if noextract:
				run_cmd = "%s -o %s --noextract -t %d -k %d -q %s" %(fastqc_prg, target, nthreads, kmer, fastq)
			else:
				run_cmd = "%s -o %s -t %d -k %d -q %s" %(fastqc_prg, target, nthreads, kmer, fastq)
			sys.stdout.write("[Progress] %s\trunning FastQC on %s\n" %(multiprocessing.current_process().name, fastq))
			try:
				p = subprocess.Popen(shlex.split(run_cmd), stdout = subprocess.PIPE)
				pout, perr = p.communicate()
				if p.returncode != 0:
					raise RuntimeError("[Error] %r failed, status code %s\n[Erorr] STDOUT: %r\n[Error] STDOUT: %r\n" %(cmd, p.returncode, pout, perr))
			except (OSError, ValueError), e:
				sys.stderr.write("[Erorr] %s\n" %(e))
				sys.exit()
			else:
				sys.stdout.write("[Progress] %s\tFastQC successfully finishes on %s\n" %(multiprocessing.current_process().name, fastq))
		finally:
			task_queue.task_done()
	for i in range(len(files_to_qc)):
		fastqc_cmd = "fastqc -o %s --noextract -t 1 %s" %(out_dir, files_to_qc[i])
		print fastqc_cmd
		system(fastqc_cmd)

def handle_cmd():
	''' parse what passed in from command line '''
	parser = argparse.ArgumentParser(description="running fastqc to given fastq files")
	parser.add_argument("-target", metavar="DIR", dest="target", required=True, help="target directory where results will be put")
	parser.add_argument("-nproc", metavar="INT", type=int, dest="nproc", default=1, help="number of cores to request [1]")
	fastqc_group = parser.add_argument_group(description="Arguments for running FastQC")
	fastqc_group.add_argument("-nthreads", metavar="INT", type=int, dest="nthreads", default=1, help="number of threads for each fastqc process [1]")
	fastqc_group.add_argument("-noextract", dest="noextract", action="store_true", help="specify to ask FastQC not to uncompress the output files [false/uncompress]")
	fastqc_group.add_argument("-kmer", metavar="INT", type=int, dest="kmer", default=5, choices=range(2,11), help="specify the length of Kmer to look for in FastQC [5]")
	fastqc_group.add_argument("-fastqc", metavar="FILE", dest="fastqc_prg", default="fastqc", help="specify a path to where FastQC is installed. If not specified, \"fastqc\" will be used [fastqc]")
	mutual_group = parser.add_mutually_exclusive_group(required=True)
	mutual_group.add_argument("-source", metavar="DIR", dest="source", help="source directory where to-be-processed-fastq files are stored. Cannot use with -list at the same time")
	mutual_group.add_argument("-list", metavar="FILE", dest="source", help="a file with a list of fastq files to be processed")

	args = parser.parse_args()
	return os.path.realpath(args.target), args.nproc, os.path.realpath(args.source), args.nthreads, args.noextract, args.kmer, args.fastqc_prg

def main():
	target, nproc, source, nthreads, noextract, kmer, fastqc_prg = handle_cmd()

	if not os.path.exists(target):
		os.makedirs(target)

	if not os.path.exists(source):
		sys.stderr.write("[Error] cannot find %s\n" %(source))
		sys.exit()

	# positive number for nproc and nthreads arguments
	if nproc < 1 or nthreads < 1:
		sys.stderr.write("[Error] integers larger than 1 are required for nproc and nthreads\n")
		sys.exit()

	# nproc needs to be larger than nthreads
	if nthreads > nproc:
		sys.stderr.write("[Error] nproc is required to be larger than nthreads\n")
		sys.exit()

	if nproc % nthreads != 0:
		sys.stdout.write("[Warning] it is highly recommended nproc can be divided by nthreads completely\n")

	fastqc_cmd = "%s -h" %(fastqc_prg)
	try:
		fastqc_attempt = subprocess.Popen(shlex.split(fastqc_cmd), stdout = subprocess.PIPE)
	except OSError, e:
		sys.stderr.write("[Error] Cannot find the executable of FastQC as provided: %s\n" %(fastqc_prg))
		sys.exit()

	fastqs = []
	if os.path.isdir(source):
		for file in glob.glob(os.path.join(source, "*.f*q")):
			fastqs.append(file)
	elif os.path.isfile(source):
		with open(source, 'r') as fSOURCE:
			for line in fSOURCE:
				if not line.startswith('#'):
					if not os.path.exists(line.strip()):
						sys.stderr.write("[Error] cannot find %s\n" %(line.strip()))
						sys.exit()
					fastqs.append(line.strip())
	sys.stdout.write("[Progress] %d fastq files parsed\n" %(len(fastqs)))

	task_queue = multiprocessing.JoinableQueue()
	create_processes(nproc, task_queue, target, nthreads, noextract, kmer, fastqc_prg)
	allocate_tasks(fastqs, task_queue)
	try:
		task_queue.join()
	except KeyboardInterrupt:
		sys.stderr.write("Terminated unexpectedly by keyboard\n")
		sys.exit()
	else:
		sys.stderr.write("Program finishes normally\n")

if __name__ == "__main__":
	main()
