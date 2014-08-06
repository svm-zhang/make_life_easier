#!/N/soft/mason/python/2.7.3/bin/python

import os
import sys
import multiprocessing
import argparse
import subprocess
import shlex
import tempfile

class Sanitizer:
	def __init__(self, args):
		self.infile = os.path.realpath(args.infile)
		self.outp = os.path.realpath(args.outp)
		self.num_chunk = args.num_chunk
		self.db = os.path.realpath(args.db)
		self.nthreads = args.nthreads
		self.evalue = args.evalue
		self.print_cmd()

	def print_cmd(self):
		sys.stdout.write("[CMD PARSER] infile: %s\n" %(os.path.realpath(self.infile)))
		sys.stdout.write("[CMD PARSER] output prefix: %s\n" %(os.path.realpath(self.outp)))
		sys.stdout.write("[CMD PARSER] num chunks: %d\n\n" %(self.num_chunk))

		sys.stdout.write("[CMD PARSER] BLAST database: %s\n" %(os.path.realpath(self.db)))
		sys.stdout.write("[CMD PARSER] num threads: %d\n" %(self.nthreads))
		sys.stdout.write("[CMD PARSER] min E-Value: %f\n\n" %(self.evalue))

	def sanitize(self, task_q, blast_q, err_q, fmt):
		while True:
			try:
				chunk_start, chunk_end, chunk_num = task_q.get()
				sys.stdout.write("%s\t%s\t%s\t%s\t%d\n" %(multiprocessing.current_process().name, chunk_start, chunk_end, self.infile, chunk_num))
				with open(self.infile, 'r') as fIN:
					fIN.seek(chunk_start)
					if chunk_end:
						seq = fIN.read(chunk_end)
					else:
						seq = fIN.read()
				out_blast = self.outp + ".%d.blast" %(chunk_num)
				blast_cmd = "blastn -db %s -outfmt 6 -out %s -num_threads %d -evalue %d " %(self.db, out_blast, self.nthreads, self.evalue)
				if fmt == "fa":
					sys.stdout.write(multiprocessing.current_process().name + "\t" + blast_cmd + "\n")
					p = subprocess.Popen(shlex.split(blast_cmd), stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
					bout, berr = p.communicate(seq)
				elif fmt == "fq":
					fa_file = tempfile.mkstemp(suffix=".fasta", dir=os.path.dirname(self.outp), prefix=os.path.basename(self.outp))[1]
					self.fq2fa(seq, fa_file)
					sys.stdout.write(multiprocessing.current_process().name + "\t" + blast_cmd + "\n")
					p = subprocess.Popen(shlex.split(blast_cmd), stdin = open(fa_file, 'r'), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
					bout, berr = p.communicate()
					os.remove(fa_file)
				if p.returncode == 1:
					err_q.put((multiprocessing.current_process().name, out_blast, p.returncode, "Error in query sequence(s) or BLAST options"))
				elif p.returncode == 2:
					err_q.put((multiprocessing.current_process().name, out_blast, p.returncode, "Error in Error in BLAST database"))
				elif p.returncode == 3:
					err_q.put((multiprocessing.current_process().name, out_blast, p.returncode, "Error in BLAST engine"))
				elif p.returncode == 4:
					err_q.put((multiprocessing.current_process().name, out_blast, p.returncode, "Out of memory"))
				elif p.returncode == 255:
					err_q.put((multiprocessing.current_process().name, out_blast, p.returncode, "Unknown error"))
				else:
					blast_q.put(out_blast)
			finally:
				task_q.task_done()

	def chunkize_fasta(self, infile_size, fmt, task_q):
		''' split FASTA into chunks/blocks '''
		delim = ""
		if fmt == "fa":
			delim = '>'
		elif fmt == "fq":
			cmd = "head -n 1 %s" %(self.infile)
			p = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
			delim = p.communicate()[0].split(':')[0]
		chunk_size = infile_size/self.num_chunk
		chunk_num = 0
		sys.stdout.write("chunk size: %d Bytes\n" %(chunk_size))
		fIN = open(self.infile, 'r')
		while True:
			start = fIN.tell()
			fIN.seek(chunk_size, 1)
			line = fIN.readline()
			if not line:
				break
			while not line.startswith(delim):
				line = fIN.readline()
			else:
				fIN.seek(-len(line), 1)
				chunk_num += 1
				task_q.put((start, fIN.tell()-start, chunk_num))
		task_q.put((start, None, chunk_num+1))
		fIN.close()
		return

	def fq2fa(self, fq_seq, fa_file):
		''' convert FastQ to Fasta '''
		i = 0
		fFA = open(fa_file, 'w')
		num_seq = 0
		for line in fq_seq.splitlines():
			if i % 4 == 0:
				fFA.write(">" + line[1:] + "\n")
				i += 1
			elif i % 4 == 1:
				fFA.write(line + "\n")
				i += 1
			elif i % 4 == 2:
				i += 1
			else:
				num_seq += 1
				i = 0
		sys.stdout.write("%s\tFastq to Fasta conversion ... %d\t[Done]\n" %(multiprocessing.current_process().name, num_seq))
		fFA.close()
		return

	def detect_format(self):
		''' detect format of input file '''
		cmd = "head -n 4 %s" %(self.infile)
		p = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
		pout = p.communicate()[0]
		for i, line in enumerate(pout.splitlines()):
			if i == 0:
				if line.startswith('>'):
					return "fa"
				elif line.startswith('@'):
					continue
			elif i == 2:
				if line.startswith('+'):
					return "fq"
		return None

	def start(self):
		self.__run()

	def __run(self):
		''' start the whole process '''
		if not os.path.exists(self.infile):
			sys.stderr.write("[%s] Error: Cannot find the input file %s\n" %(os.path.basename(__file__), self.infile))
			sys.exit(1)

		outdir = os.path.dirname(self.outp)
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		# detect format of input file
		fmt = self.detect_format()
		if fmt == None:
			sys.stderr.write("[%s] Error: Cannot decide format of the input file\n" %(os.path.basename(__file__)))
			sys.exit(1)

		try:
			infile_size = os.path.getsize(self.infile)
			sys.stdout.write("total fasta file: %d Bytes\n" %(infile_size))
		except os.error as e:
			sys.stderr.write("%s\n" %(e))
			sys.exit(1)

		task_q = multiprocessing.JoinableQueue()
		blast_q = multiprocessing.Queue()
		err_q = multiprocessing.Queue()
		chunk_processes = []
		for i in range(self.num_chunk):
			process = multiprocessing.Process(target=self.sanitize, args=(task_q, blast_q, err_q, fmt))
			process.daemon = True
			process.start()
			chunk_processes.append(process)
		self.chunkize_fasta(infile_size, fmt, task_q)

		try:
			task_q.join()
		except KeyboardInterrupt:
			sys.stderr.write("Terminated unexpectedly by keyboard\n")
			sys.exit()
		else:
			if not err_q.empty():
				while not err_q.empty():
					pid, outfile, returncode, err_msg = err_q.get_nowait()
					sys.stderr.write("[Error] %s gets sick and stop to work (BLAST EXIT %d: %s)\n" %(pid, returncode, err_msg))
				sys.exit(1)
			else:
				for process in chunk_processes:
					process.terminate()
					process.join()

				blast_files = []
				cat_cmd = "cat "
				if not blast_q.empty():
					while not blast_q.empty():
						cat_cmd += "%s " %(blast_q.get_nowait())

					try:
						pcat = subprocess.Popen(shlex.split(cat_cmd), stdout = open(self.outp + ".all.blast", 'w'))
						pcatout, pcaterr = pcat.communicate()
					except (OSError, ValueError) as e:
						sys.stderr.write(e + "\n")
						sys.stderr.write("[Error] Fail to concatenate BLAST files\n")
						sys.exit(1)
					sys.stdout.write("Programs finishes\n")
				else:
					sys.stderr.write("[Error] no or incomplete BLAST results found\n")
					sys.exit(1)
