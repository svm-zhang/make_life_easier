import os
import sys
import multiprocessing
import argparse
import subprocess
import shlex
import tempfile

from io_simo import check_file_if_exists
from io_simo import make_dir_if_necessary
from seq2file import FqSeq2Fa
from filesmasher import FqSmasher
from filesmasher import FaSmasher
from ftsensor import FtSensor

class Blast:
	def __init__(self, args):
		self.infile = os.path.realpath(args.infile)
		self.outp = os.path.realpath(args.outp)
		check_file_if_exists(self.infile)
		make_dir_if_necessary(os.path.dirname(self.outp))
		self.num_chunk = args.num_chunk
		self.db = os.path.realpath(args.db)
		self.nthreads = args.nthreads
		self.evalue = args.evalue
		self.outfmt = args.outfmt
		self.print_cmd()

	def print_cmd(self):
		sys.stdout.write("[CMD PARSER] infile: %s\n" %(os.path.realpath(self.infile)))
		sys.stdout.write("[CMD PARSER] output prefix: %s\n" %(os.path.realpath(self.outp)))
		sys.stdout.write("[CMD PARSER] num chunks: %d\n" %(self.num_chunk))
		sys.stdout.write("[CMD PARSER] BLAST database: %s\n" %(os.path.realpath(self.db)))
		sys.stdout.write("[CMD PARSER] num threads: %d\n" %(self.nthreads))
		sys.stdout.write("[CMD PARSER] minimum E-Value: %s\n" %(self.evalue))
		sys.stdout.write("\n")
		return

	def _dirtspotter(self, task_q, blast_q, err_q, fmt):
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
				blast_cmd = "blastn -db %s -outfmt %d -out %s -num_threads %d -evalue %d " %(self.db, self.outfmt, out_blast, self.nthreads, self.evalue)
				if fmt == "fa":
					sys.stdout.write(multiprocessing.current_process().name + "\t" + blast_cmd + "\n")
					p = subprocess.Popen(shlex.split(blast_cmd), stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
					bout, berr = p.communicate(seq)
				elif fmt == "fq":
					fa_file = tempfile.mkstemp(suffix=".fasta", dir=os.path.dirname(self.outp), prefix=os.path.basename(self.outp))[1]
					num_seq = FqSeq2Fa(seq, fa_file).start()
					sys.stdout.write("%s\tFastq to Fasta conversion ... %d\t[Done]\n" %(multiprocessing.current_process().name, num_seq))
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

	def start(self):
		self._run()

	def _run(self):
		''' start the whole process '''
		# smell format of input file
		fmt, _, encoding = FtSensor(self.infile).bloodhound()
		if fmt == None:
			sys.stderr.write("[%s] Error: Cannot decide format of the input file\n" %(os.path.basename(__file__)))
			sys.exit(1)

		task_q = multiprocessing.JoinableQueue()
		blast_q = multiprocessing.Queue()
		err_q = multiprocessing.Queue()
		chunk_processes = []
		for i in range(self.num_chunk):
			process = multiprocessing.Process(target=self._dirtspotter, args=(task_q, blast_q, err_q, fmt))
			process.daemon = True
			process.start()
			chunk_processes.append(process)
		if fmt == "fq":
			FqSmasher(self.infile, self.num_chunk, self.outp, task_q).getchunk()
		elif fmt == "fa":
			FaSmasher(self.infile, self.num_chunk, self.outp, task_q).getchunk()

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
