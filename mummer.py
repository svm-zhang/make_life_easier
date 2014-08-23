import os
import sys
import multiprocessing
import argparse
import subprocess
import shlex
import tempfile

from io_simo import check_files_if_exist
from io_simo import make_dir_if_necessary
from filesmasher import FaSmasher
from ftsensor import FtSensor

class Mummer:
	def __init__(self, args):
		self.infile = os.path.realpath(args.infile)
		self.outp = os.path.realpath(args.outp)
		self.num_chunk = args.num_chunk
		self.qry = args.query
		self.db = os.path.realpath(args.db)
		self.maxgap = args.maxgap
		self.anchor_match = args.anchor_match
		self._check_cmd()
		self._print_cmd()

	def _check_cmd(self):
		check_files_if_exist(self.infile, self.db, self.qry)
		make_dir_if_necessary(os.path.dirname(self.outp))

		nucmer_test = "nucmer -h"
		try:
			p = subprocess.Popen(shlex.split(nucmer_test), stdout = subprocess.PIPE)
			pout, perr = p.communicate()
		except (OSError, ValueError) as e:
			sys.stderr.write("[mleFA Nucmer] Error: cannot find nucmer executable\n")
			sys.exit(1)

	def _print_cmd(self):
		sys.stdout.write("[mleFA Nucmer] infile: %s\n" %(os.path.realpath(self.infile)))
		sys.stdout.write("[mleFA Nucmer] output prefix: %s\n" %(os.path.realpath(self.outp)))
		sys.stdout.write("[mleFA Nucmer] num chunks: %d\n" %(self.num_chunk))
		sys.stdout.write("[mleFA Nucmer] reference: %s\n" %(os.path.realpath(self.db)))
		sys.stdout.write("[mleFA Nucmer] maximum gap: %d\n" %(self.maxgap))
		sys.stdout.write("[mleFA Nucmer] use of anchor matches: %s\n" %(self.anchor_match))
		sys.stdout.write("\n")
		return

	def _runucmer(self, task_q, results_q, err_q, fmt):
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
				db = self.outp + ".%d.fasta" %(chunk_num)
				fDB = open(db, 'w')
				fDB.write(seq)
				fDB.close()

				outp_delta = self.outp + ".%d" %(chunk_num)
				nucmer_cmd = "nucmer --%s -g %d -p %s %s %s" %(self.anchor_match, self.maxgap, outp_delta, db, self.qry)
				sys.stdout.write(multiprocessing.current_process().name + "\t" + nucmer_cmd + "\n")
				try:
					pnucmer= subprocess.Popen(shlex.split(nucmer_cmd), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
					pnucmer_out, pnucmer_err = pnucmer.communicate()
				except (OSError, ValueError) as e:
					sys.stdout.write("%s\n" %(e))
					err_q.put((multiprocessing.current_process().name, out_delta+".delta", pnucmer.returncode))
#				if pnucmer.returncode != 0:
				else:
					out_coords = self.outp + ".%d.coords" %(chunk_num)
					coords_cmd = "show-coords -lcoHT %s " %(outp_delta+".delta")
					sys.stdout.write(multiprocessing.current_process().name + "\t" + coords_cmd + "\n")
					try:
						pcoords = subprocess.Popen(shlex.split(coords_cmd), stdout = open(out_coords, 'w'))
						pcoords_err = pcoords.communicate()[1]
					except (OSError, ValueError) as e:
						sys.stdout.write("%s\n" %(e))
						err_q.put((multiprocessing.current_process().name, out_coords, pcoords.returncode))
#					if pcoords.returncode != 0:
#						sys.stdout.write("%s\n" %(pcoords_err))
					else:
						results_q.put(out_coords)
			finally:
				task_q.task_done()

	def start(self):
		self._run()

	def _run(self):
		''' start the whole process '''
		fmt = FtSensor(self.infile).getfiletype()
		if fmt != "fa":
			sys.stderr.write("[mleFA Nucmer] Error: input file is not in FASTA format\n")
			sys.exit(1)

		task_q = multiprocessing.JoinableQueue()
		results_q = multiprocessing.Queue()
		err_q = multiprocessing.Queue()
		chunk_processes = []
		for i in range(self.num_chunk):
			process = multiprocessing.Process(target=self._runucmer, args=(task_q, results_q, err_q, fmt))
			process.daemon = True
			process.start()
			chunk_processes.append(process)
		FaSmasher(self.infile, self.num_chunk, self.outp, task_q).getchunk()

		try:
			task_q.join()
		except KeyboardInterrupt:
			sys.stderr.write("Terminated unexpectedly by keyboard\n")
			sys.exit()
		else:
			if not err_q.empty():
				while not err_q.empty():
					pid, outfile, returncode = err_q.get_nowait()
					sys.stderr.write("[mleFA Nucmer] Error: %s gets sick and stop to work [EXIT %d]\n" %(pid, returncode))
				sys.exit(1)
			else:
				for process in chunk_processes:
					process.terminate()
					process.join()

				nucmer_files = []
				cat_cmd = "cat "
				if not results_q.empty():
					while not results_q.empty():
						cat_cmd += "%s " %(results_q.get_nowait())

					try:
						pcat = subprocess.Popen(shlex.split(cat_cmd), stdout = open(self.outp + ".all.coords", 'w'))
						pcatout, pcaterr = pcat.communicate()
					except (OSError, ValueError) as e:
						sys.stderr.write(e + "\n")
						sys.stderr.write("[mleFA Nucmer] Error: fail to concatenate COORDs files\n")
						sys.exit(1)
					sys.stdout.write("[mleFA Nucmer] Programs finishes\n")
				else:
					sys.stderr.write("[mleFA Nucmer] Error: no or incomplete COORDs results found\n")
					sys.exit(1)
