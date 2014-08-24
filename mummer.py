import os
import sys
import multiprocessing
import subprocess
import shlex
import glob

from io_simo import check_files_if_exist
from io_simo import make_dir_if_necessary
from filesmasher import FaSmasher
from ftsensor import FtSensor

class Mummer:
	def __init__(self, args):
		self.outp = os.path.realpath(args.outp)
		self.num_chunk = args.num_chunk
		self.qry = os.path.realpath(args.query)
		self.db = os.path.realpath(args.db)
		self.maxgap = args.maxgap
		self.mincluster = args.mincluster
		self.anchor_match = args.anchor_match
		if args.smashdb:
			self.smashwhich = "db"
		elif args.smashqry:
			self.smashwhich = "qry"
		self._check_cmd()
		self._print_cmd()

	def _check_cmd(self):
		check_files_if_exist(self.db, self.qry)
		make_dir_if_necessary(os.path.dirname(self.outp))

		nucmer_test = "nucmer -h"
		try:
			p = subprocess.Popen(shlex.split(nucmer_test), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
			pout, perr = p.communicate()
		except (OSError, ValueError) as e:
			sys.stderr.write("[mleFA Nucmer] Error: cannot find nucmer executable\n")
			sys.exit(1)

	def _print_cmd(self):
		sys.stdout.write("[mleFA Nucmer] output prefix: %s\n" %(os.path.realpath(self.outp)))
		sys.stdout.write("[mleFA Nucmer] num chunks: %d\n" %(self.num_chunk))
		sys.stdout.write("[mleFA Nucmer] reference: %s\n" %(os.path.realpath(self.db)))
		sys.stdout.write("[mleFA Nucmer] maximum gap: %d\n" %(self.maxgap))
		sys.stdout.write("[mleFA Nucmer] use of anchor matches: %s\n" %(self.anchor_match))
		sys.stdout.write("\n")
		return

	def _smash(self, file, chunk_start, chunk_end, smashed_db):
		with open(file, 'r') as fIN:
			fIN.seek(chunk_start)
			if chunk_end:
				seq = fIN.read(chunk_end)
			else:
				seq = fIN.read()
		fDB = open(smashed_db, 'w')
		fDB.write(seq)
		fDB.close()

	def _runucmer(self, task_q, results_q, err_q):
		while True:
			try:
				chunk_start, chunk_end, chunk_num = task_q.get()
				sys.stdout.write("%s\t%s\t%s\t%d\n" %(multiprocessing.current_process().name, chunk_start, chunk_end, chunk_num))
				smashed_db = self.outp + ".%d.fasta" %(chunk_num)
				if self.smashwhich == "db":
					self._smash(self.db, chunk_start, chunk_end, smashed_db)
				elif self.smashwhich == "qry":
					self._smash(self.qry, chunk_start, chunk_end, smashed_db)
				outp_delta = self.outp + ".%d" %(chunk_num)
				nucmer_cmd = "nucmer --%s -g %d -c %d -p %s %s %s" %(self.anchor_match, self.maxgap, self.mincluster, outp_delta, smashed_db, self.qry)
				sys.stdout.write(multiprocessing.current_process().name + "\t" + nucmer_cmd + "\n")
				try:
					pnucmer= subprocess.Popen(shlex.split(nucmer_cmd), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
					pnucmer_out, pnucmer_err = pnucmer.communicate()
				except (OSError, ValueError) as e:
					sys.stderr.write("%s\n" %(pnucmer_err))
					sys.stderr.write("%s\n" %(e))
					err_q.put((multiprocessing.current_process().name, out_delta+".delta", pnucmer.returncode))
				else:
					if os.path.exists(outp_delta+".delta") and os.path.getsize(outp_delta+".delta") > 0:
						out_coords = self.outp + ".%d.coords" %(chunk_num)
						coords_cmd = "show-coords -lcoHT %s " %(outp_delta+".delta")
						sys.stdout.write(multiprocessing.current_process().name + "\t" + coords_cmd + " %s" %(out_coords) + "\n")
						try:
							pcoords = subprocess.Popen(shlex.split(coords_cmd), stdout = open(out_coords, 'w'))
							pcoords_err = pcoords.communicate()[1]
						except (OSError, ValueError) as e:
							sys.stderr.write("%s\n" %(e))
							err_q.put((multiprocessing.current_process().name, out_coords, pcoords.returncode))
						else:
								results_q.put(out_coords)
					else:
						sys.stderr.write("[mleFA Nucmer] Error: no delta file found after Nucmer\n")
						err_q.put((multiprocessing.current_process().name, outp_delta+".delta", 1))
			finally:
				task_q.task_done()

	def _check_fmt(self, file):
		fmt = FtSensor(file).getfiletype()
		if fmt != "fa":
			sys.stderr.write("[mleFA Nucmer] Error: %s is not in FASTA format\n" %(file))
			sys.exit(1)

	def start(self):
		self._run()

	def _run(self):
		''' start the whole process '''
		if self.smashwhich == "db":
			self._check_fmt(self.db)
		elif self.smashwhich == "qry":
			self._check_fmt(self.qry)

		task_q = multiprocessing.JoinableQueue()
		results_q = multiprocessing.Queue()
		err_q = multiprocessing.Queue()
		chunk_processes = []
		for i in range(self.num_chunk):
			process = multiprocessing.Process(target=self._runucmer, args=(task_q, results_q, err_q))
			process.daemon = True
			process.start()
			chunk_processes.append(process)
		if self.smashwhich == "db":
			FaSmasher(self.db, self.num_chunk, self.outp, task_q).getchunk()
		elif self.smashwhich == "qry":
			FaSmasher(self.qry, self.num_chunk, self.outp, task_q).getchunk()

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
					else:
						sys.stdout.write("[mleFA Nucmer]\tcleaning up ...\n")
						for file in glob.glob(os.path.join(os.path.dirname(self.outp), "*.mgaps")):
							try:
								os.remove(file)
							except OSError as e:
								sys.stderr.write("%s\n" %(e))
								sys.exit(1)
						sys.stdout.write("[mleFA Nucmer] Programs finishes\n")
				else:
					sys.stderr.write("[mleFA Nucmer] Error: no COORDs files found\n")
					sys.exit(1)
