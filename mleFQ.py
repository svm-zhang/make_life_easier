
import argparse
import os
import sys


class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
	def _format_action(self, action):
		flag = 0
		parts = super(argparse.RawDescriptionHelpFormatter, self)._format_action(action)
		if action.nargs == argparse.PARSER:
			sub_cmd = "\n"
			for i, part in enumerate(parts.split("\n")):
				if i == 0:
					continue
				else:
					if flag == 1:
						sub_cmd += 4*" "+ " ".join(filter(None, part.split(" "))) + "\n"
						flag = 0
						continue
					if len(part.split(" ")) > 4:
						if len(part.split(" ")[4]) > 7:
							sub_cmd += " ".join(part.split(" ")[0:5])
							flag = 1
						else:
							sub_cmd += " ".join(part.split(" ")[0:5]) + (9-len(part.split(" ")[4])+4)*" " + " ".join(filter(None, part.split(" "))) + "\n"
			return sub_cmd
		else:
			return parts

def handle_cmd():
	main_parser = argparse.ArgumentParser(description="FastQ main portal", formatter_class=SubcommandHelpFormatter)
	sub_parsers = main_parser.add_subparsers(title="Commands", metavar="[command]", dest="command")

	# Arguments for blast
	runblast_parser = sub_parsers.add_parser("blast", help="sharply spot dirts in your data")
	runblast_parser.add_argument("-in", metavar="FILE", dest="infile", required=True, help="input file in either FASTQ or FASTA format")
	runblast_parser.add_argument("-outp", metavar="PREFIX", dest="outp", required=True, help="output prefix string")
	runblast_parser.add_argument("-nchunks", metavar="NUM", dest="num_chunk", type=int, required=True, help="number of chunks FASTA file to split. This also specifies the number of BLAST jobs running simultaneously. nchunks * nthreads = total number of processors to be requested")
	blast_group = runblast_parser.add_argument_group(description="Arugments for running BLAST")
	blast_group.add_argument("-db", metavar="FILE", dest="db", required=True, help="database for running BLAST")
	blast_group.add_argument("-nthreads", metavar="NUM", dest="nthreads", type=int, default=1, help="number of threads for running BLAST")
	blast_group.add_argument("-evalue", metavar="NUM", dest="evalue", default=10, help="E-value threshold for running BLAST. 1e-2 is supported")
	blast_group.add_argument("-ofmt", metavar="NUM", dest="outfmt", default=6, type=int, help="specify an integer for BLAST -outfmt option [6]")
	runblast_parser.set_defaults(func=runblast)

	# Arguments for clean
	clean_parser = sub_parsers.add_parser("clean", help="remove reads with BLAST hits obtained from mleFQ.py blast command")
	clean_parser.add_argument("-fq1", metavar="FILE", dest="fq1", required=True, help="FastQ file of first-end to be cleaned")
	clean_parser.add_argument("-fq2", metavar="FILE", dest="fq2", required=True, help="FastQ file of second-end to be cleaned")
	clean_parser.add_argument("-b1", metavar="FILE", dest="b1", required=True, help="BLAST result in m6 format for fq1")
	clean_parser.add_argument("-b2", metavar="FILE", dest="b2", required=True, help="BLAST result in m6 format for fq2")
	clean_parser.add_argument("-o1", metavar="FILE", dest="o1", required=True, help="output FastQ file for the first-end")
	clean_parser.add_argument("-o2", metavar="FILE", dest="o2", required=True, help="output FastQ file for the second-end")

	# Arguments for fq2fa
	fq2fa_parser = sub_parsers.add_parser("fq2fa", help="convert given FastQ to Fasta")
	fq2fa_parser.add_argument("-f", metavar="(Fq Fa)", dest="files", nargs='+', help="specify at least one FastQ and one Fasta files, separated by single space. For example: a.fastq a.fasta b.fastq b.fasta.")
	fq2fa_parser.add_argument("-nproc", metavar="INT", dest="nproc", type=int, default=1, help="specify number of conversion jobs running simultaneously [1]")
	fq2fa_parser.set_defaults(func=fq2fa)

	# Arguments for smash
	smash_parser = sub_parsers.add_parser("smash", help="smash file into pieces")
	smash_parser.add_argument("-in", metavar="FILE", dest="infile", required=True, help="specify FastQ file to be split")
	smash_parser.add_argument("-nchunks", metavar="INT", dest="num_chunk", type=int, default=2, help="specify number of chunks/blocks to split [2]")
	smash_parser.add_argument("-nproc", metavar="INT", dest="nproc", type=int, default=1, help="specify number of smashing jobs running simultaneously [1]")
	smash_parser.add_argument("-p", metavar="PREFIX", dest="outp", help="output prefix")
	smash_parser.set_defaults(func=smash)

	# Arguments for order
	fixp_parser = sub_parsers.add_parser("order", help="fixing reads order")
	fixp_parser.add_argument("-save", dest="mem_save", action="store_true", help="specify whether use memory saving mode. Warning: this runs very slow")
	fixp_parser.set_defaults(func=fixpairing)

	fastqc_parser = sub_parsers.add_parser("fastqc", help="run FastQC on give file(s)")

	interleave_parser = sub_parsers.add_parser("mergeFQ", help="create interleaved Fastq file from pair-end data")

	return main_parser.parse_args()

def fq2fa(args):
	if args.nproc < 1:
		sys.stderr.write("[CMD PARSER] number of processes must be larger than 0\n")
		sys.exit(1)
	elif args.nproc == 1:
		from fx2fy import Fq2Fa
		i = 0
		while i < len(args.files)-1:
			Fq2Fa(args.files[i], args.files[i+1]).start()
			i += 2
	else:
		from fx2fy import Fqs2Fas
		Fqs2Fas(args.files, args.nproc).start()

def runblast(args):
	from sharpeye import SharpEye
	SharpEye(args).start()

def fixpairing(args):
	pass

def smash(args):
	from filesmasher import FqSmasher
	import multiprocessing
	task_q = multiprocessing.JoinableQueue()
	FqSmasher(args.infile, args.num_chunk, args.outp, task_q).start(args.nproc)

def interleave(args):
	pass

def call_fastqc(args):
	pass

def main():
	args = handle_cmd()
	args.func(args)

if __name__ == "__main__":
	main()
