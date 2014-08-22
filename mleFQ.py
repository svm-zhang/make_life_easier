import argparse
import os
import sys
import multiprocessing

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
	clean_parser.add_argument("-in", metavar="FQ,BLAST", dest="inputs", action="append", help="input files. Can be specified multiple times. Example -in a1.fq,a1.blast -in a2.fq,a2.blast")
	clean_parser.add_argument("-nproc", metavar="INT", dest="nproc", type=int, default=1, help="Specify the number cleaning jobs running simultaneously")
	clean_parser.set_defaults(func=clean)

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
	fixp_parser = sub_parsers.add_parser("fixpe", help="fixing reads order")
	fixp_parser.add_argument("-fq1", metavar="FILE", dest="fq1", required=True, help="specify FastQ file of the first-end of a read pair")
	fixp_parser.add_argument("-fq2", metavar="FILE", dest="fq2", required=True, help="specify FastQ file of the second-end of a read pair")
	fixp_parser.add_argument("-outp", metavar="PREFIX", dest="outp", required=True, help="specify prefix for output file. PREFIX_r1.fastq, PREFIX_r2.fastq, and PREFIX_orphan.fastq will be output")
	fixp_parser.add_argument("-save", dest="mem_save", action="store_true", help="specify whether use memory saving mode. Warning: this runs very slow")
	fixp_parser.set_defaults(func=fixpairing)

	# Arguments for fastqc
	fastqc_parser = sub_parsers.add_parser("fastqc", help="run FastQC on give file(s)")
	fastqc_parser.add_argument("-target", metavar="DIR", dest="target", required=True, help="target directory where results will be put")
	fastqc_parser.add_argument("-nproc", metavar="INT", type=int, dest="nproc", default=1, help="number of FastQC jobs to run simultaneously. nproc*nthreads = tot_num_cores_to_request [1]")
	fastqc_group = fastqc_parser.add_argument_group(description="Arguments for running FastQC")
	fastqc_group.add_argument("-nthreads", metavar="INT", type=int, dest="nthreads", default=1, help="number of threads for each fastqc process [1]")
	fastqc_group.add_argument("-noextract", dest="noextract", action="store_true", help="specify to ask FastQC not to uncompress the output files [false/uncompress]")
	fastqc_group.add_argument("-kmer", metavar="INT", type=int, dest="kmer", default=5, choices=range(2,11), help="specify the length of Kmer to look for in FastQC [5]")
	fastqc_group.add_argument("-fastqc", metavar="FILE", dest="fastqc_prg", default="fastqc", help="specify a path to where FastQC is installed. If not specified, \"fastqc\" will be used [fastqc]")
	mutual_group = fastqc_parser.add_mutually_exclusive_group(required=True)
	mutual_group.add_argument("-source", metavar="DIR", dest="source", help="source directory where to-be-processed-fastq files are stored. Cannot use with -list at the same time")
	mutual_group.add_argument("-list", metavar="FILE", dest="source", help="a file with a list of fastq files to be processed")
	fastqc_parser.set_defaults(func=runfastqc)

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
	from blast import Blast
	Blast(args).start()

def clean(args):
	from clean import Clean
	task_q = multiprocessing.JoinableQueue()
	Clean(args.inputs, args.nproc, task_q).start()

def fixpairing(args):
	if not args.mem_save:
		from fixpe import FixPE
		FixPE(args).start()
	else:
		from fixpe import FixPE_Save
		FixPE_Save(args).start()

def smash(args):
	from filesmasher import FqSmasher
	task_q = multiprocessing.JoinableQueue()
	FqSmasher(args.infile, args.num_chunk, args.outp, task_q).start(args.nproc)

def interleave(args):
	pass

def runfastqc(args):
	from fastqc import FastQC
	FastQC(args).start()

def main():
	args = handle_cmd()
	args.func(args)

if __name__ == "__main__":
	main()
