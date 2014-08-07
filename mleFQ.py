
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

	# Arguments for sharpeye
	sharpeye_parser = sub_parsers.add_parser("sharpeye", help="sharply spot dirts in your data")
	sharpeye_parser.add_argument("-in", metavar="FILE", dest="infile", required=True, help="input file in either FASTQ or FASTA format")
	sharpeye_parser.add_argument("-outp", metavar="PREFIX", dest="outp", required=True, help="output prefix string")
	sharpeye_parser.add_argument("-nchunks", metavar="NUM", dest="num_chunk", type=int, required=True, help="number of chunks FASTA file to split. This also specifies the number of BLAST jobs running simultaneously. nchunks * nthreads = total number of processors to be requested")
	blast_group = sharpeye_parser.add_argument_group(description="Arugments for running BLAST")
	blast_group.add_argument("-db", metavar="FILE", dest="db", required=True, help="database for running BLAST")
	blast_group.add_argument("-nthreads", metavar="NUM", dest="nthreads", type=int, default=1, help="number of threads for running BLAST")
	blast_group.add_argument("-evalue", metavar="NUM", dest="evalue", default=10, help="E-value threshold for running BLAST. 1e-2 is supported")
	blast_group.add_argument("-ofmt", metavar="NUM", dest="outfmt", default=6, type=int, help="specify an integer for BLAST -outfmt option [6]")

	fixp_parser = sub_parsers.add_parser("order", help="fixing reads order")

	fq2fa_parser = sub_parsers.add_parser("fq2fa", help="convert given FastQ to Fasta")

	smash_parser = sub_parsers.add_parser("smash", help="smash file into pieces")

	fastqc_parser = sub_parsers.add_parser("fastqc", help="run FastQC on give file(s)")

	interleave_parser = sub_parsers.add_parser("mergeFQ", help="create interleaved Fastq file from pair-end data")

	return main_parser.parse_args()

def main():
	args = handle_cmd()

	if args.command == "sharpeye":
		sys.path.insert(0, os.getcwd())
		from sharpeye import SharpEye
		SharpEye(args).start()

if __name__ == "__main__":
	main()
