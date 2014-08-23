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
	runucmer_parser = sub_parsers.add_parser("nucmer", help="sharply spot dirts in your data")
	runucmer_parser.add_argument("-in", metavar="FILE", dest="infile", required=True, help="input file in FASTA format")
	runucmer_parser.add_argument("-outp", metavar="PREFIX", dest="outp", required=True, help="output prefix string")
	runucmer_parser.add_argument("-nchunks", metavar="NUM", dest="num_chunk", type=int, required=True, help="number of chunks FASTA file to split. This also specifies the number of BLAST jobs running simultaneously. nchunks * nthreads = total number of processors to be requested")
	nucmer_group = runucmer_parser.add_argument_group(description="Arugments for running BLAST")
	nucmer_group.add_argument("-qry", metavar="FILE", dest="query", required=True, help="query sequence")
	nucmer_group.add_argument("-db", metavar="FILE", dest="db", required=True, help="database for running BLAST")
	nucmer_group.add_argument("-maxgap", metavar="INT", dest="maxgap", type=int, default=90, help="specify -g|maxgap argument for nucmer")
	nucmer_group.add_argument("-anchor", dest="anchor_match", choices=["maxmatch"], default="maxmatch", help="specify how nucmer use anchor matches")
	runucmer_parser.set_defaults(func=runucmer)

	return main_parser.parse_args()

def runucmer(args):
	from mummer import Mummer
	Mummer(args).start()


def main():
	args = handle_cmd()
	args.func(args)

if __name__ == "__main__":
	main()
