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
	runucmer_parser = sub_parsers.add_parser("nucmer", help="split either large database or query and run nucmer on each")
	runucmer_parser.add_argument("-outp", metavar="PREFIX", dest="outp", required=True, help="output prefix string")
	runucmer_parser.add_argument("-nchunks", metavar="NUM", dest="num_chunk", type=int, required=True, help="number of chunks FASTA file to split. This also specifies the number of BLAST jobs running simultaneously. nchunks * nthreads = total number of processors to be requested")
	runnucmer_mutual_group = runucmer_parser.add_mutually_exclusive_group(required=True)
	runnucmer_mutual_group.add_argument("-smashdb", action="store_true", help="specify to smash database")
	runnucmer_mutual_group.add_argument("-smashqry", action="store_true", help="specify to smash query")
	nucmer_group = runucmer_parser.add_argument_group(description="Arugments for running Nucmer")
	nucmer_group.add_argument("-qry", metavar="FILE", dest="query", required=True, help="query sequence")
	nucmer_group.add_argument("-db", metavar="FILE", dest="db", required=True, help="database file in Fasta")
	nucmer_group.add_argument("-g", metavar="INT", dest="maxgap", type=int, default=90, help="specify -g|maxgap argument for nucmer [90]")
	nucmer_group.add_argument("-c", metavar="INT", dest="mincluster", type=int, default=65, help="specify -c|mincluster argument for nucmer [65]")
	nucmer_group.add_argument("-a", metavar="STR", dest="anchor_match", choices=["maxmatch", "mum", "mumreference"], default="maxmatch", help="specify how nucmer uses anchor matches. choices: maxmatch, mum, mumreference [maxmatch]")
	runucmer_parser.set_defaults(func=runucmer)

	# Argument for smash
	smash_parser = sub_parsers.add_parser("smash", help="smash files into chunks")
	smash_parser.add_argument("-in", metavar="FILE", dest="infile", required=True, help="specify input file in Fasta format")
	smash_parser.add_argument("-nchunks", metavar="INT", dest="num_chunk", type=int, default=2, help="specify number of chunks the given file will be smashed [2]")
	smash_parser.add_argument("-p", metavar="PREFIX", dest="outp", required=True, help="specify the prefix of output file")
	smash_parser.add_argument("-nproc", metavar="INT", dest="nproc", type=int, default=1, help="specify number of smashing jobs running simultaneously [1]")
	smash_parser.set_defaults(func=smash)

	return main_parser.parse_args()

def smash(args):
	from filesmasher import FaSmasher
	task_q = multiprocessing.JoinableQueue()
	FaSmasher(args.infile, args.num_chunk, args.outp, task_q).start(args.nproc)

def runucmer(args):
	from mummer import Mummer
	Mummer(args).start()

def main():
	args = handle_cmd()
	args.func(args)

if __name__ == "__main__":
	main()
