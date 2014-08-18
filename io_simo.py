import os
import sys

def check_file_if_exists(path):
	abs_path = os.path.realpath(path)
	if not os.path.exists(abs_path):
		sys.stdout.write("[IO Error] Cannot find file %s\n" %(abs_path))
		sys.exit(1)

def check_files_if_exist(*paths):
	for path in paths:
		check_file_if_exists(path)

def make_dir_if_necessary(path):
	if not os.path.exists(path):
		os.makedirs(path)

def make_dirs_if_necessary(*paths):
	for path in paths:
		make_dir_if_necessary(path)
