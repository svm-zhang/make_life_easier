import os
import sys

def check_file_if_exists(path):
	if not os.path.exists(path):
		sys.stdout.write("[IO Error] Cannot find file %s\n" %(path))
		sys.exit(1)

def check_files_if_exist(*paths):
	for path in paths:
		check_file_if_exists(path)

def make_dir_if_necessary(path):
	if not os.path.exists(path):
		os.makedirs(path)

def make_dirs_if_necessary(*paths):
	for dir in paths:
		make_dir_if_necessary(dir)
