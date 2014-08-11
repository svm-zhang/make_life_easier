import os
import sys

#class IO_Simo(object) :
#	def __init__(self, *paths) :
#		self.paths = paths

#class FileIO(IO_Simo):
#	def __init__(self, *paths) :
#		super(FileIO, self).__init__(paths)
#		if len(self.paths) > 1:
#			self.check_files_if_exist(self.paths)
#		else:
#			self.check_file_if_exists(self.paths)

def check_file_if_exists(path):
	if not os.path.exists(path):
		sys.stdout.write("[IO Error] Cannot find file %s\n" %(path))
#		sys.exit(1)

def check_files_if_exist(paths):
	for path in self.paths:
		self.check_file_if_exists(path)

def make_dir_if_necessary(path):
	if not os.path.exists(path):
		os.makedirs(path)

def make_dirs_if_necessary(paths):
	for dir in self.paths:
		self.make_dir_if_necessary(dir)

#class DirIO(IO_Simo):
#	def __init__(self, *paths) :
#		super(FileIO, self).__init__(paths)
#		if len(self.paths) > 1:
#			self.make_dirs_if_necessary(self.paths)
#		else:
#			self.make_dir_if_necessary(self.paths)
