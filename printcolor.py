class PrintColor:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'

	def disable(self):
		self.HEADER = ''
		self.OKBLUE = ''
		self.OKGREEN = ''
		self.WARNING = ''
		self.FAIL = ''
		self.ENDC = ''

print PrintColor.WARNING + "Warning: No active frommets remain. Continue?" + PrintColor.ENDC

print "\033[1;34mGREEN TEXT\033[0m"
print "\033[4;34mGREEN TEXT\033[0m"
print "\033[5;34mGREEN TEXT\033[0m"
print "\033[7;34mGREEN TEXT\033[0m"
