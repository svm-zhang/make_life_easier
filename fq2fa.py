import sys

class Fq2Fa:
	def __init__(self, fq="", fa=""):
		self.fa = fa
		self.fq = fq

	def start(self):
		return self.__run()

	def __run(self):
		print "parent run() function"
		num_seq = 0
		fFA = open(self.fa, 'w')
		with open(self.fq, 'r') as fFQ:
			for i, line in enumerate(fFQ):
				if i % 4 == 0:
					fFA.write(">" + line[1:])
				elif i % 4 == 1:
					fFA.write(line)
				else:
					num_seq += 1
				fFA.flush()
		return num_seq

class FqSeq2Fa(Fq2Fa):
	def __init__(self, fq_seq="", fa=""):
		self.fq_seq = fq_seq
		self.fa = fa
		self.num_seq = 0

	def start(self):
		self.__run()
		return self.num_seq

	def __run(self):
		i = 0
		fFA = open(self.fa, 'w')
		for line in self.fq_seq.splitlines():
			if i % 4 == 0:
				fFA.write(">" + line[1:] + "\n")
				i += 1
			elif i % 4 == 1:
				fFA.write(line + "\n")
				i += 1
			elif i % 4 == 2:
				i += 1
			else:
				self.num_seq += 1
				i = 0
			fFA.flush()
		fFA.close()
		return

#fq = sys.argv[1]
#fa = sys.argv[2]
#Fq2Fa(fq, fa).start()
