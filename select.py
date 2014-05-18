#!/usr/bin/python
import sys
import re
from string import *

f = open(sys.argv[1], 'r')

while True:
	line = f.readline()
	if len(line) == 0: # Zero length indicates EOF
		break
	line = re.sub(r'.:', "", line)
	arr = line.split(" ")
	sys.stdout.write(arr[0])
	for i in range(2, len(sys.argv)):
		sys.stdout.write(" %d:%.17g" %(i - 1, float(arr[int(sys.argv[i])])))
	print
		# Notice comma to avoid automatic newline added by Python
f.close() # close the fileprint sys.argv[1]


