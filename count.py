#!/usr/bin/python
import sys
from string import *

count = 0
if len(sys.argv) >= 2:
	f = sys.argv[1]
else:
	f = '\t1'
for line in sys.stdin:
	count += line.count(f)
print count
#print sys.argv[1].count(sys.argv[2])

