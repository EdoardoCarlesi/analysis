#!/usr/bin/python
import sys, getopt

narg = len(sys.argv)

Years = float(sys.argv[1])
Total = int(sys.argv[2])
fname = sys.argv[3]
fname_out = sys.argv[4]

file_in = open(fname, "r")
file_out = open(fname_out, "w")

def z2Gyrs(z):
	x = float(z)
	result = (2 * Years) / (1 + (1+x)*(1+x))
	return result

file_out.write('#number(1)\tz(2)\tGYrs(3)\n')

lines = file_in.readlines()

file_in.close()

time = []
num = []
zeta = []

for i in range (0, Total):
	columns =  (lines[i].strip().split())
	time.append((z2Gyrs(columns[1])))
	zeta.append(float(columns[1]))
	num.append(columns[0])

for i in range (0, Total):
	print >> file_out, "\t%s\t" % num[i], "%.3f\t" % zeta[i], "%.3f\t" % time[i]
