#!/usr/bin/python
import sys, getopt
from scipy.interpolate import interp1d

narg = len(sys.argv)

fname_out = sys.argv[1]

# redshift list
fz_list = sys.argv[2]

# general statistics file
f2_list = sys.argv[4]

# average profiles file
#f3_list = sys.argv[5]

f1_in = open(f1_list, "r")
#f2_in = open(f2_list, "r")
fz_in = open(fz_list, "r")
f_out = open(fname_out, "w")

# This reads in from the file list and the corresponding redshifts
lines_1 = f1_in.readlines()
lines_z = fz_in.readlines()

f1_in.close()
fz_in.close()

TotFiles = len(lines_1)

zetas = []

for i in range (0, TotFiles):
	columns_z =  (lines_z[i].strip().split())
	zetas.append(float(columns_z[0]))
	fname_temp = str(lines_1[i].strip())
	f_temp = open(fname_temp, "r")
	
	lines_temp = f_temp.readlines()
	TotLinesTemp = len(lines_temp)
	print "File = ", fname_temp, "z = ", zetas[i], " Lines= ", TotLinesTemp



"""
time = []
num = []
zeta = []

for i in range (0, Total):
	columns =  (lines[i].strip().split())
	time.append((z2Gyrs(columns[1])))
	zeta.append(float(columns[1]))
	num.append(columns[0])

file_out.write('#number(1)\tz(2)\tGYrs(3)\n')
for i in range (0, Total):
	print >> file_out, "\t%s\t" % num[i], "%.3f\t" % zeta[i], "%.3f\t" % time[i]

"""
