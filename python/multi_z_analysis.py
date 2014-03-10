#!/usr/bin/python
import sys, getopt
from scipy.interpolate import interp1d

narg = len(sys.argv)

fname_out = sys.argv[1]

# redshift list
fz_list = sys.argv[2]

# general statistics file
f1_list = sys.argv[3]

# File type to be read in
list_type = sys.argv[4]

header_lines = 1
rvir_gas = 0.2

f1_in = open(f1_list, "r")
fz_in = open(fz_list, "r")
f_out = open(fname_out, "w")

# This reads in from the file list and the corresponding redshifts
lines_1 = f1_in.readlines()
lines_z = fz_in.readlines()

f1_in.close()
fz_in.close()

TotFiles = len(lines_1)

zetas = []

if (list_type=="avg_profiles"):
	f_out.write('#r_vir(1)\tgas_f(2)\trho_dm(3)\n')
	rvir = []
	rhodm = []
	gasf = []

	rvir_temp = []
	rhodm_temp = []
	gasf_temp = []

elif (list_type=="numerical_mass_function"):
	mass = []
	num = []

	mass_temp = []
	num_temp = []

for i in range (0, TotFiles):
	columns_z =  (lines_z[i].strip().split())
	zetas.append(float(columns_z[0]))
	fname_temp = str(lines_1[i].strip())
	f_temp = open(fname_temp, "r")

	lines_temp = f_temp.readlines()
	TotLinesTemp = len(lines_temp)
	print "File = ", fname_temp, "z = %.3f" % zetas[i], " Lines= ", TotLinesTemp

	if (list_type=="avg_profiles"):
		for j in range (header_lines, TotLinesTemp):
			columns_f =  (lines_temp[j].strip().split())
			rvir_temp.append(float(columns_f[0]))
			rhodm_temp.append(float(columns_f[1]))
			gasf_temp.append(float(columns_f[4]))
		
		interp_r = interp1d(rvir_temp, gasf_temp)
		interp_rho = interp1d(rvir_temp, rhodm_temp)
		rvir.append(interp_r(rvir_gas))
		rhodm.append(interp_rho(rvir_gas))
	
		print "Rvir = %.2f" % rvir_gas, " g_f = %.3f " % rvir[i], " rhodm = %.3f " %rhodm[i] 
		print >> f_out, "%.3f\t\t" % zetas[i], "%.3f\t\t" % rvir[i], "%.3f\n" % rhodm[i]

	elif (list_type=="numerical_mass_function"):
		for j in range (header_lines, TotLinesTemp):
			columns_f =  (lines_temp[j].strip().split())
			mass_temp.append(float(columns_f[1]))
			num_temp.append(float(columns_f[2]))

print "Files read in and interpolated, output saved to ", fname_out


