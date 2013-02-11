#!/usr/bin/env python
import sys
import numpy


#defines variables
the_infile = ""
the_value = 0
luminthreshold = 1  #checked if luminousity of the star is above this value. This is to make sure the module returns a giant star, not a main sequence one.


#global mass, lumin, temp, matchlumin, matchtemp

matchlumin = numpy.empty(0)
matchtemp = numpy.empty(0)


def readmod(infile):            #Reads in model data
#	mass = numpy.empty(0)
#	lumin = numpy.empty(0)
#	temp = numpy.empty(0)
	a = numpy.loadtxt(infile,delimiter=',',skiprows=13)
	mass = a[:,2]
	print mass
	lumin=a[:,4]
	temp=a[:,5]

#	with open(infile, 'rb') as csvfile1:                       
#		datain = csv.reader(csvfile1, delimiter=',')
#		raise UserWarning
#		for line in datain:
#			if "#" not in line[0]:
#				mass=numpy.append(mass,float(line[3]))     #Reads in the .csv with the isochrone model data.
#				lumin=numpy.append(lumin,float(line[4]))
#				temp=numpy.append(temp,float(line[5]))
	return mass, lumin, temp

def findall(value, mass, lumin, temp):
	matchlumin = numpy.empty(0)
	matchtemp = numpy.empty(0)
	nearest = find_nearest(mass, value)
	matchlumin = numpy.append(matchlumin, lumin[nearest])
	matchtemp = numpy.append(matchtemp, temp[nearest])
	if lumin[nearest] < luminthreshold:
		print "Not a giant"
	return matchlumin, matchtemp


def check(luminThreshold = 1):
	for i in numpy.arange(len(matchlumin)):
		if matchlumin[i] >= luminThreshold:
			return matchlumin[i], matchtemp[i]
		elif i+1 == len(matchlumin):
			print "Not a giant star"
			return matchlumin[0], matchtemp[0]


def set_data(the_file, the_mass):
	global the_infile
	global the_value
	the_infile = the_file
	the_value = the_mass


def find_nearest(array,value):
	idx = (numpy.abs(array-value)).argmin()
	return idx
#	return array[idx]


def run():
	mass, lumin, temp=readmod(the_infile)
	matchlumin, matchtemp = findall(the_value, mass, lumin, temp)
 	return matchlumin,matchtemp

'''
set_data("test_iso2.csv", 1.37)
#readmod(the_infile)
#findall(the_value)
result = run()
print result
'''

