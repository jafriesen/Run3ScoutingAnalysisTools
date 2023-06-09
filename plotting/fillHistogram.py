import ROOT
from ROOT import *
from math import *
#!/usr/bin/python3                                                                                                                                            
#-----------------------------------------------                                                                                                             
import sys, os, pwd, subprocess, glob, fnmatch
import optparse, shlex, re
import time
from time import gmtime, strftime
import math
from array import array

#define function for parsing options
def parseOptions():

	usage = ('usage: %prog [options]\n'
			+ '%prog -h for help')
	parser = optparse.OptionParser(usage)

	# input options
	parser.add_option('-i', '--input', dest='INPUT', type='string', help='input file')
	parser.add_option('-o', '--output', dest='OUTPUT', type='string', help='output file')
	parser.add_option('-n', '--njobs', dest='NJOBS', type=int, help='njobs')
	parser.add_option('-j', '--job', dest='JOB', type=int, help='job')
	parser.add_option('-r', '--run', dest='RUN', type='string', help='run2 or run3')
	parser.add_option('-c', '--condition', dest='CONDITION', type='string', help='full or jpsi')

	# store options and arguments as global variables
	global opt, args
	(opt, args) = parser.parse_args()

# define function for processing the external os commands
def processCmd(cmd, quite = 0):
	status, output = subprocess.getstatusoutput(cmd)
	if (status !=0 and not quite):
		print('Error in processing command:\n 	['+cmd+']')
		print('Output:\n   ['+output+'] \n')
	return output

def fillHistogram():
	
	global opt, args
	parseOptions()

	ROOT.gROOT.SetBatch()

	print(opt.INPUT)
	listDir = "/afs/cern.ch/user/j/jfriesen/CMSSW_12_4_2/src/Run3ScoutingAnalysisTools/FillHistogram"
	if ( opt.RUN == "run2" ) : 
		files = [ (opt.INPUT + "/000" + str(folder) + "/" + line) for folder in range(2) for line in open(listDir+"/list000" + str(folder) + "run2.txt","r") ]
	else :
		files = [ (opt.INPUT + "/000" + str(folder) + "/" + line) for folder in range(6) for line in open(listDir+"/list000" + str(folder) + ".txt","r") ]
	N = len(files)

	first = int(float(N)/float(opt.NJOBS)*float(opt.JOB-1))
	last = int(float(N)/float(opt.NJOBS)*float(opt.JOB))

	print(first, last)

	t = TChain("scoutingTree/tree")
	for i in range(len(files)):
		if (i<first or i>=last): continue
		print(files[i])
		t.Add(files[i])
		print(t.GetEntries())

	bins = 1200
	conditions = {
		"jpsi" : "vtxMatch==1 && mass > 2.95 && mass < 3.25 && vtxXError < 0.05 && vtxYError < 0.05 && vtxZError < 0.10 && vtxChi2/vtxNdof < 5",
		"full" : "vtxMatch==1 && vtxXError < 0.05 && vtxYError < 0.05 && vtxZError < 0.10 && vtxChi2/vtxNdof < 5",
		"id" : "vtxMatch==1 && vtxXError < 0.05 && vtxYError < 0.05 && vtxZError < 0.10 && vtxChi2/vtxNdof < 5 && muonID1 == 1 && muonID2 == 1"
	}

	h_lxy = TH1F("h_lxy","h_lxy",bins,0,120)
	print("drawing " + conditions[opt.CONDITION])
	t.Draw("Lxy>>h_lxy",conditions[opt.CONDITION])
	print("saving as "+str(opt.OUTPUT)+str(opt.JOB)+".root")
	h_lxy.SaveAs(str(opt.OUTPUT)+str(opt.JOB)+".root")

if __name__ == "__main__":
	fillHistogram()


