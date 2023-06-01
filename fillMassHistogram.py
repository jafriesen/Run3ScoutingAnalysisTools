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
	parser.add_option('-c', '--condition', dest='CONDITION', type='string', help='selection conditions')

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
		l1condition = "(l1Result[0] > 0 || l1Result[1] > 0 || l1Result[2] > 0 || l1Result[3] > 0)"
	else :
		files = [ (opt.INPUT + "/000" + str(folder) + "/" + line) for folder in range(6) for line in open(listDir+"/list000" + str(folder) + ".txt","r") ]
		l1condition = "(l1Result[0] > 0 || l1Result[1] > 0 || l1Result[2] > 0 || l1Result[3] > 0 || l1Result[4] > 0 || l1Result[5] > 0)"

	condition = l1condition if (opt.CONDITION == "l1") else ""
	vtxSelection = "vtxMatch==1 && vtxXError < 0.05 && vtxYError < 0.05 && vtxZError < 0.10 && vtxChi2/vtxNdof < 5"
	if (opt.CONDITION == "lxy20"): condition = "Lxy > 20 && " + vtxSelection
	if (opt.CONDITION == "lxy11to20"): condition = "Lxy > 11 && Lxy < 20 && " + vtxSelection
	if (opt.CONDITION == "lxy7to11"): condition = "Lxy > 7 && Lxy < 11 && " + vtxSelection

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

	mu_mass = 105.658/1000;
	nbins = 1500
	bin_edges = [0] + [mu_mass*(1.01)**i for i in range(nbins+1)]
	h = TH1F("h_mass","h_mass",nbins,array('d',bin_edges))
	
	print("drawing ")
	t.Draw("mass>>h_mass", condition)
	print("saving as "+str(opt.OUTPUT)+str(opt.JOB)+".root")
	h.SaveAs(str(opt.OUTPUT)+str(opt.JOB)+".root")

if __name__ == "__main__":
	fillHistogram()


