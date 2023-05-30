#!/usr/bin/env python3

import ROOT, array, random, copy
from ROOT import TCanvas, TFile, TH1, TH1F, TF1, gSystem
from ROOT import *
import ROOT, array, CMSGraphics, CMS_lumi, random, copy
from array import array
import math
import os

ROOT.gROOT.SetBatch()

print("folder?")
folder = str(input())
hdir = "/eos/user/j/jfriesen/fillHistogram/" + folder

print("name?")
name = str(input())

file = TFile.Open(hdir +"/"+ name+"1.root","READ")
histo = file.Get(name)

for i in range(400-2) :
	#print("adding "+name+str(i+2)+".root")
	#print(histo.GetEntries())
	f = TFile.Open(hdir +"/"+ name+str(i+2)+".root","READ")
	histo.Add(f.Get(name))

print(histo.GetEntries())
histo.SaveAs(hdir+"/"+folder+".root")