#!/usr/bin/python3                                                                                                                               
#-----------------------------------------------                                                                                                             
import sys, os, pwd, subprocess, glob, fnmatch
import optparse, shlex, re
import time
from time import gmtime, strftime
import math
from array import array


JOB_PREFIX = """#!/bin/sh
cd %(CMSSW_BASE)s/src
export SCRAM_ARCH=%(SCRAM_ARCH)s
eval `scramv1 runtime -sh`
cd -
cp %(SCRIPTNAME)s .
"""

CONDOR_TEMPLATE = """executable = %(EXE)s
arguments = $(ProcId)
output                = %(TASK)s.$(ClusterId).$(ProcId).out
error                 = %(TASK)s.$(ClusterId).$(ProcId).err
log                   = %(TASK)s.$(ClusterId).log
use_x509userproxy = true
transfer_output_files   = ""
# Send the job to Held state on failure.
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
# Periodically retry the jobs every 10 minutes, up to a maximum of 5 retries.
periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)
+JobFlavour = "longlunch"
queue %(NUMBER)s
"""

#define function for parsing options
def parseOptions():

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-i', '--input', dest='INPUT', type='string', help='input file')
    parser.add_option('-o', '--output', dest='OUTPUT', type='string', help='output file')
    parser.add_option('-e', '--outdir', dest='OUTDIR', type='string', help='output directory')
    parser.add_option('-r', '--run', dest='RUN', type='string', help='run2 or run3')
    parser.add_option('-c', '--condition', dest='CONDITION', type='string', help='full or jpsi')
    parser.add_option('-t', '--taskname', dest='TASKNAME', type='string', help='task name')    
    parser.add_option('-s', '--scriptname', dest='SCRIPTNAME', type='string', help='script name')    
    parser.add_option('-n', '--njobs', dest='NJOBS', type=int, help='njobs')
    parser.add_option('-d', '--dryrun', dest='DRY', action='store_true', help='dry run only')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# define function for processing the external os commands
def processCmd(cmd, quite = 0):
    status, output = subprocess.getstatusoutput(cmd)
    if (status !=0 and not quite):
        print('Error in processing command:\n   ['+cmd+']')
        print('Output:\n   ['+output+'] \n')
    return output

def find_files(directory, pattern):
  for root, dirs, files in os.walk(directory):
    for basename in files:
      if fnmatch.fnmatch(basename, pattern):
        filename = os.path.join(root, basename)
        yield filename

def submitFillHistogram():

  global opt, args
  parseOptions()

  startdir = os.getcwd()
  os.mkdir(opt.TASKNAME)
  os.chdir(startdir+"/"+opt.TASKNAME)

  outscriptname = 'condor_%s.sh' % opt.TASKNAME
  subfilename = 'condor_%s.sub' % opt.TASKNAME
  print('>> condor job script will be %s' % outscriptname)
  
  outscript = open(outscriptname, "w")
  job_settings = JOB_PREFIX % {
      'CMSSW_BASE': os.environ['CMSSW_BASE'],
      'SCRAM_ARCH': os.environ['SCRAM_ARCH'],
      'SCRIPTNAME': opt.SCRIPTNAME
      }
  outscript.write(job_settings)
  for i in range(1, opt.NJOBS+1):
    outscript.write('\nif [ $1 -eq %i ]; then\n' % i)
    line = "python3 "+str(opt.SCRIPTNAME)+" -i "+str(opt.INPUT)+" -o "+str(opt.OUTPUT)+" -n "+str(opt.NJOBS)+" -j "+str(i)+" -r "+str(opt.RUN)+" -c "+str(opt.CONDITION)
    outscript.write('  ' + line + '\n')
    outscript.write('fi \n')
  line = "cp "+opt.OUTPUT+"* "+opt.OUTDIR
  outscript.write(line + '\n')
  line = "rm "+opt.OUTPUT+"*"
  outscript.write(line + '\n')
  outscript.close()

  subfile = open(subfilename, "w")
  condor_settings = CONDOR_TEMPLATE % {
      'EXE': outscriptname,
      'TASK': opt.TASKNAME,
      'SCRIPTNAME': opt.SCRIPTNAME,
      'NUMBER': opt.NJOBS
      }

  subfile.write(condor_settings)
  subfile.close()

  if (not opt.DRY):
    processCmd('condor_submit %s' % (subfilename))

  os.chdir(startdir)

      
      
if __name__ == "__main__":
  submitFillHistogram()
