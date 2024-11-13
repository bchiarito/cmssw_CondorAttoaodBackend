#!/usr/bin/env python
#from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.fmk_atto.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
import glob

# command line options
parser = argparse.ArgumentParser(description="", usage="./%(prog)s INPUT OUTPUT")
# input/output
io_args = parser.add_argument_group('input/output options')
io_args.add_argument("input", metavar='INPUT',help="")
io_args.add_argument("output", metavar='OUTPUT', help="")
parser.add_argument("--dataset", default="MyDatasetName", help="name of dataset")
parser.add_argument("--proc", default=0, type=int, help="integer flag")
parser.add_argument("--drop", dest="branchsel", default=None, help=".txt file to drop branches")
parser.add_argument("--filter", dest="selection", default="None", metavar='CHOICE', help="")
parser.add_argument("-n", "--numEvents", dest="numEvents", default=-1, type=int, help="")
parser.add_argument("--outfile", default="out.root", help="name of final rootfile")
parser.add_argument("--report", default="report.txt", help="name of intermediate report txt file")
parser.add_argument("--analyzer", default="", help='')
datamc_options = parser.add_mutually_exclusive_group()
datamc_options.add_argument("--data", action="store_true", default=False, help="running on data")
datamc_options.add_argument("--mc", action="store_true", default=False, help="running on bkg mc")
datamc_options.add_argument("--sigRes", action="store_true", default=False, help="running on resonant signal mc")
datamc_options.add_argument("--sigNonRes", action="store_true", default=False, help="running on nonresonant signal mc")
args = parser.parse_args()


# import modules
from PhysicsTools.NanoAODTools.fmk_atto.simpleCounter import simpleCounter
from PhysicsTools.NanoAODTools.fmk_atto.simpleSelector import simpleSelector
from PhysicsTools.NanoAODTools.fmk_atto.analysisModule import analysisModule
from PhysicsTools.NanoAODTools.fmk_atto.mcHatModule import mcHatModule
from PhysicsTools.NanoAODTools.fmk_atto.baselineModule import baselineModule
from PhysicsTools.NanoAODTools.fmk_atto.PDFUncertaintyProducer import PDFUncertaintyProducer
from PhysicsTools.NanoAODTools.fmk_atto.triggerModule import triggerModule
from PhysicsTools.NanoAODTools.fmk_atto.ZttModule_muon_muon import zttModule

# temporary
tags = [0,2018,-1,1] ## data=0 mc=1 ,run year, mc number (-1 if data), antiiso = 1 iso = 0 (iso is normal run + same sign , anti iso is to get the same and opposite sign candidates in the anti iso region for ratio)

if os.path.isfile(args.input) and args.input[-5:] == '.root':
  files = [args.input]
elif os.path.isfile(args.input):
  files = []
  with open(args.input) as fi:
    for line in fi:
      # process line if .dat
      if args.input[-4:] == '.dat':
        newline = line.strip()
        i = newline.rfind('/')
        newline = newline[i+1:len(line)]
        files.append(newline)
      # dont process for .txt
      if args.input[-4:] == '.txt':
        files.append(line.strip())

# process arguments
datasetname = args.dataset
flag = args.proc
if args.data: datamc = 'data'
elif args.mc: datamc = 'mc'
elif args.sigRes: datamc = 'sigRes'
elif args.sigNonRes: datamc = 'sigNonRes'
else: raise SystemExit('ERROR: Must specify one of --data / --mc / --sigRes / --sigNonRes !')

# modules
modules = []
#modules += [PDFUncertaintyProducer(False, False)]
modules += [simpleCounter(args.report, "TotalEventsProcessed")]
modules += [baselineModule(datamc, datasetname, flag)]
modules += [simpleCounter(args.report, "TotalEventsPassDataFilters")]
if args.mc: modules += [mcHatModule()]
if args.analyzer == 'main': modules += [analysisModule(), analysisModule(cutbased=True)]
if args.analyzer == 'ztt': modules += [zttModule(tags)]
if args.analyzer == 'trigger': modules += [triggerModule(year="2018", filters_off=False)]
if not args.selection=="None":
  modules += [simpleSelector(args.selection)]
modules += [simpleCounter(args.report, "TotalEventsWritten")]

p = PostProcessor(args.output,
                  files,
                  modules=modules,
                  maxEntries=None,
                  totalEntries=args.numEvents,
                  outputbranchsel=args.branchsel,
                  fwkJobReport=True,
                  haddFileName=args.outfile
                  )
p.run()
