#!/usr/bin/env python3
from __future__ import print_function
import sys
import os
import hashlib
from array import array
import subprocess
import argparse
import ROOT

# command line options
parser = argparse.ArgumentParser(description="", usage="./%(prog)s")
parser.add_argument("report", help=".txt file with metadata report")
parser.add_argument("rootfile", help="file with Events tree")
parser.add_argument("dataset", default="MyDatasetName", help="name of dataset")
parser.add_argument("proc", type=int, help="integer flag")
parser.add_argument("--xs", type=float, default=1.0, help="integer flag")
args = parser.parse_args()

report_filename = args.report
tree_filename = args.rootfile
arg_flag = int(args.proc)
datasetname = args.dataset
datasetname_id = int(hashlib.sha256(datasetname.encode('utf-8')).hexdigest(), 16) % 10**8

in_key =  "---#---#---"
out_key = "===#===#==="

def parse(filename):
  fi = open(filename)
  collection = []
  i = False
  for line in fi:
    line = line.strip()
    #print(line)
    if line == in_key:
      d = {}
      i = True
    if line == out_key:
      collection.append(d)
      i = False
    if i:
      line = line.split()
      #print(line)
      if len(line) == 1: d['name'] = line[0]
      if len(line) == 2:
        num = line[1]
        try: num = int(num)
        except ValueError:
          num = 0
          print("got a ValueError! on line", line)
        d[line[0]] = num
  return collection

reports = parse(report_filename)

print("this many reports:", len(reports))
for i, d in enumerate(reports):
  print("  report", i+1)
  for key in d:
    print(key, d[key])

for report in reports:
  if report['name'] == 'TotalEventsWritten':
    report_evtWritten = report['count']
  if report['name'] == 'TotalEventsProcessed':
    report_evtProcessed = report['count']

print("adding metadata")
fi = ROOT.TFile(tree_filename, "update")
metadata = ROOT.TTree('Metadata', 'Metadata')

dataset = ROOT.string('')
dataset_id = array('i', [ 0 ])
flag = array('i', [ 0 ])
evtWritten = array('i', [ 0 ])
evtProcessed = array('i', [ 0 ])
xs = array('f', [ 0.0 ])

metadata.Branch('dataset', dataset)
metadata.Branch('dataset_id', dataset_id, 'dataset_id/I')
metadata.Branch('flag', flag, 'flag/I')
metadata.Branch('evtWritten', evtWritten, 'evtWritten/I')
metadata.Branch('evtProcessed', evtProcessed, 'evtProcessed/I')
metadata.Branch('xs', xs, 'xs/F')

dataset.assign(datasetname)
dataset_id[0] = datasetname_id
flag[0] = arg_flag
evtWritten[0] = report_evtWritten
evtProcessed[0] = report_evtProcessed
xs[0] = args.xs

metadata.Fill()
metadata.BuildIndex("dataset_id", "flag")
#metadata.Print()
metadata.Write()
