#!/usr/bin/env python3
from __future__ import print_function
import sys
import os
import hashlib
from array import array
import subprocess
import argparse

# constants
convert = {'I' : 'i','F' : 'f',}
in_key =  "---#---#---"
out_key = "===#===#==="

# command line options
parser = argparse.ArgumentParser(description="")
parser.add_argument("report", help=".txt file with metadata report")
parser.add_argument("rootfile", help="file with Events tree")
parser.add_argument("dataset", default="MyDatasetName", help="name of dataset")
parser.add_argument("proc", type=int, help="integer flag")
parser.add_argument("--xs", type=float, default=1.0, help="integer flag")
args = parser.parse_args()
import ROOT

datasetname_id = int(hashlib.sha256(args.dataset.encode('utf-8')).hexdigest(), 16) % 10**8
name_to_branchname = {
'TotalEventsWritten' : 'evtWritten/I',
'TotalEventsProcessed' : 'evtProcessed/I',
'TotalEventsPassDataFilters' : 'evtPassDatafilter/I',
}
prefixes = ['AnaTau_', 'AnaTp_', 'AnaTpm_']
tags = [
'ztt_total',
'ztt_flow_trig',
'ztt_flow_allveto',
'ztt_flow_muAndTau',
'ztt_flow_fullsel',
'ztt_passOnly_muVeto',
'ztt_passOnly_elVeto',
'ztt_passOnly_diMuVeto',
'ztt_passOnly_bVeto',
'ztt_nminusone_mt',
'ztt_nminusone_pzeta',
'ztt_nminusone_eveto',
'ztt_nminusone_muveto',
'ztt_nminusone_dimuveto',
'ztt_nminusone_bjetveto',
]
for prefix in prefixes:
  for tag in tags:
    name_to_branchname[prefix+tag] = 'ztt'+prefix+tag+'/I'

branchname_to_value = {
'xs/F' : 2.0,
'dataset_id/I' : datasetname_id,
'dataset/S' : args.dataset,
}

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

def get_values(reports, names):
    name_to_value = {}
    for report in reports:
        for name in names:
            if report['name'] == name:
                name_to_value[name] = report['count']
    return name_to_value

def add_value(tree, branch, value):
    t = branch[-1]
    branchname = branch[:-2]
    branchstr = branch
    if t == 'S':
      globals()[branch] = ROOT.string('')
      globals()[branch].assign(value)
      tree.Branch(branchname, globals()[branch])
    else:
      globals()[branch] = array(convert[t], [ 0 ])
      globals()[branch][0] = value
      tree.Branch(branchname, globals()[branch], branchstr)

### run ###

reports = parse(args.report)
print("this many reports:", len(reports))
print("adding metadata")
fi = ROOT.TFile(args.rootfile, "update")
metadata = ROOT.TTree('Metadata', 'Metadata')

name_to_value = get_values(reports, name_to_branchname.keys())
for name in name_to_value:
    add_value(metadata, name_to_branchname[name], name_to_value[name])
for branchname in branchname_to_value:
    add_value(metadata, branchname, branchname_to_value[branchname])

metadata.Fill()
metadata.BuildIndex("dataset_id")
metadata.Scan("*")
print("Dump")
for entry in metadata:
  for b in entry.GetListOfBranches():
    print(b.GetName(), entry.__getattr__(b.GetName()))

metadata.Write()
