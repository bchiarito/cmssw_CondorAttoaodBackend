from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import os
import sys
import ROOT
from ROOT import TMath
ROOT.PyConfig.IgnoreCommandLineOptions = True
import array
import math
import hashlib
from math import fabs


class triggerModule(Module):
    def __init__(self, year="2018", filters_off=False):
        self.year = year
        self.filters_off = filters_off
        pass     

    def beginJob(self):
        pass
                      
    def endJob(self):
        pass
                             
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
 
    def mygetattr(self, my_obj, my_branch, default_bool):
        try:
            return getattr(my_obj, my_branch)
        except RuntimeError:
            return default_bool

    def analyze(self, event):
        return True

