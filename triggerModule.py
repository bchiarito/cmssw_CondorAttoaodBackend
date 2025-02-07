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
    def __init__(self, year="2018", filters_off=True):
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

        if self.year == 'UL18' or self.year == 'UL17':
          pass_filters = (
          self.mygetattr(flags, 'goodVertices', True)
          and self.mygetattr(flags, 'HBHENoiseFilter', True)
          and self.mygetattr(flags, 'HBHENoiseIsoFilter', True)
          and self.mygetattr(flags, 'EcalDeadCellTriggerPrimitiveFilter', True)
          and self.mygetattr(flags, 'BadPFMuonFilter', True)
          and self.mygetattr(flags, 'BadChargedCandidateFilter', True)
          and self.mygetattr(flags, 'ecalBadCalibFilter', True)
          and self.mygetattr(flags, 'globalSuperTightHalo2016Filter', True)
          and self.mygetattr(flags, 'eeBadScFilter', True)
          )
          if not (pass_filters): return False
          self.countPassingEvents += 1
          return True
        elif self.year == 'UL16':
          pass_filters = (
          self.mygetattr(flags, 'goodVertices', True)
          and self.mygetattr(flags, 'HBHENoiseFilter', True)
          and self.mygetattr(flags, 'HBHENoiseIsoFilter', True)
          and self.mygetattr(flags, 'EcalDeadCellTriggerPrimitiveFilter', True)
          and self.mygetattr(flags, 'BadPFMuonFilter', True)
          and self.mygetattr(flags, 'BadChargedCandidateFilter', True)
          and self.mygetattr(flags, 'globalSuperTightHalo2016Filter', True)
          and self.mygetattr(flags, 'eeBadScFilter', True)
          )
          if not (pass_filters): return False
          self.countPassingEvents += 1
          return True
        else: return True
