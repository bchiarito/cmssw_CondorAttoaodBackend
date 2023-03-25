from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
import math

class simpleSelector(Module):
    def __init__(self, selection=''):
        self.sel = selection
        self.total = 0
        self.passed = 0

    def beginJob(self):
        pass

    def endJob(self):
        print "====================="
        print "simpleSelector Report"
        print self.sel
        print "====================="
        print "{:<30} {}".format("Total Events", self.total)
        print "{:<30} {}".format("Passed", self.passed)
        print "====================="

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")
        photons = Collection(event, "Photon")
        id_photons = Collection(event, "HighPtIdPhoton")
        twoprongs = Collection(event, "TwoProng")
        recophi = Object(event, "RecoPhi")
        cb_recophi = Object(event, "CutBased_RecoPhi")
        pass_trigger = event.HLT_Photon200
        self.total += 1
        if self.sel == 'one_muon':
          muons = Collection(event, "Muon")
          if len(muons) == 0 : return False
        if self.sel == 'one_photon':
          photons = Collection(event, "Photon")
          if len(photons) == 0 : return False
        if self.sel == 'one_hpid_photon':
          id_photons = Collection(event, "HighPtIdPhoton")
          if len(id_photons) == 0 : return False
        if self.sel == 'hpid_photon_ptcut':
          id_photons = Collection(event, "HighPtIdPhoton")
          if len(id_photons) == 0 or id_photons[0].pt<30: return False
        if self.sel == 'bkgEst_allRegions':
          if len(id_photons) == 0 : return False
          if len(twoprongs) == 0 : return False
        if self.sel == 'full_selection_hpid':
          if not event.Region == 1: return False
          if not pass_trigger: return False
          if not (id_photons[recophi.photonindex].pt > 220) : return False
        if self.sel == 'full_selection_cbl':
          if not (event.CutBased_Region == 1): return False
          if not pass_trigger: return False
          if not (photons[cb_recophi.photonindex].pt > 220) : return False
        self.passed += 1
        return True
