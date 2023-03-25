from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import hashlib
from math import fabs

PHOTON_CUTBASED_ID = 0 # loose

def get_vec(obj):
  return ROOT.Math.PtEtaPhiMVector(obj.pt, obj.eta, obj.phi, obj.mass)

class analysisModule(Module):
    def __init__(self, selection_type="", cutbased=False):
        self.selection_type = selection_type # for future, high pt / low pt selection
        self.cutbased = cutbased
        if cutbased: self.cutbasedstr = "CutBased_"
        else: self.cutbasedstr = ""

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.cutbasedstr+"HT", "F")
        self.out.branch(self.cutbasedstr+"NJets", "I")
        self.out.branch(self.cutbasedstr+"Region", "I")
        self.out.branch(self.cutbasedstr+"RecoPhi_pt", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_eta", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_phi", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_mass", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_photonindex", "I")
        self.out.branch(self.cutbasedstr+"RecoPhi_twoprongindex", "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        twoprongs = Collection(event, "TwoProng")
        jets = Collection(event, "Jet")
        if not self.cutbased: photons = Collection(event, "HighPtIdPhoton")
        else: photons = Collection(event, "Photon")

        # FIXME sorting should be done upstream
        twoprongs = sorted(twoprongs, reverse=True, key=lambda obj : obj.pt)

        # recophi
        pass_photon = False
        pass_twoprong = ''
        photon_index = -1
        twoprong_index = -1
        if self.cutbased:
          for i, photon in enumerate(photons):
            if photon.cutBased <= PHOTON_CUTBASED_ID: continue
            pass_photon = True
            photon_index = i
            cutbased_photon_index = i
            photon_vec = ROOT.Math.PtEtaPhiMVector(photon.pt, photon.eta, photon.phi, photon.mass)
            break
        else:
          if len(photons) >= 1:
            pass_photon = True
            photon_index = 0
            photon_vec = get_vec(photons[0])
        if pass_photon:
          # look for tight
          for i, twoprong in enumerate(twoprongs):
            try:
              tight = twoprong.isTight
            except RuntimeError:
              tight = True
            if not (tight and twoprong.pt > 20 and fabs(twoprong.eta)<2.5): continue
            if ROOT.Math.VectorUtil.DeltaR(photon_vec, get_vec(twoprong)) < 0.1: continue
            twoprong_index = i
            pass_twoprong = 'tight'
            break
          # look for loose
          if twoprong_index == -1:
            for i, twoprong in enumerate(twoprongs):
              try:
                tight = twoprong.isTight
              except RuntimeError:
                tight = True
              if not (not tight and twoprong.pt > 20 and fabs(twoprong.eta)<2.5): continue
              if ROOT.Math.VectorUtil.DeltaR(photon_vec, get_vec(twoprong)) < 0.1: continue
              twoprong_index = i
              pass_twoprong = 'loose'
              break
        if pass_photon and pass_twoprong == 'tight':
          region = 1
          recophi_vec = get_vec(twoprongs[twoprong_index]) + photon_vec
        elif pass_photon and pass_twoprong == 'loose':
          region = 2
          recophi_vec = get_vec(twoprongs[twoprong_index]) + photon_vec
        else:
          region = 0
          recophi_vec = ROOT.Math.PtEtaPhiMVector(0, 0, 0, 0)

        # event wide quantities
        HT = 0
        NJets = 0
        if pass_photon: photon_vec = get_vec(photons[photon_index])
        if pass_twoprong: twoprong_vec = get_vec(twoprongs[twoprong_index])
        for jet in jets:
          if jet.pt < 30: continue
          if fabs(jet.eta) > 2.4: continue
          jet_vec = get_vec(jet)
          if pass_photon and ROOT.Math.VectorUtil.DeltaR(photon_vec, jet_vec) < 0.3: continue
          if pass_twoprong and ROOT.Math.VectorUtil.DeltaR(twoprong_vec, jet_vec) < 0.3: continue
          HT += jet.pt
          NJets += 1
        if pass_photon: HT += photons[photon_index].pt
        if pass_twoprong: HT += twoprongs[twoprong_index].pt

        # fill branches
        self.out.fillBranch(self.cutbasedstr+"HT", HT)
        self.out.fillBranch(self.cutbasedstr+"NJets", NJets)
        self.out.fillBranch(self.cutbasedstr+"Region", region)
        self.out.fillBranch(self.cutbasedstr+"RecoPhi_pt", recophi_vec.Pt())
        self.out.fillBranch(self.cutbasedstr+"RecoPhi_eta", recophi_vec.Eta())
        self.out.fillBranch(self.cutbasedstr+"RecoPhi_phi", recophi_vec.Phi())
        self.out.fillBranch(self.cutbasedstr+"RecoPhi_mass", recophi_vec.M())
        self.out.fillBranch(self.cutbasedstr+"RecoPhi_photonindex", photon_index)
        self.out.fillBranch(self.cutbasedstr+"RecoPhi_twoprongindex", twoprong_index)
        return True
