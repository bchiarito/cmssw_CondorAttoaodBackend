from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
from math import fabs

const_PhotonCutBasedIDMin = 2
const_TwoProngMaxEta = 2.5
const_PhotonMaxEta = 2.5
const_PairMinDr = 0.1

def deltaR(obj1, obj2):
  phi1 = obj1.phi
  phi2 = obj2.phi
  eta1 = obj1.eta
  eta2 = obj2.eta
  dphi = fabs(phi2-phi1)
  if dphi > math.pi: dphi -= 2*math.pi
  if dphi < math.pi: dphi += 2*math.pi
  deta = fabs(eta2-eta1)
  return math.sqrt(dphi**2 + deta**2)

class recoPhiModule(Module):
    def __init__(self, photon):
        self.photon = photon
        # cutflow vars
        self.total = 0
        self.cutflow_pass_photon = 0
        self.cutflow_pass_twoprong = 0

    def beginJob(self):
        pass

    def endJob(self):
        print "======================================"
        print "Phi Reconstruction Event-based Cutflow"
        print "======================================"
        print "{:<30} {}".format("Total Events", self.total)
        print "{:<30} {}".format(">=1 photon", self.cutflow_pass_photon)
        print "{:<30} {}".format(">=1 twoprong", self.cutflow_pass_twoprong)
        print "======================================\n"

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("HT", "F")
        self.out.branch("NJets", "I")
        if self.photon == "cutBased": self.cutbasedstr = "CutBased"
        if self.photon == "HPID": self.cutbasedstr = ""
        self.out.branch(self.cutbasedstr+"RecoPhi_pass", "B")
        self.out.branch(self.cutbasedstr+"RecoPhi_pt", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_eta", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_phi", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_mass", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_photonLeg_pt", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_photonLeg_eta", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_photonLeg_phi", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_photonLeg_mass", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_twoprongLeg_pt", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_twoprongLeg_eta", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_twoprongLeg_phi", "F")
        self.out.branch(self.cutbasedstr+"RecoPhi_twoprongLeg_mass", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        twoprongs = Collection(event, "TwoProng")
        jets = Collection(event, "Jet")
        if self.photon == "cutBased": photons = Collection(event, "Photon")
        if self.photon == "HPID": photons = Collection(event, "HighPtIdPhoton")

        # cutflow vars
        self.total += 1
        pass_photon = False
        pass_twoprong = False

        # reconstruct phi
        ii = -1
        for i, photon in enumerate(photons):
          pass_id = False
          if self.photon == "HPID": pass_id = True
          if self.photon == "cutBased":
            if photon.cutBased < 2: pass_id = True
          if not pass_id: continue
          ii = i
          break
        if not ii == -1:
          pass_photon = True
          self.cutflow_pass_photon += 1
          photon_cand = photons[ii]
          photon_vec = ROOT.Math.PtEtaPhiMVector(photon_cand.pt, photon_cand.eta, photon_cand.phi, photon_cand.mass)
        jj = -1
        for i, twoprong in enumerate(twoprongs):
          pass_id = False
          if twoprong.pt > 20 and fabs(twoprong.eta)<2.5: pass_id = True
          if not pass_id: continue
          if pass_photon and deltaR(twoprong, photons[ii]) < 0.1: continue # only check if too close to the leading photon
          jj = i
          break
        if not jj == -1:
          pass_twoprong = True
          self.cutflow_pass_twoprong += 1
          twoprong_cand = twoprongs[jj]
          twoprong_vec = ROOT.Math.PtEtaPhiMVector(twoprong_cand.pt, twoprong_cand.eta, twoprong_cand.phi, twoprong_cand.mass)
        if pass_photon and pass_twoprong:
          recophi_vec = twoprong_vec + photon_vec

        # compute event wide quantities
        HT = 0
        if pass_photon: HT += photon_cand.pt
        if pass_twoprong: HT += twoprong_cand.pt
        NJets = 0
        for jet in jets:
          if jet.pt < 30: continue
          if fabs(jet.eta) > 2.4: continue
          if pass_photon and deltaR(jet, photon_cand) < 0.3: continue
          if pass_twoprong and deltaR(jet, twoprong_cand) < 0.3: continue
          HT += jet.pt
          NJets += 1
        
        recophi_pass = pass_twoprong and pass_photon

        # fill branches
        self.out.fillBranch("HT", HT)
        self.out.fillBranch("NJets", NJets)
        self.out.fillBranch(self.cutbasedstr+"RecoPhi_pass", recophi_pass)
        if recophi_pass:
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_pt", recophi_vec.Pt())
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_eta", recophi_vec.Eta())
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_phi", recophi_vec.Phi())
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_mass", recophi_vec.M())
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_photonLeg_pt", photon_vec.Pt())
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_photonLeg_eta", photon_vec.Eta())
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_photonLeg_phi", photon_vec.Phi())
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_photonLeg_mass", photon_vec.M())
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_twoprongLeg_pt", twoprong_vec.Pt())
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_twoprongLeg_eta", twoprong_vec.Eta())
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_twoprongLeg_phi", twoprong_vec.Phi())
          self.out.fillBranch(self.cutbasedstr+"RecoPhi_twoprongLeg_mass", twoprong_vec.M())
        return True

recoPhiConstr_cutBased = lambda: recoPhiModule(photon='cutBased')
recoPhiConstr_HPID = lambda: recoPhiModule(photon='HPID')
