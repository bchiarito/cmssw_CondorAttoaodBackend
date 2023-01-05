from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math

class mcHatModule(Module):
    def __init__(self):
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("htHat_lhe", "F")
        self.out.branch("htHat_genPart", "F")
        self.out.branch("htHat_genJet", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        genparticles = Collection(event, "GenPart")
        genjets = Collection(event, "GenJet")
        lheparticles = Collection(event, "LHEPart")

        htHat_genJet = 0
        for genjet in genjets:
          if genjet.pt < 30: continue
          if math.fabs(genjet.eta) > 2.5: continue
          htHat_genJet += genjet.pt

        htHat_genPart = 0
        for genparticle in genparticles:
          if genparticle.status < 20 or genparticle.status > 29: continue
          i = abs(genparticle.pdgId)
          if i == 1 or i == 2 or i == 3 or i == 4 or i == 5 or i == 6 or i == 21:
            htHat_genPart += genparticle.pt

        htHat_lhe = 0
        for lheparticle in lheparticles:
          i = abs(lheparticle.pdgId)
          if i == 1 or i == 2 or i == 3 or i == 4 or i == 5 or i == 6 or i == 21:
            htHat_lhe += lheparticle.pt

        # fill branches
        self.out.fillBranch("htHat_genJet", htHat_genJet)
        self.out.fillBranch("htHat_genPart", htHat_genPart)
        self.out.fillBranch("htHat_lhe", htHat_lhe)
        return True

#recoPhiConstr_HPID = lambda: recoPhiModule(photon='HPID')
