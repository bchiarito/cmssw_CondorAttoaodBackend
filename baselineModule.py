from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import hashlib
from math import fabs

def deltaR(obj1, obj2):
  phi1 = obj1.phi
  phi2 = obj2.phi
  eta1 = obj1.eta
  eta2 = obj2.eta
  dphi = fabs(phi2-phi1)
  if dphi > math.pi: dphi -= 2*math.pi
  if dphi < -math.pi: dphi += 2*math.pi
  deta = fabs(eta2-eta1)
  return math.sqrt(dphi**2 + deta**2)

class baselineModule(Module):
    def __init__(self, datamc, dataset='', flag=0):
        self.dataset = dataset
        self.flag = flag
        self.datamc = datamc

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("dataset_id", "I")
        self.out.branch("flag", "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets = Collection(event, "Jet")
        flags = Object(event, "Flag")

        # fill branches
        dataset_id = int(hashlib.sha256((self.dataset).encode('utf-8')).hexdigest(), 16) % 10**8
        flag = self.flag
        self.out.fillBranch("dataset_id", dataset_id)
        self.out.fillBranch("flag", flag)

        # filter descion for data and signal mc
        pass_filters = False
        if self.datamc == 'data' or self.datamc == 'sigRes' or self.datamc == 'sigNonRes':
          if flags.goodVertices and \
            flags.globalSuperTightHalo2016Filter and \
            flags.HBHENoiseFilter and \
            flags.HBHENoiseIsoFilter and \
            flags.EcalDeadCellTriggerPrimitiveFilter and \
            flags.BadPFMuonFilter and \
            flags.eeBadScFilter and \
            flags.ecalBadCalibFilter:
              pass_filters = True
          return pass_filters
        else:
          return True

