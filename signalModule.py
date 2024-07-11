from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math

class signalModule(Module):
    def __init__(self):
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("weight_HEM", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        genomegas = Collection(event, "GenOmega")

        weight_HEM = 1
        for genomega in genomegas:
          if genomega.phi > 1: weight_HEM = 2
          else: weight_HEM = 1

        # fill branches
        self.out.fillBranch("weight_HEM", weight_HEM)
        return True

