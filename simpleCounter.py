from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
import math

class simpleCounter(Module):
    def __init__(self, name):
        self.total = 0
        self.name = name

    def beginJob(self):
        pass

    def endJob(self):
        print("===#===#===")
        print self.name
        print "count", self.total
        print "---#---#---\n"

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        self.total += 1
        return True
