from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import sys
import ROOT
import math

class simpleCounter(Module):
    def __init__(self, filename, name):
        self.total = 0
        self.filename = filename
        self.name = name

    def beginJob(self):
        pass

    def endJob(self):
        with open(self.filename, "a") as f:
          f.write("---#---#---\n")
          f.write(self.name+"\n")
          f.write("count "+ str(self.total)+"\n")
          f.write("===#===#===\n\n")

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        self.total += 1
        return True
