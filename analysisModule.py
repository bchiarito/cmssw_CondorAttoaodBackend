from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import hashlib
from math import fabs

PHOTON_CUTBASED_ID = 1 # loose
PHOTON_BARREL_ETA = 1.4442
PHOTON_MIN_PT = 35 # b/c HLT_Photon35_Twoprongs35
PHOTON_HoverE_CUT = 0.04596

MIN_DR = 0.1 # between photon and twoprong

TWOPRONG_MIN_PT = 20

def get_vec(obj, corr=""):
  if not corr=="": et = obj.pt * ROOT.TMath.CosH(obj.eta)
  if corr == "scale_up":
      et_scale_up = et + abs(obj.dEscaleUp)
      pt_scale_up = et_scale_up / ROOT.TMath.CosH(obj.eta)
      return ROOT.Math.PtEtaPhiMVector(pt_scale_up, obj.eta, obj.phi, obj.mass)
  elif corr == "scale_down":
      et_scale_down = et - abs(obj.dEscaleDown)
      pt_scale_down = et_scale_down / ROOT.TMath.CosH(obj.eta)
      return ROOT.Math.PtEtaPhiMVector(pt_scale_down, obj.eta, obj.phi, obj.mass)
  elif corr == "smear_up":
      et_sigma_up = et + abs(obj.dEsigmaUp)
      pt_sigma_up = et_sigma_up / ROOT.TMath.CosH(obj.eta)
      return ROOT.Math.PtEtaPhiMVector(pt_sigma_up, obj.eta, obj.phi, obj.mass)
  elif corr == "smear_down":
      et_sigma_down = et - abs(obj.dEsigmaDown)
      pt_sigma_down = et_sigma_down / ROOT.TMath.CosH(obj.eta)
      return ROOT.Math.PtEtaPhiMVector(pt_sigma_down, obj.eta, obj.phi, obj.mass)
  else: return ROOT.Math.PtEtaPhiMVector(obj.pt, obj.eta, obj.phi, obj.mass)

class analysisModule(Module):
    def __init__(self, selection_type="", cutbased=True):
        self.selection_type = selection_type # for future, high pt / low pt selection
        self.cutbased = cutbased
        if cutbased: self.photon_type = "CBL_"
        else: self.photon_type = "HPID_"

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.photon_type+"HT", "F")
        self.out.branch(self.photon_type+"NJets", "I")
        self.out.branch(self.photon_type+"Region", "I")
        self.out.branch(self.photon_type+"RecoPhi_pt", "F")
        self.out.branch(self.photon_type+"RecoPhi_eta", "F")
        self.out.branch(self.photon_type+"RecoPhi_phi", "F")
        self.out.branch(self.photon_type+"RecoPhi_mass", "F")
        self.out.branch(self.photon_type+"RecoPhi_photonindex", "I")
        self.out.branch(self.photon_type+"RecoPhi_twoprongindex", "I")
        self.out.branch(self.photon_type+"RecoPhi_TwoProng_pt", "F")
        self.out.branch(self.photon_type+"RecoPhi_TwoProng_eta", "F")
        self.out.branch(self.photon_type+"RecoPhi_TwoProng_phi", "F")
        self.out.branch(self.photon_type+"RecoPhi_TwoProng_mass", "F")
        self.out.branch(self.photon_type+"RecoPhi_TwoProng_massPi0", "F")
        self.out.branch(self.photon_type+"RecoPhi_TwoProng_massEta", "F")
        self.out.branch(self.photon_type+"RecoPhi_Photon_pt", "F")
        self.out.branch(self.photon_type+"RecoPhi_Photon_eta", "F")
        self.out.branch(self.photon_type+"RecoPhi_Photon_phi", "F")
        self.out.branch(self.photon_type+"RecoPhi_Photon_mass", "F")
        self.out.branch(self.photon_type+"RecoPhi_dr", "F")
        # Add scale/smear recophi mass branches
        self.out.branch(self.photon_type+"RecoPhi_mass_scale_up", "F")
        self.out.branch(self.photon_type+"RecoPhi_mass_scale_down", "F")
        self.out.branch(self.photon_type+"RecoPhi_mass_smear_up", "F")
        self.out.branch(self.photon_type+"RecoPhi_mass_smear_down", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        twoprongs = Collection(event, "TwoProng")
        jets = Collection(event, "Jet")
        if not self.cutbased: photons = Collection(event, "HighPtIdPhoton")
        else: photons = Collection(event, "Photon")

        photon_corr = True
        try:
            if not len(photons) == 0: test_var = photons[0].dEscaleUp
        except:
            photon_corr = False

        # FIXME sorting should be done upstream
        twoprongs = sorted(twoprongs, reverse=True, key=lambda obj : obj.pt)

        # find photon
        pass_photon = False
        photon_region = ''
        photon_index = -1
        if self.cutbased:
          # priority for tight
          for i, photon in enumerate(photons):
            if photon.cutBased < PHOTON_CUTBASED_ID: continue
            if abs(photon.scEta) > PHOTON_BARREL_ETA: continue
            if photon.pt <= PHOTON_MIN_PT: continue
            if photon.hadTowOverEm > PHOTON_HoverE_CUT: continue
            photon_region = 'tight'
            photon_index = i
            pass_photon = True
            break
          # try to find loose
          if photon_index == -1:
            for i, photon in enumerate(photons):
              if abs(photon.scEta) > PHOTON_BARREL_ETA: continue
              if photon.pt <= PHOTON_MIN_PT: continue
              if photon.hadTowOverEm > PHOTON_HoverE_CUT: continue
              photon_region = 'loose'
              photon_index = i
              pass_photon = True
              break
        else:
          for i, photon in enumerate(photons):
            if abs(photon.scEta) > PHOTON_BARREL_ETA: continue
            if photon.pt <= PHOTON_MIN_PT: continue
            photon_region = 'tight'
            photon_index = i
            pass_photon = True
            break
        
        # find twoprong
        pass_twoprong = False
        twoprong_region = ''
        twoprong_index = -1
        TT = []
        TL = []
        LT = []
        LL = []
        if not photon_index == -1:
          for i, twoprong in enumerate(twoprongs):
            if twoprong.pt < TWOPRONG_MIN_PT: continue
            if ROOT.Math.VectorUtil.DeltaR(get_vec(twoprong), get_vec(photons[photon_index])) < MIN_DR: continue
            try:
              if twoprong.isTight: TT.append((i, twoprong))
              if twoprong.passIso and not twoprong.passSym: TL.append((i, twoprong))
              if not twoprong.passIso and twoprong.passSym: LT.append((i, twoprong))
              if not twoprong.passIso and not twoprong.passSym: LL.append((i, twoprong))
            except RuntimeError:
              TT.append((i, twoprong))
        if len(TT) > 0:
          pass_twoprong = True
          twoprong_region = 'tight'
          twoprong_index = TT[0][0]
        elif len(TL) > 0:
          pass_twoprong = True
          twoprong_region = 'iso_asym'
          twoprong_index = TL[0][0]
        elif len(LT) > 0:
          pass_twoprong = True
          twoprong_region = 'noniso_sym'
          twoprong_index = LT[0][0]
        elif len(LL) > 0:
          pass_twoprong = True
          twoprong_region = 'noniso_asym'
          twoprong_index = LL[0][0]

        if photon_region == '' or twoprong_region == '': region = 0
        elif photon_region == 'tight' and twoprong_region == 'tight': region = 1
        elif photon_region == 'tight' and twoprong_region == 'iso_asym': region = 2
        elif photon_region == 'tight' and twoprong_region == 'noniso_sym': region = 3
        elif photon_region == 'tight' and twoprong_region == 'noniso_asym': region = 4
        elif photon_region == 'loose' and twoprong_region == 'tight': region = 11
        elif photon_region == 'loose' and twoprong_region == 'iso_asym': region = 12
        elif photon_region == 'loose' and twoprong_region == 'noniso_sym': region = 13
        elif photon_region == 'loose' and twoprong_region == 'noniso_asym': region = 14


        if region > 0:
          recophi_vec = get_vec(twoprongs[twoprong_index]) + get_vec(photons[photon_index])  # uncorrected photon pt
          if photon_corr:  # Apply scaling and resolution corrections
              recophi_vec_scale_up = get_vec(twoprongs[twoprong_index]) + get_vec(photons[photon_index], "scale_up")
              recophi_vec_scale_down = get_vec(twoprongs[twoprong_index]) + get_vec(photons[photon_index], "scale_down")
              recophi_vec_smear_up = get_vec(twoprongs[twoprong_index]) + get_vec(photons[photon_index], "smear_up")
              recophi_vec_smear_down = get_vec(twoprongs[twoprong_index]) + get_vec(photons[photon_index], "smear_down")
        else:
          recophi_vec = ROOT.Math.PtEtaPhiMVector(0, 0, 0, 0)
          if photon_corr: 
              recophi_vec_scale_up = ROOT.Math.PtEtaPhiMVector(0, 0, 0, 0)
              recophi_vec_scale_down = ROOT.Math.PtEtaPhiMVector(0, 0, 0, 0)
              recophi_vec_smear_up = ROOT.Math.PtEtaPhiMVector(0, 0, 0, 0)
              recophi_vec_smear_down = ROOT.Math.PtEtaPhiMVector(0, 0, 0, 0)
        
        # compute dr
        if region == 1:
            pi = ROOT.TMath.Pi()
            dEta = abs(photons[photon_index].eta - twoprongs[twoprong_index].eta)
            dPhi = photons[photon_index].phi - twoprongs[twoprong_index].phi
            if dPhi > pi: dPhi -= 2.0*pi
            elif dPhi <= -pi: dPhi += 2.0*pi
            dr = ROOT.TMath.Sqrt(dEta**2 + dPhi**2)

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
        self.out.fillBranch(self.photon_type+"HT", HT)
        self.out.fillBranch(self.photon_type+"NJets", NJets)
        self.out.fillBranch(self.photon_type+"Region", region)
        self.out.fillBranch(self.photon_type+"RecoPhi_pt", recophi_vec.Pt())
        self.out.fillBranch(self.photon_type+"RecoPhi_eta", recophi_vec.Eta())
        self.out.fillBranch(self.photon_type+"RecoPhi_phi", recophi_vec.Phi())
        self.out.fillBranch(self.photon_type+"RecoPhi_mass", recophi_vec.M())
        self.out.fillBranch(self.photon_type+"RecoPhi_photonindex", photon_index)
        self.out.fillBranch(self.photon_type+"RecoPhi_twoprongindex", twoprong_index)
        if photon_corr:
            # scale and smear corrections
            self.out.fillBranch(self.photon_type+"RecoPhi_mass_scale_up", recophi_vec_scale_up.M())
            self.out.fillBranch(self.photon_type+"RecoPhi_mass_scale_down", recophi_vec_scale_down.M())
            self.out.fillBranch(self.photon_type+"RecoPhi_mass_smear_up", recophi_vec_smear_up.M())
            self.out.fillBranch(self.photon_type+"RecoPhi_mass_smear_down", recophi_vec_smear_down.M())
        if not twoprong_index == -1:
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_pt", twoprongs[twoprong_index].pt)
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_eta", twoprongs[twoprong_index].eta)
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_phi", twoprongs[twoprong_index].phi)
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_mass", twoprongs[twoprong_index].mass)
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_massPi0", twoprongs[twoprong_index].massPi0)
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_massEta", twoprongs[twoprong_index].massEta)
        else:
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_pt", -1)
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_eta", -1)
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_phi", -1)
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_mass", -1)
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_massPi0", -1)
          self.out.fillBranch(self.photon_type+"RecoPhi_TwoProng_massEta", -1)
        if not photon_index == -1:
          self.out.fillBranch(self.photon_type+"RecoPhi_Photon_pt", photons[photon_index].pt)
          self.out.fillBranch(self.photon_type+"RecoPhi_Photon_eta", photons[photon_index].eta)
          self.out.fillBranch(self.photon_type+"RecoPhi_Photon_phi", photons[photon_index].phi)
          self.out.fillBranch(self.photon_type+"RecoPhi_Photon_mass", photons[photon_index].mass)
        else:
          self.out.fillBranch(self.photon_type+"RecoPhi_Photon_pt", -1)
          self.out.fillBranch(self.photon_type+"RecoPhi_Photon_eta", -1)
          self.out.fillBranch(self.photon_type+"RecoPhi_Photon_phi", -1)
          self.out.fillBranch(self.photon_type+"RecoPhi_Photon_mass", -1)
        if region == 1:
            if photon_index == -1 or twoprong_index == -1:
              self.out.fillBranch(self.photon_type+"RecoPhi_dr", -1)
            else:
              self.out.fillBranch(self.photon_type+"RecoPhi_dr", dr)

          
        return True
