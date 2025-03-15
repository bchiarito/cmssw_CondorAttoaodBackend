from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np
import os
import ctypes
import json

CUT_ZPAIR_DR = 0.5
CUT_PZETA = -25
CUT_MT = 30

class zttModule(Module):
    def __init__(self, report_filename, tags, datamc):
      self.report = report_filename
      self.metadata = {}
      self.m4_6func= None
      self.centralfunc= None
      self.p4_6func= None
      self.kSpreadMC = None
      self.kSmearMC = None
      self.kScaleDT = None
      self.roccin = None
      self.DMC = None
      self.yeartag= None
      self.mcnum= None
      self.tags = tags
      self.this_record = None
      self.datasetscalefactor= None
      self.nevents=None
      self.antiiso= None
      self.datamc = datamc

      prefixes = ['AnaTau_', 'AnaTp_', 'AnaTpm_']
      for prefix in prefixes:
        self.metadata[prefix+'ztt_total'] = 0
        self.metadata[prefix+'ztt_flow_trig'] = 0
        self.metadata[prefix+'ztt_flow_allveto'] = 0
        self.metadata[prefix+'ztt_flow_muAndTau'] = 0
        self.metadata[prefix+'ztt_flow_fullsel'] = 0
        self.metadata[prefix+'ztt_passOnly_muVeto'] = 0
        self.metadata[prefix+'ztt_passOnly_elVeto'] = 0
        self.metadata[prefix+'ztt_passOnly_diMuVeto'] = 0
        self.metadata[prefix+'ztt_passOnly_bVeto'] = 0
        self.metadata[prefix+'ztt_nminusone_mt'] = 0
        self.metadata[prefix+'ztt_nminusone_pzeta'] = 0
        self.metadata[prefix+'ztt_nminusone_eveto'] = 0
        self.metadata[prefix+'ztt_nminusone_muveto'] = 0
        self.metadata[prefix+'ztt_nminusone_dimuveto'] = 0
        self.metadata[prefix+'ztt_nminusone_bjetveto'] = 0

    def beginJob(self):
	self.DMC = self.tags[0]
	self.yeartag= self.tags[1]
	self.mcnum= self.tags[2]
	self.antiiso=self.tags[3]
	self.nevents=0


	if self.DMC == 1: ##if is MC
		with open('/cms/smd376/CMSSW_11_1_0/src/Ztt.json', 'r') as json_file:
			loaded_data = json.load(json_file)
		for record in loaded_data:
			if record["Number"] == self.mcnum:
				self.this_record = record
				break
		if self.this_record:
    			print("Dataset Located")
    			print("Number: ", self.this_record["Number"])
    			print("Name: ", self.this_record["Name"])
    			print("Cross Section: ", self.this_record["Cross Section"])
			print("Luminosity: ",self.this_record["Luminosity"])
			print("Generated Events: ",self.this_record["Generated Events"])
    			print("Pileup File: ", self.this_record["Pile Up File Name"])
        		self.datasetscalefactor=(self.this_record["Cross Section"]*self.this_record["Luminosity"])/self.this_record["Generated Events"]
			print("Dataset Scale Factor Calculated: ", self.datasetscalefactor)
		else:
			print("No MC Dataset with number {} found. Run will not have proper scale factors.".format(self.mcnum))
			self.datasetscalefactor=1
	else: self.datasetscalefactor = 1
	rc = ctypes.CDLL('/cms/smd376/CMSSW_11_1_0/src/RocWrapper.so')
	print(rc)
	roccorfile = "/cms/smd376/CMSSW_11_1_0/src/PhysicsTools/NanoAODTools/python/postprocessing/data/roccor.Run2.v3/RoccoR2018.txt"
	RoccoR_instance = rc.RoccoRinit
	RoccoR_instance.restype = ctypes.c_void_p
	Roccin = RoccoR_instance(roccorfile)
	self.roccin = Roccin
	print("Instance of RoccoR created")
	rc.kScaleDT.argtypes = [ctypes.c_void_p,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_int]

	rc.kScaleDT.restype = ctypes.c_double

	rc.kSpreadMC.argtypes = [ctypes.c_void_p,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_int]

	rc.kSpreadMC.restype = ctypes.c_double

	rc.kSmearMC.argtypes = [ctypes.c_void_p,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_double,ctypes.c_int,ctypes.c_int]

	rc.kSmearMC.restype = ctypes.c_double
	self.kSpreadMC = rc.kSpreadMC
	self.kSmearMC = rc.kSmearMC
	self.kScaleDT = rc.kScaleDT

	if self.DMC ==1: ##if is MC
		pileup_location = '/cms/smd376/CMSSW_11_1_0/src/PileupFiles/'
		pileup_filename = self.this_record["Pile Up File Name"]
		pileup_path = pileup_location + pileup_filename
		pileup_file=ROOT.TFile(pileup_path)
		tree_name = "Events"
		Truepileup=ROOT.TFile('/cms/smd376/CMSSW_11_1_0/src/TruePileupHistogram.root')
		m4_6Truepileup=ROOT.TFile('/cms/smd376/CMSSW_11_1_0/src/Minus_TruePileupHistogram_66000.root')
		p4_6Truepileup=ROOT.TFile('/cms/smd376/CMSSW_11_1_0/src/Plus_TruePileupHistogram_72400.root')

		centralTrue=Truepileup.Get('pileup')
		centralTrue.SetDirectory(0)

		m4_6True=m4_6Truepileup.Get('pileup')
		m4_6True.SetDirectory(0)

		p4_6True=p4_6Truepileup.Get('pileup')
		p4_6True.SetDirectory(0)

		inputpileup = pileup_file.Get('hist_sum_1')
		inputpileup.SetDirectory(0)


		InthT=centralTrue.Integral(0,100)
		Intm4_6T=m4_6True.Integral(0,100)
		Intp4_6T=p4_6True.Integral(0,100)

		centralTrue.Scale(1/InthT)
		m4_6True.Scale(1/Intm4_6T)
		p4_6True.Scale(1/Intp4_6T)
		intmc=inputpileup.Integral(0,100)
		inputpileup.Scale(1/intmc)


		p4_6Truepileup.Close()
		m4_6Truepileup.Close()
		Truepileup.Close()
		pileup_file.Close()

		self.m4_6func=ROOT.TH1D('minus4_6','-4.6% Scale Factor',100,0,100)
		self.centralfunc=ROOT.TH1D('central','Central Scale Factor',100,0,100)
		self.p4_6func=ROOT.TH1D('plus4_6','4.6% Scale Factor',100,0,100)

		for n in range(0,100):
			try: self.m4_6func.Fill(n,m4_6True.GetBinContent(n+1)/inputpileup.GetBinContent(n+1))
			except: pass
		for n in range(0,100):
        		try: self.centralfunc.Fill(n,centralTrue.GetBinContent(n+1)/inputpileup.GetBinContent(n+1))
        	        except: pass
		for n in range(0,100):
        		try: self.p4_6func.Fill(n,p4_6True.GetBinContent(n+1)/inputpileup.GetBinContent(n+1))
        	 	except: pass
		print("Pileup Calculation Done")

    def endJob(self):
        with open(self.report, "a") as f:
          for key in self.metadata:
            tag = key
            val = self.metadata[key]
            f.write("---#---#---\n")
            f.write(tag+"\n")
            f.write("count "+ str(val)+"\n")
            f.write("===#===#===\n\n")

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("GnTwoProng", "I")
        self.out.branch("GTwoProng_pt", "F")
        self.out.branch("GTwoProng_eta", "F")
        self.out.branch("GTwoProng_phi", "F")
        self.out.branch("GTwoProng_mass", "F")
        self.out.branch("GTwoProng_massl", "F")
        self.out.branch("GTwoProng_massPi0", "F")
        self.out.branch("GTwoProng_massEta", "F")

 	self.out.branch("GnMuon", "I")
        self.out.branch("GMuon_pt", "F")
        self.out.branch("GMuon_eta", "F")
        self.out.branch("GMuon_phi", "F")
        self.out.branch("GMuon_mass", "F")
	self.out.branch("GMuon_charge", "F")
	self.out.branch("GMuon_genpt","F")
	self.out.branch("GMuon_ntrackerlayers","I")
	self.out.branch("GMuon_pfisoid","D")

        self.out.branch("GnTau", "I")
        self.out.branch("GTau_pt", "F")
        self.out.branch("GTau_eta", "F")
        self.out.branch("GTau_phi", "F")
        self.out.branch("GTau_mass", "F")
        self.out.branch("GTau_charge","F")
	self.out.branch("GNPV","D")
	self.out.branch("GnElectron", "I")

	self.out.branch("GEventsRunOver","I")

	self.out.branch("GChpospt","D")
	self.out.branch("GChnegpt","D")
	self.out.branch("GChextra_charge","D")

	self.out.branch("GTPnTracks","D")
	self.out.branch("GMET","D")
	self.out.branch("GPileup_nTrueInt","D")

	self.out.branch("nMuons","D")
	self.out.branch("nTaus","D")
	self.out.branch("nJets","D")

	self.out.branch("GMuon_pt_tpcheck","D")

	self.out.branch("GJet_pt","D")
	self.out.branch("GJet_eta","D")
	self.out.branch("GJet_phi","D")
	self.out.branch("GJet_mass","D")

        ###

	self.out.branch("GPZeta_Tau","D")
	self.out.branch("GPZeta_TP","D")
	self.out.branch("Z_massTau","D")
    	self.out.branch("Z_ptTau","D")
    	self.out.branch("Z_phiTau","D")
    	self.out.branch("Z_etaTau","D")
	self.out.branch("Z_massTP","D")
    	self.out.branch("Z_ptTP","D")
    	self.out.branch("Z_phiTP","D")
    	self.out.branch("Z_etaTP","D")
	self.out.branch("TPTauDeltaR","D")
	self.out.branch("GnJets_corrected","D")
    	self.out.branch("Ght","D")
	self.out.branch("GTransMass","D")

	self.out.branch("GSF_central","D")
	self.out.branch("GSF_plus","D")
	self.out.branch("GSF_minus","D")
	self.out.branch("GSF_pileup_central","D")
        self.out.branch("GSF_pileup_plus","D")
        self.out.branch("GSF_pileup_minus","D")
        self.out.branch("GSF_id_iso_trig","D")
	self.out.branch("GSF_rochester","D")
	self.out.branch("GSF_dataset","D")
    	self.out.branch("GSF_HEM","D")

        ###

        def add_branches(ver):
            self.out.branch(ver+"RegionIso","D")
            self.out.branch(ver+"RegionOSSS","D")
            self.out.branch(ver+"Zvis_pt","D")
            self.out.branch(ver+"Zvis_eta","D")
            self.out.branch(ver+"Zvis_phi","D")
            self.out.branch(ver+"Zvis_mass","D")
            self.out.branch(ver+"Zvis_dR","D")
            self.out.branch(ver+"Index_tauobj","I")
            self.out.branch(ver+"Index_muon","I")
            self.out.branch(ver+"cut_Pzeta","D")
            self.out.branch(ver+"cut_MT","D")
            self.out.branch(ver+"HT","D")
            self.out.branch(ver+"NJets","D")
            self.out.branch(ver+"PassPzeta","I")
            self.out.branch(ver+"PassMT","I")
            self.out.branch(ver+"PassSel","I")
            if self.datamc == 'sigRes': 
              self.out.branch(ver+"DyDecayType","I")
        add_branches("AnaTau_")
        add_branches("AnaTp_")
        add_branches("AnaTpm_")

        # code snippet in case array branches are needed
        #self.out.branch("n"+"", "I")
        #self.out.branch(""+"_pt", "F", lenVar="n"+"")

#	self.out.branch("GMCTC","I")
#	self.out.branch("Z_massTPIso","D")
#	self.out.branch("Z_massTP3track","D")
#	self.out.branch("Z_massMuonMuon","D")

#	self.out.branch("GTau_iso_pt","D")
#	self.out.branch("GTau_noiso_pt","D")
#	self.out.branch("GTwoProng_iso_pt","D")
#	self.out.branch("GTwoProng_noiso_pt","D")
#	self.out.branch("GTwoProng_3track_pt","D")

#	self.out.branch("GTau_iso_phi","D")
#        self.out.branch("GTau_noiso_phi","D")
#        self.out.branch("GTwoProng_iso_phi","D")
#        self.out.branch("GTwoProng_noiso_phi","D")
#	self.out.branch("GTwoProng_3track_phi","D")

#	self.out.branch("GTau_iso_eta","D")
#        self.out.branch("GTau_noiso_eta","D")
#        self.out.branch("GTwoProng_iso_eta","D")
#        self.out.branch("GTwoProng_noiso_eta","D")
#	self.out.branch("GTwoProng_3track_eta","D")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        #self.out.fillBranch("GEventsRunOver",self.nevents)
        pass

    def analyze(self, event):
        jets = Collection(event, "Jet")
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        taus = Collection(event, "Tau")
        twoprongs = Collection(event, "TwoProng")
        #twoprongs = sorted(twoprongs, reverse = True, key = lambda obj : obj.pt)
        twoprongsmod = Collection(event, "TwoProngModified")
        #twoprongsmod = sorted(twoprongsmod, reverse = True, key = lambda obj : obj.pt)
        pv = Object(event, "PV")
        met = Object(event, "MET")
        hlt = Object(event, "HLT")
        pileup = Object(event, "Pileup")
        run = Object(event, "run")
        flags = Object(event, "Flag")
        #if DMC==1: genpart = Collection(event, "GenPart")
        year = 2018
        if self.datamc == 'sigRes': genparticles = Collection(event, "GenPart")

        if year == 2018:
            passEventFilter = flags.goodVertices and flags.HBHENoiseFilter and \
                    flags.HBHENoiseIsoFilter and flags.EcalDeadCellTriggerPrimitiveFilter \
                    and flags.BadPFMuonFilter and flags.globalSuperTightHalo2016Filter
            if self.datamc == "data":
                passEventFilter = passEventFilter and flags.eeBadScFilter
            else:
                pass
        '''
        if self.era == "2017" or self.era == "2018":
            ## Common filters
            ## Missing the latest ecalBadCalibReducedMINIAODFilter, not in MiniAOD
            ## But still using the old ecalBadCalibFilter from MiniAOD
            passEventFilter = flags.goodVertices and flags.HBHENoiseFilter and \
                    flags.HBHENoiseIsoFilter and flags.EcalDeadCellTriggerPrimitiveFilter \
                    and flags.BadPFMuonFilter and flags.ecalBadCalibFilter  
            ## Only data
            if self.isData:
                passEventFilter = passEventFilter and flags.globalSuperTightHalo2016Filter and flags.eeBadScFilter
            elif not self.isFastSim:
                passEventFilter = passEventFilter and flags.globalSuperTightHalo2016Filter
        '''

        def get_dy_decay_type(genparticles):
          decay = 0
          lep = 0
          leps = []
          for index, gp in enumerate(genparticles):
            pid = abs(gp.pdgId)
            if not (pid==11 or pid==13 or pid==15): continue
            if not gp.genPartIdxMother == -1: momid = genparticles[gp.genPartIdxMother].pdgId
            else: momid = -1
            sf = gp.statusFlags
            ispromt = bool(sf & 0b00000000000001)
            ishard =  bool(sf & 0b00000010000000)
            isfhard = bool(sf & 0b00000100000000)
            islast =  bool(sf & 0b10000000000000)
            if momid == 23 or ishard:
              decay += 1
              lep = pid
            if momid == 23 or ishard or isfhard:
              if islast: leps.append(index)
          def hadronic(l):
            if 111 in l or 211 in l or \
               130 in l or 310 in l or 311 in l or 321 in l or \
               323 in l or 221 in l:
              return True
            else: return False
          tp = -1
          if lep == 15:
            dp1 = []
            dp2 = []
            for gp in genparticles:
              pid = abs(gp.pdgId)
              sf = gp.statusFlags
              ispromt = bool(sf & 0b00000000000001)
              ishard =  bool(sf & 0b00000010000000)
              isfhard = bool(sf & 0b00000100000000)
              if not gp.genPartIdxMother == -1: momid = genparticles[gp.genPartIdxMother].pdgId
              else: momid = -1
              if gp.genPartIdxMother == leps[0]:
                dp1.append(pid)
              if gp.genPartIdxMother == leps[1]:
                dp2.append(pid)
            if 11 in dp1 and 11 in dp2: tp = 1
            elif 13 in dp1 and 13 in dp2: tp = 2
            elif 11 in dp1 and 13 in dp2: tp = 3
            elif 11 in dp2 and 13 in dp1: tp = 3
            elif 11 in dp1 and hadronic(dp2): tp = 4
            elif 11 in dp2 and hadronic(dp1): tp = 4
            elif 13 in dp1 and hadronic(dp2): tp = 5
            elif 13 in dp2 and hadronic(dp1): tp = 5
            elif hadronic(dp1) and hadronic(dp2): tp = 6
            else: tp = 0
          if lep == 15 and tp == 5: return 1 # tau_mu tau_had
          elif lep == 11: return 2 # el el
          elif lep == 13: return 3 # mu mu
          elif lep == 15 and tp == 1: return 4 # tau_e tau_e
          elif lep == 15 and tp == 2: return 5 # tau_mu tau_mu
          elif lep == 15 and tp == 3: return 6 # tau_e tau_mu
          elif lep == 15 and tp == 4: return 7 # tau_e tau_had
          elif lep == 15 and tp == 6: return 8 # tau_had tau_had
          else: return -1 # unknown

        def passTrigger():
            if year == 2018 or year == 2016: return hlt.IsoMu24
            elif year == 2017: return hlt.IsoMu27
            else: return False

        def passEVeto():
            for electron in electrons:
                if self.Electron_Veto(electron): return False
            return True

        def passBJetVeto():
            for jet in jets:
              if self.BJet_Veto(jet): return False
            return True

        def passMuVeto(index):
            if index == -1: return False
            for i, muon in enumerate(muons):
              if i == index: continue
              if self.Muon_Veto(muon): return False
            return True

        def passDiMuVeto():
          for i, muon1 in enumerate(muons):
            for j, muon2 in enumerate(muons):
              if j<=i: continue
              if self.DiMuon_Veto(muon1,muon2): return False
          return True

        def findMuon():
          index = -1
          for i, muon in enumerate(muons):
            if self.MuonID_sansIso(muon) and (index == -1 or muon.pt > muons[index].pt): index = i
          if index == -1: isoBit = -1
          elif self.MuonIsIso(muons[index]): isoBit = 1
          else: isoBit = 2 
          return index, isoBit
    
        def TauCandID(obj, ver):
          obj_vec = ROOT.Math.PtEtaPhiMVector(obj.pt, obj.eta, obj.phi, obj.mass)
          globaloverlap = False
          for muon in muons:
            if muon.pt < 5: continue
            if not muon.isGlobal: continue
            muon_vec = ROOT.Math.PtEtaPhiMVector(muon.pt, muon.eta, muon.phi, muon.mass)
            dr = ROOT.Math.VectorUtil.DeltaR(muon_vec,obj_vec)
            if dr < 0.1:
              globaloverlap = True
              break
          if ver == 'tau': additionalID = obj.idDecayModeOldDMs and (obj.leadTkPtOverTauPt*obj.pt) >5 and (obj.idAntiEleDeadECal &0x1) == 0x1 and (obj.idAntiMu &0x2) == 0x2
          elif ver == 'tp' or ver == 'tpm': additionalID = obj.isTight
          if obj.pt > 20 and abs(obj.eta)< 2.3 and globaloverlap == False and additionalID: return True

        def TauCandCharge(obj, ver):
            if ver == 'tau': return obj.charge
            elif ver == 'tp':
                if obj.CHpos_pt > obj.CHneg_pt: return 1
                else: return -1
            elif ver == 'tpm':
                if obj.nTracks == 3: return int(obj.CHextra_charge)
                else:
                    if obj.CHpos_pt > obj.CHneg_pt: return 1
                    else: return -1    
            else: return 0

        def findTauCand(mu_index, ver):
            if mu_index == -1: return -1, -1
            index = -1
            if ver == 'tau': collection = taus
            elif ver == 'tp': collection = twoprongs
            elif ver == 'tpm': collection = twoprongsmod
            for i, obj in enumerate(collection):
              if TauCandID(obj, ver) and (index == -1 or obj.pt > collection[index].pt): index = i
            if index == -1: return -1, -1
            obj = collection[index]
            mu_charge = muons[mu_index].charge
            taucand_charge = TauCandCharge(obj, ver)
            if mu_charge == 0 or taucand_charge == 0: osssBit = -1
            elif mu_charge * taucand_charge < 0: osssBit = 1
            else: osssBit = 2
            return index, osssBit

        def mydot(vec1, vec2):
            return vec1.X() * vec2.X() + vec1.Y() * vec2.Y()

        def mydphi(phi1, phi2):
            dphi = abs(phi2-phi1)
            if dphi > math.pi: dphi -= 2*math.pi
            if dphi < -math.pi: dphi += 2*math.pi
            return dphi

        def computeEventWideVars(mu_index, tau_index, ver):
            if mu_index == -1 or tau_index == -1: return 0, 0, 0
            mu_obj = muons[mu_index]
            if ver == 'tau': tau_obj = taus[tau_index]
            if ver == 'tp': tau_obj = twoprongs[tau_index]
            if ver == 'tpm': tau_obj = twoprongsmod[tau_index]
            tau_vec = ROOT.Math.PtEtaPhiMVector(tau_obj.pt, tau_obj.eta, tau_obj.phi, tau_obj.mass)
            mu_vec = ROOT.Math.PtEtaPhiMVector(mu_obj.pt, mu_obj.eta, mu_obj.phi, mu_obj.mass)

            dr = ROOT.Math.VectorUtil.DeltaR(mu_vec, tau_vec)

            mu_Tvec = ROOT.Math.Polar2DVector(mu_vec.Pt(), mu_vec.Phi())
            tau_Tvec = ROOT.Math.Polar2DVector(tau_vec.Pt(), tau_vec.Phi())
            met_Tvec = ROOT.Math.Polar2DVector(met.pt, met.phi)
            muhat = ROOT.Math.Polar2DVector(1.0, mu_vec.Phi())
            tauhat = ROOT.Math.Polar2DVector(1.0, tau_vec.Phi())
            zeta = muhat + tauhat
            zeta.SetR(1.0)
            pall = mydot((mu_Tvec + tau_Tvec + met_Tvec), zeta)
            pvis = mydot(mu_Tvec + tau_Tvec, zeta)
            pzeta = pall - 0.85*pvis

            dphi = mydphi(mu_vec.Phi(), met.phi)
            mt = pow(2*mu_vec.Pt()*met.pt*(1-math.cos(dphi)), 0.5)

            return dr, pzeta, mt

        def passEventWideCuts(deltar, pzeta, mt):
            return (deltar > CUT_ZPAIR_DR) and (pzeta > CUT_PZETA) and (mt < CUT_MT)

        def fill_branches(ver):
            if ver == 'tau': prefix = 'AnaTau_'
            if ver == 'tp': prefix = 'AnaTp_'
            if ver == 'tpm': prefix = 'AnaTpm_'
            Index_muonobj, RegionIso = findMuon()
            Index_tauobj, RegionOSSS = findTauCand(Index_muonobj, ver)
            haveMuonAndTauCand = Index_tauobj != -1
            deltar, pzeta, mt = computeEventWideVars(Index_muonobj, Index_tauobj, ver)
            passDeltaR = deltar > CUT_ZPAIR_DR
            passPzeta = pzeta > CUT_PZETA
            passMT = mt < CUT_MT
            passEVetoBit = passEVeto()
            passBJetVetoBit = passBJetVeto()
            passMuVetoBit = passMuVeto(Index_muonobj)
            passDiMuVetoBit = passDiMuVeto()
            passTriggerBit = passTrigger()

            passVetosBit = passEVetoBit and passBJetVetoBit and passMuVetoBit and passDiMuVetoBit
            passAllButVetosBit = passEventFilter and passTriggerBit and haveMuonAndTauCand and passDeltaR and passMT and passPzeta
            passAllButEventWideBit = passEventFilter and passTriggerBit and haveMuonAndTauCand and passDeltaR
            passBit = passAllButVetosBit and passVetosBit
        
            # metadata
            self.metadata[prefix+'ztt_total'] += 1

            if passMuVetoBit:    self.metadata[prefix+'ztt_passOnly_muVeto'] += 1
            if passEVetoBit:     self.metadata[prefix+'ztt_passOnly_elVeto'] += 1
            if passDiMuVetoBit:  self.metadata[prefix+'ztt_passOnly_diMuVeto'] += 1
            if passBJetVetoBit:  self.metadata[prefix+'ztt_passOnly_bVeto'] += 1

            if passAllButVetosBit and                  passBJetVetoBit and passMuVetoBit and passDiMuVetoBit: self.metadata[prefix+'ztt_nminusone_eveto'] += 1
            if passAllButVetosBit and passEVetoBit and                     passMuVetoBit and passDiMuVetoBit: self.metadata[prefix+'ztt_nminusone_bjetveto'] += 1
            if passAllButVetosBit and passEVetoBit and passBJetVetoBit and                   passDiMuVetoBit: self.metadata[prefix+'ztt_nminusone_muveto'] += 1
            if passAllButVetosBit and passEVetoBit and passBJetVetoBit and passMuVetoBit                    : self.metadata[prefix+'ztt_nminusone_dimuveto'] += 1

            if passVetosBit and passEventFilter and passTriggerBit and haveMuonAndTauCand and passDeltaR and            passPzeta: self.metadata[prefix+'ztt_nminusone_pzeta'] += 1
            if passVetosBit and passEventFilter and passTriggerBit and haveMuonAndTauCand and passDeltaR and passMT              : self.metadata[prefix+'ztt_nminusone_mt'] += 1

            if passTriggerBit: 
                self.metadata[prefix+'ztt_flow_trig'] += 1
                if passVetosBit:
                    self.metadata[prefix+'ztt_flow_allveto'] += 1
                    if passEventFilter and passTriggerBit and haveMuonAndTauCand and passDeltaR:
                        self.metadata[prefix+'ztt_flow_muAndTau'] += 1
                        if passMT and passPzeta:
                            self.metadata[prefix+'ztt_flow_fullsel'] += 1

            if passAllButEventWideBit:
              if ver == 'tau': taucand_obj = taus[Index_tauobj]
              if ver == 'tp': taucand_obj = twoprongs[Index_tauobj]
              if ver == 'tpm': taucand_obj = twoprongsmod[Index_tauobj]
              muon_obj = muons[Index_muonobj]
              taucand_vec = ROOT.Math.PtEtaPhiMVector(taucand_obj.pt, taucand_obj.eta, taucand_obj.phi, taucand_obj.mass)
              muon_vec = ROOT.Math.PtEtaPhiMVector(muon_obj.pt, muon_obj.eta, muon_obj.phi, muon_obj.mass)
              Zvis_vec = taucand_vec + muon_vec
              Zvis_pt = Zvis_vec.Pt()
              Zvis_eta = Zvis_vec.Eta()
              Zvis_phi = Zvis_vec.Phi()
              Zvis_mass = Zvis_vec.M()
              Zvis_dR = deltar
              cut_Pzeta = pzeta
              cut_MT = mt
              HT = muon_vec.Pt() + taucand_vec.Pt()
              nJets = 0
              for jet in jets:
                if jet.pt < 30: continue
                if abs(jet.eta) > 2.5: continue
                if not (jet.jetId & 0b10): continue
                jet_vec = ROOT.Math.PtEtaPhiMVector(jet.pt, jet.eta, jet.phi, jet.mass)
                if ROOT.Math.VectorUtil.DeltaR(muon_vec, jet_vec) < 0.3: continue
                if ROOT.Math.VectorUtil.DeltaR(taucand_vec, jet_vec) < 0.3: continue
                HT += jet.pt
                nJets += 1

              self.out.fillBranch(prefix+"RegionIso", RegionIso) # 1 = isolated, 2 = anti-isolated, -1 = fail
              self.out.fillBranch(prefix+"RegionOSSS", RegionOSSS) # 1 = OS, 2 = SS, -1 = fail
              self.out.fillBranch(prefix+"Zvis_pt", Zvis_pt)
              self.out.fillBranch(prefix+"Zvis_eta", Zvis_eta)
              self.out.fillBranch(prefix+"Zvis_phi", Zvis_phi)
              self.out.fillBranch(prefix+"Zvis_mass", Zvis_mass)
              self.out.fillBranch(prefix+"Zvis_dR", Zvis_dR)
              self.out.fillBranch(prefix+"Index_tauobj", Index_tauobj)
              self.out.fillBranch(prefix+"Index_muon", Index_muonobj)
              self.out.fillBranch(prefix+"cut_Pzeta", cut_Pzeta)
              self.out.fillBranch(prefix+"cut_MT", cut_MT)
              self.out.fillBranch(prefix+"HT", HT)
              self.out.fillBranch(prefix+"NJets", nJets)
              self.out.fillBranch(prefix+"PassPzeta", passPzeta)
              self.out.fillBranch(prefix+"PassMT", passMT)
              self.out.fillBranch(prefix+"PassSel", passBit)
              if self.datamc == 'sigRes':
                self.out.fillBranch(prefix+"DyDecayType", get_dy_decay_type(genparticles))
              return True
            else:
              self.out.fillBranch(prefix+"RegionIso", -2) # 1 = isolated, 2 = anti-isolated, -1 = fail
              self.out.fillBranch(prefix+"RegionOSSS", -2) # 1 = OS, 2 = SS, -1 = fail
              self.out.fillBranch(prefix+"Zvis_pt", -1)
              self.out.fillBranch(prefix+"Zvis_eta", 10)
              self.out.fillBranch(prefix+"Zvis_phi", 10)
              self.out.fillBranch(prefix+"Zvis_mass", -1)
              self.out.fillBranch(prefix+"Zvis_dR", -1)
              self.out.fillBranch(prefix+"Index_tauobj", -2)
              self.out.fillBranch(prefix+"Index_muon", -2)
              self.out.fillBranch(prefix+"cut_Pzeta", -1000)
              self.out.fillBranch(prefix+"cut_MT", -1)
              self.out.fillBranch(prefix+"HT", -1)
              self.out.fillBranch(prefix+"NJets", -1)
              self.out.fillBranch(prefix+"PassPzeta", -1)
              self.out.fillBranch(prefix+"PassMT", -1)
              if self.datamc == 'sigRes':
                self.out.fillBranch(prefix+"DyDecayType", -2)
              return False
        
        bit1 = fill_branches('tau')
        bit2 = fill_branches('tp')
        bit3 = fill_branches('tpm')
        return bit1 or bit2 or bit3

########

	###Data or MC Tag and other Tags###
	DMC= self.DMC ##0 for data, 1 for MC
	yeartag = self.yeartag ##2018,2017, or 2016
	antiiso= self.antiiso #antiiso if set to 1, iso if set to 0 (iso is normal run)
	self.nevents+=1
	mcnum = self.mcnum ##will correspond to a JSON file to pull the MC to data Scalefactor

	sf_pileup_c=1
        sf_pileup_p=1
        sf_pileup_m=1
        sf_id_iso_trig=1

	EVeto = False
	BJetVeto = False
	MVeto = False
	DiMVeto = False
	PassMuonCuts = False
#	PassMuonCuts2 = False
	PassTauCuts = False
	PassTwoProngCuts = False
	HaveGoodMuon = False
	HaveGoodTau = False
	HaveGoodTP = False
	HaveTwoGoodMuon = False
	IsoTP = False
	pzetaTau = False
	pzetaTP = False
        GoodMuon_pt = 0.0
	GoodMuon_eta = 0.0
	GoodMuon_phi = 0.0
	GoodMuon_mass = 0.0
	GoodMuon_charge = 0.0
        GoodMuon_number = 0
	GoodMuon_pfisoid=0
	GoodMuon_genpt= -50
	GoodMuon = None

	GMET= 0
	GoodTau_pt = 0.0
	GoodTau_eta = 0.0
	GoodTau_phi = 0.0
	GoodTau_mass =  0.0
	GoodTau_charge = 0.0
	GoodTau_number = 0
	GoodTau = None
    	hem_sf=1
#	GoodMuon_pt2 = 0.0
#        GoodMuon_eta2 = 0.0
#        GoodMuon_phi2 = 0.0
#        GoodMuon_mass2 = 0.0
#        GoodMuon_charge2 = 0.0
#        GoodMuon_number2 = 0
#        GoodMuon_genpt2= -50
#        GoodMuon2 = None

	gnpv=0

	GoodTwoProng_pt = 0.0
        GoodTwoProng_eta = 0.0
        GoodTwoProng_phi = 0.0
        GoodTwoProng_number = 0
        GoodTwoProng = None
	GoodTwoProng_ntracks= 0
	GoodMuon_ntrackerlayer=0
	InvarMassZVTau = 0
	InvarMassZVTP = 0

	GnJets=0
	jetskips=0
	sf_central = 1
	sf_plus = 1
	sf_minus = 1
	sf_roc = 1


	"""process event, return True (go to next module) or False (fail, go to next event)"""

	jets = Collection(event, "Jet")
	electrons = Collection(event, "Electron")
	muons = Collection(event, "Muon")
	taus = Collection(event, "Tau")
	twoprongs = Collection(event, "TwoProngModified")
	twoprongs = sorted(twoprongs, reverse=True, key=lambda obj : obj.pt)
	if DMC==1: genpart = Collection(event, "GenPart")
	pv= Object(event, "PV")
	met = Object(event, "MET")
	hlt = Object(event,"HLT")
	pileup = Object(event,"Pileup")
	run = Object(event,"run")

	if DMC==1: self.out.fillBranch("GPileup_nTrueInt", pileup.nTrueInt)

	else: self.out.fillBranch("GPileup_nTrueInt", -50)
#    	if antiiso==0:
	if 1 ==1:
        	if yeartag == 2018 or yeartag == 2016:
		        if hlt.IsoMu24 == False:
			        return False
		if yeartag == 2017:
			if hlt.IsoMu27 == False:
				return False

	gnpv=pv.npvs

	#electron veto
	#electron_branch_name = 'Electron'
	for i, electron in enumerate(electrons):
		EVeto = self.Electron_Veto(electron)
		if EVeto: break
	if EVeto:
		return False

	#jet veto
	for i, jet in enumerate (jets):
        	BJetVeto= self.BJet_Veto(jet)
		if BJetVeto: break
	if BJetVeto:
		return False


	#DiMuon Veto
	for i, muon1 in enumerate (muons):
		for j, muon2 in enumerate (muons):
			if j>i:
				DiMVeto = self.DiMuon_Veto(muon1,muon2)
			if DiMVeto: break
		if DiMVeto: break
	if DiMVeto:
		return False


    	#Muon Cuts
	for i, muon in enumerate(muons):
        	if antiiso ==1: PassMuonCuts = self.AntiIso_Cuts(muon)
		else: PassMuonCuts = self.Muon_Cuts(muon)
		if PassMuonCuts:
			if muon.pt > GoodMuon_pt:
				GoodMuon_pt = muon.pt
				GoodMuon_eta = muon.eta
				GoodMuon_phi = muon.phi
				GoodMuon_mass = muon.mass
				GoodMuon_charge = math.copysign(1,muon.charge)
				GoodMuon_number = i
				GoodMuon = muon
				GoodMuon_pfisoid = muon.pfIsoId
				if DMC==1:
					if muon.genPartFlav !=0:
						GoodMuon_genpt= genpart[muon.genPartIdx].pt
					else:
						GoodMuon_genpt = -1
				GoodMuon_ntrackerlayer= muon.nTrackerLayers
				HaveGoodMuon = True


	if HaveGoodMuon==False: return False


	#Double Muon Zmass
#	for i, muon in enumerate(muons):
#                PassMuonCuts2 = self.Muon_Cuts(muon)
#                if PassMuonCuts2:
#			gmV2 =  ROOT.Math.PtEtaPhiMVector(muon.pt, muon.eta, muon.phi, muon.mass)
#        		gmV = ROOT.Math.PtEtaPhiMVector(GoodMuon.pt, GoodMuon.eta, GoodMuon.phi, GoodMuon.mass)
#        		delrmm= ROOT.Math.VectorUtil.DeltaR(gmV,gmV2)
#			if muon != GoodMuon and muon.pt != GoodMuon_pt and delrmm >=0.1 and self.MCTC_Cut(muon,GoodMuon):
#                        	if muon.pt > GoodMuon_pt2:
#                                	GoodMuon_pt2 = muon.pt
#                                	GoodMuon_eta2 = muon.eta
#                                	GoodMuon_phi2 = muon.phi
#                                	GoodMuon_mass2 = muon.mass
#                                	GoodMuon_charge2 = math.copysign(1,muon.charge)
#                                	GoodMuon_number2 = i
#                                	GoodMuon2 = muon
#                                	if DMC==1:
#                                        	if muon.genPartFlav !=0:
#                                                	GoodMuon_genpt2= genpart[muon.genPartIdx].pt
#                                        	else:
#                                            	GoodMuon_genpt2 = -1
#                                	GoodMuon_ntrackerlayer2= muon.nTrackerLayers
#					if GoodMuon_pt2 > 0: HaveTwoGoodMuon = True

	#Muon Veto if there is a muon that passes cuts
	for i, muon in enumerate(muons):
		if i!= GoodMuon_number and GoodMuon_pt >0:
			MVeto = self.Muon_Veto(muon)
		if MVeto: break
	if MVeto:
		return False



	#Tau Cuts
	for i, tau in enumerate(taus):
		if GoodMuon != None:
			PassTauCuts= self.Tau_Cuts(tau, GoodMuon, muons)
			if PassTauCuts:
				if tau.pt >GoodTau_pt:
					GoodTau_eta = tau.eta
					GoodTau_pt = tau.pt
					GoodTau_phi= tau.phi
					GoodTau_mass = tau.mass
					GoodTau_charge = math.copysign(1,tau.charge)
					GoodTau_number = i
					GoodTau = tau
					HaveGoodTau = True

	#Transmass Cut

	transmass = self.TransMass_Cut(GoodMuon_pt, met, GoodMuon_phi)
	if transmass >= 30: return False


	#TwoProng Cuts
	for i, twoprong in enumerate(twoprongs):
		if GoodMuon != None:
			PassTwoProngCuts = self.TwoProng_Cuts(twoprong, GoodMuon, muons,run)
			if PassTwoProngCuts:
				if twoprong.pt > GoodTwoProng_pt:
					GoodTwoProng_eta = twoprong.eta
                                        GoodTwoProng_pt = twoprong.pt
                                        GoodTwoProng_phi= twoprong.phi
                                        GoodTwoProng_number = i
                                        GoodTwoProng = twoprong
					try: GoodTwoProng_ntracks = twoprong.nTracks
        				except RuntimeError: GoodTwoProng_ntracks = 2
					if -3 < twoprong.eta < -1.4 and -1.57 < twoprong.phi < -0.87 and DMC==1: hem_sf = 0.35
                    			else: hem_sf = 1
					HaveGoodTP = True
    	PassJetCuts= False
    	for i, jet in enumerate(jets):
        	PassJetCuts= self.Jet_Cuts(jet)
        	if PassJetCuts:
        	    break
	if HaveGoodTau == False and HaveGoodTP == False and PassJetCuts ==False:
		return False
#	if HaveGoodTP==False: return False

	#rochester corrections
	nrandom=np.random.random()

    	if DMC==1:
            if GoodMuon_genpt != -1 : sf_roc= self.kSpreadMC(self.roccin, int(GoodMuon_charge), GoodMuon_pt, GoodMuon_eta ,GoodMuon_phi,GoodMuon_genpt, 0,0 )
            else: sf_roc= self.kSmearMC(self.roccin, int(GoodMuon_charge), GoodMuon_pt, GoodMuon_eta ,GoodMuon_phi ,int(GoodMuon_ntrackerlayer), nrandom, 0,0 )
    	else: sf_roc= self.kScaleDT(self.roccin, int(GoodMuon_charge), GoodMuon_pt, GoodMuon_eta, GoodMuon_phi, 0, 0)



	if HaveGoodTau:
		pzetaTau = self.PZeta_Cut(GoodMuon, GoodTau, met,sf_roc)
	if HaveGoodTP:
		pzetaTP = self.PZeta_Cut(GoodMuon, GoodTwoProng, met,sf_roc)
	#if pzetaTau != True  and pzetaTP != True:
	#	return False

	#Invar Masses
	#if HaveGoodTau and pzetaTau:
    	if HaveGoodTau:
		ZInvarMassVTau = self.ZInvarMass(GoodMuon, GoodTau,sf_roc).M()
        	ZInvarPtVTau = self.ZInvarMass(GoodMuon, GoodTau,sf_roc).Pt()
        	ZInvarPhiVTau = self.ZInvarMass(GoodMuon, GoodTau,sf_roc).Phi()
        	ZInvarEtaVTau = self.ZInvarMass(GoodMuon, GoodTau,sf_roc).Eta()
	# HaveGoodTP and pzetaTP:
    	if HaveGoodTP:
		ZInvarMassVTP =  self.ZInvarMass(GoodMuon, GoodTwoProng,sf_roc).M()
        	ZInvarPtVTP =  self.ZInvarMass(GoodMuon, GoodTwoProng,sf_roc).Pt()
        	ZInvarPhiVTP =  self.ZInvarMass(GoodMuon, GoodTwoProng,sf_roc).Phi()
        	ZInvarEtaVTP =  self.ZInvarMass(GoodMuon, GoodTwoProng,sf_roc).Eta()
    	else: hem_sf = 1
#	if HaveTwoGoodMuon: ZInvarMassVMuMu = self.ZInvarMass(GoodMuon, GoodMuon2,sf_roc)

	if HaveGoodTau and HaveGoodTP:
		goodtpV =  ROOT.Math.PtEtaPhiMVector(GoodTwoProng.pt, GoodTwoProng.eta, GoodTwoProng.phi, GoodTwoProng.mass)
        	goodtauV = ROOT.Math.PtEtaPhiMVector(GoodTau.pt, GoodTau.eta, GoodTau.phi, GoodTau.mass)
        	deltargood= ROOT.Math.VectorUtil.DeltaR(goodtauV,goodtpV)
		if deltargood < 0.1:
			IsoTP = True

		else:
			IsoTP = False

	GMET=met.pt

        #Jet Correction
        ht=0
        for i, jet in enumerate(jets):
                if jet.pt >30 and abs(jet.eta) < 2.5 and (jet.jetId & 0b10):
                        PassJetCleaner= True
                        PassJetCleaner = self.Jet_Cleaner(jet,GoodMuon)
                        if HaveGoodTau and PassJetCleaner == True: PassJetCleaner = self.Jet_Cleaner(jet,GoodTau)
			if HaveGoodTP and PassJetCleaner == True: PassJetCleaner = self.Jet_Cleaner(jet,GoodTwoProng)
                        if PassJetCleaner==False:
				jetskips+=1
				self.out.fillBranch("GJet_pt",-100)
                                self.out.fillBranch("GJet_eta",-100)
                                self.out.fillBranch("GJet_phi",-100)
                                self.out.fillBranch("GJet_mass",-100)
                        else:
				self.out.fillBranch("GJet_pt",jet.pt)
				self.out.fillBranch("GJet_eta",jet.eta)
                                self.out.fillBranch("GJet_phi",jet.phi)
                                self.out.fillBranch("GJet_mass",jet.mass)
                                ht+=jet.pt
        	else: jetskips+=1
        GnJets = event.nJet - jetskips
	if GnJets <0 : GnJets=0
    	self.out.fillBranch("Ght",ht)
	###Scale Factor Section###
	if DMC==1:
		sf_pileup_c= self.centralfunc.GetBinContent(int(pileup.nTrueInt)+1)
		sf_pileup_p= self.p4_6func.GetBinContent(int(pileup.nTrueInt)+1)
		sf_pileup_m= self.m4_6func.GetBinContent(int(pileup.nTrueInt)+1)
	#	sf_id_iso_trig= self.MediumID_MediumIso_ScaleFactor(GoodMuon_eta, GoodMuon_pt,yeartag)
		sf_id_iso_trig= self.TightID_TightIso_ScaleFactor(GoodMuon_eta, GoodMuon_pt,yeartag)
		sf_central=sf_id_iso_trig*sf_pileup_c*self.datasetscalefactor
		sf_plus=sf_id_iso_trig*sf_pileup_p*self.datasetscalefactor
		sf_minus=sf_id_iso_trig*sf_pileup_m*self.datasetscalefactor





	###Branch Filling Section###

	self.out.fillBranch("GSF_central",sf_central)
	self.out.fillBranch("GSF_plus",sf_plus)
	self.out.fillBranch("GSF_minus",sf_minus)
	self.out.fillBranch("GSF_pileup_central",sf_pileup_c)
        self.out.fillBranch("GSF_pileup_plus",sf_pileup_p)
        self.out.fillBranch("GSF_pileup_minus",sf_pileup_m)
	self.out.fillBranch("GSF_id_iso_trig",sf_id_iso_trig)
	self.out.fillBranch("GSF_rochester",sf_roc)
        self.out.fillBranch("GSF_HEM",hem_sf)



#	if MCTC == True: self.out.fillBranch("GMCTC",1)
#	else: self.out.fillBranch("GMCTC",0)
	self.out.fillBranch("nMuons",event.nMuon)
	self.out.fillBranch("nTaus",event.nTau)
	self.out.fillBranch("nJets",event.nJet)

	self.out.fillBranch("GMET", GMET)
	self.out.fillBranch("GNPV",gnpv)
	self.out.fillBranch("GSF_dataset",self.datasetscalefactor)
        self.out.fillBranch("GnJets_corrected",GnJets)

	self.out.fillBranch("GMuon_pt", GoodMuon_pt)
	self.out.fillBranch("GMuon_eta", GoodMuon_eta)
	self.out.fillBranch("GMuon_phi", GoodMuon_phi)
	self.out.fillBranch("GMuon_mass",GoodMuon_mass)
	self.out.fillBranch("GMuon_charge", GoodMuon_charge)
	self.out.fillBranch("GMuon_genpt", GoodMuon_genpt)
	self.out.fillBranch("GMuon_ntrackerlayers",GoodMuon_ntrackerlayer)
	self.out.fillBranch("GMuon_pfisoid",GoodMuon_pfisoid)
	self.out.fillBranch("GTransMass",transmass)

    	if HaveGoodTau:
            self.out.fillBranch("GTau_pt", GoodTau_pt)
            self.out.fillBranch("GTau_eta", GoodTau_eta)
            self.out.fillBranch("GTau_phi", GoodTau_phi)
            self.out.fillBranch("GTau_mass",GoodTau_mass)
            self.out.fillBranch("Z_massTau", ZInvarMassVTau)
            self.out.fillBranch("Z_ptTau", ZInvarPtVTau)
            self.out.fillBranch("Z_phiTau", ZInvarPhiVTau)
            self.out.fillBranch("Z_etaTau", ZInvarEtaVTau)
            self.out.fillBranch("GTau_charge",GoodTau_charge)
	    self.out.fillBranch("GPZeta_Tau",pzetaTau)
    	else:
		self.out.fillBranch("GTau_pt", -100)
        	self.out.fillBranch("GTau_eta", -100)
        	self.out.fillBranch("GTau_phi", -100)
        	self.out.fillBranch("GTau_mass", -100)
        	self.out.fillBranch("Z_massTau", -100)
            	self.out.fillBranch("Z_phiTau", -100)
            	self.out.fillBranch("Z_ptTau", -100)
            	self.out.fillBranch("Z_etaTau", -100)
        	self.out.fillBranch("GTau_charge",-100)
                self.out.fillBranch("GPZeta_Tau",-100)

	if HaveGoodTP:
        	self.out.fillBranch("GTwoProng_pt", GoodTwoProng_pt)
        	self.out.fillBranch("GTwoProng_eta", GoodTwoProng_eta)
        	self.out.fillBranch("GTwoProng_phi", GoodTwoProng_phi)
        	self.out.fillBranch("GTwoProng_mass", GoodTwoProng.mass)
        	self.out.fillBranch("GTwoProng_massl", GoodTwoProng.massl)
        	self.out.fillBranch("GTwoProng_massPi0", GoodTwoProng.massPi0)
        	self.out.fillBranch("GTwoProng_massEta",GoodTwoProng.massEta)
        	self.out.fillBranch("Z_massTP", ZInvarMassVTP)
            	self.out.fillBranch("Z_ptTP", ZInvarPtVTP)
            	self.out.fillBranch("Z_phiTP", ZInvarPhiVTP)
            	self.out.fillBranch("Z_etaTP", ZInvarEtaVTP)
            	self.out.fillBranch("GChpospt",GoodTwoProng.CHpos_pt)
            	self.out.fillBranch("GChnegpt",GoodTwoProng.CHneg_pt)
		self.out.fillBranch("GChextra_charge",GoodTwoProng.CHextra_charge)
        	self.out.fillBranch("GTPnTracks",GoodTwoProng_ntracks)
                self.out.fillBranch("GPZeta_TP",pzetaTP)

	else:
	 	self.out.fillBranch("GTwoProng_pt", -100)
        	self.out.fillBranch("GTwoProng_eta", -100)
        	self.out.fillBranch("GTwoProng_phi", -100)
        	self.out.fillBranch("GTwoProng_mass", -100)
        	self.out.fillBranch("GTwoProng_massl", -100)
        	self.out.fillBranch("GTwoProng_massPi0", -100)
        	self.out.fillBranch("GTwoProng_massEta", -100)
            	self.out.fillBranch("GChpospt",-100)
            	self.out.fillBranch("GChnegpt",-100)
		self.out.fillBranch("GChextra_charge",-100)
        	self.out.fillBranch("Z_massTP", -100)
            	self.out.fillBranch("Z_ptTP", -100)
            	self.out.fillBranch("Z_phiTP", -100)
            	self.out.fillBranch("Z_etaTP", -100)
        	self.out.fillBranch("GTPnTracks",-100)
		self.out.fillBranch("GPZeta_TP",-100)

	if HaveGoodTP and HaveGoodTau:
		self.out.fillBranch("TPTauDeltaR",deltargood)
	else:
		self.out.fillBranch("TPTauDeltaR",-100)


	if HaveGoodTP:
		self.out.fillBranch("GMuon_pt_tpcheck",GoodMuon.pt)
	else:
		self.out.fillBranch("GMuon_pt_tpcheck",-100)












#	if HaveGoodTP and HaveGoodTau and IsoTP and pzetaTau and pzetaTP and HaveGoodMuon:
#		self.out.fillBranch("Z_massTPIso", ZInvarMassVTP)
#		self.out.fillBranch("GTwoProng_iso_pt", GoodTwoProng.pt)
#		self.out.fillBranch("GTwoProng_iso_phi",GoodTwoProng.phi)
#		self.out.fillBranch("GTwoProng_iso_eta",GoodTwoProng.eta)
#		self.out.fillBranch("GTau_iso_pt",GoodTau.pt)
#		self.out.fillBranch("GTau_iso_phi", GoodTau.phi)
#		self.out.fillBranch("GTau_iso_eta", GoodTau.eta)
#		self.out.fillBranch("TPTauDeltaR",deltargood)
#	else:
#		self.out.fillBranch("Z_massTPIso", -100)

#		if HaveGoodTP and HaveGoodTau and IsoTP == False and pzetaTP and pzetaTau and HaveGoodMuon:
#			self.out.fillBranch("TPTauDeltaR",deltargood)
#			self.out.fillBranch("GTwoProng_iso_pt", -100)
#                	self.out.fillBranch("GTwoProng_iso_phi",-100)
#                	self.out.fillBranch("GTwoProng_iso_eta",-100)
#			self.out.fillBranch("GTau_iso_pt",-100)
#                	self.out.fillBranch("GTau_iso_phi", -100)
#	               	self.out.fillBranch("GTau_iso_eta",-100 )
#			self.out.fillBranch("GTwoProng_noiso_pt", GoodTwoProng.pt)
#			self.out.fillBranch("GTwoProng_noiso_phi",GoodTwoProng.phi)
#                	self.out.fillBranch("GTwoProng_noiso_eta",GoodTwoProng.eta)
#                	self.out.fillBranch("GTau_noiso_pt",GoodTau.pt)
#                	self.out.fillBranch("GTau_noiso_phi", GoodTau.phi)
#                	self.out.fillBranch("GTau_noiso_eta", GoodTau.eta)

#		else:
#			self.out.fillBranch("GTwoProng_iso_pt", -100)
#                        self.out.fillBranch("GTwoProng_iso_phi",-100)
#                        self.out.fillBranch("GTwoProng_iso_eta",-100)
#                        self.out.fillBranch("GTau_iso_pt",-100)
#                        self.out.fillBranch("GTau_iso_phi", -100)
#                        self.out.fillBranch("GTau_iso_eta",-100 )
#			self.out.fillBranch("TPTauDeltaR",-500)
#                        self.out.fillBranch("GTwoProng_noiso_pt", -100)
#                        self.out.fillBranch("GTwoProng_noiso_phi",-100)
#                        self.out.fillBranch("GTwoProng_noiso_eta",-100)
#                        self.out.fillBranch("GTau_noiso_pt",-100)
#                        self.out.fillBranch("GTau_noiso_phi",-100)
#                        self.out.fillBranch("GTau_noiso_eta", -100)

#	if HaveGoodTP and GoodTwoProng_ntracks==3 and pzetaTP and HaveGoodMuon:
#			self.out.fillBranch("GTwoProng_3track_eta",GoodTwoProng.eta)
#			self.out.fillBranch("GTwoProng_3track_pt",GoodTwoProng.pt)
#			self.out.fillBranch("GTwoProng_3track_phi",GoodTwoProng.phi)
#			self.out.fillBranch("Z_massTP3track",ZInvarMassVTP)
#			self.out.fillBranch("GTPnTracks",GoodTwoProng_ntracks)

#	else:
#			self.out.fillBranch("GTwoProng_3track_eta",-100)
#                        self.out.fillBranch("GTwoProng_3track_pt",-100)
#                        self.out.fillBranch("GTwoProng_3track_phi",-100)
#                        self.out.fillBranch("Z_massTP3track",-100)
#			self.out.fillBranch("GTPnTracks",-100)

#	if HaveTwoGoodMuon and HaveGoodMuon:
#		self.out.fillBranch("Z_massMuonMuon",ZInvarMassVMuMu)
#	        self.out.fillBranch("GMuon_pt2", GoodMuon_pt2)
#        	self.out.fillBranch("GMuon_eta2", GoodMuon_eta2)
#        	self.out.fillBranch("GMuon_phi2", GoodMuon_phi2)
#        	self.out.fillBranch("GMuon_mass2",GoodMuon_mass2)
#        	self.out.fillBranch("GMuon_charge2", GoodMuon_charge2)
#        	self.out.fillBranch("GMuon_genpt2", GoodMuon_genpt2)
#        	self.out.fillBranch("GMuon_ntrackerlayers2",GoodMuon_ntrackerlayer2)
#	else:
#		self.out.fillBranch("Z_massMuonMuon",-100)
#                self.out.fillBranch("GMuon_pt2", -100)
#                self.out.fillBranch("GMuon_eta2", -100)
#                self.out.fillBranch("GMuon_phi2", -100)
#                self.out.fillBranch("GMuon_mass2",-100)
#                self.out.fillBranch("GMuon_charge2", -100)
#                self.out.fillBranch("GMuon_genpt2", -100)
#                self.out.fillBranch("GMuon_ntrackerlayers2",-100)


	return True

















    def Electron_Veto(self, electron):
	if electron.pt >10 and abs(electron.eta) < 2.5 and abs(electron.dxy) <0.045 and abs(electron.dz) < 0.2 and electron.mvaFall17V2Iso_WP90 and electron.lostHits <= 1 and electron.pfRelIso03_all <0.3: return True
	else: return False


    def BJet_Veto(self, jet):
	if jet.btagCSVV2 >0.89 and jet.pt > 20 and abs(jet.eta) <2.4 : return True
	else: return False
	#if jet.pt > 20 and abs(jet.eta) <2.4 : return True

#    def Mom_Organizer(self, GenPart):
#	pdgid = 0
#	genpartz = 0
#	mom = GenPart.genPartIdxMother

#	if GenPart.pdgId == 23: genpartz = 1 #if genpart is Z
#	if GenPart.pdgId == 11:
#		if mom == 23: pdgid = 11 # electron event
#	if GenPart.pdgId == 13:
#		if mom == 23: pdgid = 13 # mu event
#        if GenPart.pdgId == 15:
#		if mom == 23: pdgid = 15 #tau event
#	return pdgid, genpartz

    def DiMuon_Veto(self, muon1, muon2):
	muon1V= ROOT.Math.PtEtaPhiMVector(muon1.pt, muon1.eta, muon1.phi, muon1.mass)
        muon2V= ROOT.Math.PtEtaPhiMVector(muon2.pt, muon2.eta, muon2.phi, muon2.mass)
        deltarm= ROOT.Math.VectorUtil.DeltaR(muon1V,muon2V)
        if deltarm>0.15 and muon2.pt>15 and  muon1.pt> 15 and abs(muon2.eta)<2.4 and abs(muon1.eta)<2.4 and np.signbit(muon2.charge)!= np.signbit(muon1.charge) and muon2.isGlobal and muon1.isGlobal and muon2.isTracker and muon1.isTracker and muon2.isPFcand and muon1.isPFcand and abs(muon2.dz)<0.2 and abs(muon1.dz)<0.2 and abs(muon2.dxy)<0.045 and abs(muon1.dxy) <0.045 and muon2.pfIsoId >=1 and muon1.pfIsoId >=1: return True
	else: return False

    def MuonID_sansIso(self, muon):
	    #if muon.pt > 28 and abs(muon.eta) < 2.1 and muon.mediumId and abs(muon.dz) < 0.2 and abs(muon.dxy) < 0.045: return True
	    if muon.pt > 28 and abs(muon.eta) < 2.1 and muon.tightId and abs(muon.dz) < 0.2 and abs(muon.dxy) < 0.045: return True
            else: return False
     
    def MuonIsIso(self, muon):
      return muon.pfIsoId >=4 


    def Muon_Cuts(self, muon):
	if muon.pt >28 and abs(muon.eta) <2.1 and muon.tightId and muon.pfIsoId >=4 and muon.dz < 0.2 and abs(muon.dxy) <0.045: return True
	##changing to tight temp
	else: return False

    def AntiIso_Cuts(self, muon):
        if muon.pt >28 and abs(muon.eta) <2.1 and muon.tightId and muon.pfIsoId <4 and muon.dz < 0.2 and abs(muon.dxy) <0.045: return True
    #    if muon.pfIsoId <3 and muon.dz < 0.2 : return True
        else: return False



    def Tau_Cuts(self, tau, goodmuon, muons):
	globaloverlap = False
	goodmuonV = ROOT.Math.PtEtaPhiMVector(goodmuon.pt, goodmuon.eta, goodmuon.phi, goodmuon.mass)
	tauV =  ROOT.Math.PtEtaPhiMVector(tau.pt, tau.eta, tau.phi, tau.mass)
	deltarmt= ROOT.Math.VectorUtil.DeltaR(goodmuonV,tauV)
	if deltarmt < 0.5: return False
	for i, muon in enumerate(muons):
		muonVector= ROOT.Math.PtEtaPhiMVector(muon.pt, muon.eta, muon.phi, muon.mass);
            	deltart= ROOT.Math.VectorUtil.DeltaR(muonVector,tauV);
		if muon.isGlobal and deltart <0.1 and muon.pt >5 : globaloverlap= True
	if tau.pt > 20 and abs(tau.eta)< 2.3 and tau.idDecayModeOldDMs and (tau.leadTkPtOverTauPt*tau.pt) >5 and (tau.idAntiEleDeadECal &0x1) == 0x1 and (tau.idAntiMu &0x2) == 0x2 and globaloverlap == False: return True
	else: return False

#    def Tau_Cuts(self, tau, goodmuon, muons):
#	return True

    def Muon_Veto(self, muon):
	if muon.pt > 10 and abs(muon.eta) < 2.4 and abs(muon.dxy) < 0.045 and abs(muon.dz) < 0.2 and muon.mediumId and muon.pfIsoId >= 1: return True
	else: return False



    def TwoProng_Cuts(self, twoprong, goodmuon, muons,run):
	globaloverlap = False
	tpV =  ROOT.Math.PtEtaPhiMVector(twoprong.pt, twoprong.eta, twoprong.phi, twoprong.mass)
	goodmuonV = ROOT.Math.PtEtaPhiMVector(goodmuon.pt, goodmuon.eta, goodmuon.phi, goodmuon.mass)
	deltargood= ROOT.Math.VectorUtil.DeltaR(goodmuonV,tpV)
	if deltargood <0.5: return False
#	if deltargood <0.1: return False
	for i, muon in enumerate(muons):
		muonV = ROOT.Math.PtEtaPhiMVector(muon.pt, muon.eta, muon.phi, muon.mass)
		deltart= ROOT.Math.VectorUtil.DeltaR(muonV,tpV)
		if muon.isGlobal and deltart <0.1 and muon.pt > 5: globaloverlap = True
	if globaloverlap: return False
	try: tight = twoprong.isTight
	except RuntimeError: tight = True

	if -3 < twoprong.eta < -1.4 and -1.57 < twoprong.phi < -0.87:
        	if self.DMC==0 and run >=319077: return False
	if twoprong.pt >20 and abs(twoprong.eta)< 2.3 and tight: return True

	return False

    def Jet_Cuts(self, jet):
	if jet.pt >20 and abs(jet.eta)< 2.3: return True

	return False

    def Jet_Cleaner(self, jet,part):
	partV = ROOT.Math.PtEtaPhiMVector(part.pt, part.eta, part.phi, part.mass)
	jetV =  ROOT.Math.PtEtaPhiMVector(jet.pt, jet.eta, jet.phi, jet.mass)
	deltarj= ROOT.Math.VectorUtil.DeltaR(partV,jetV)
	if deltarj<0.3: return False
        else: return True



    def TransMass_Cut(self, MuPtV, met, MuPhiV):
	transversemass =  pow(2*MuPtV*met.pt*(1-math.cos(MuPhiV-met.phi)),0.5)
	return transversemass


    def PZeta_Cut(self, muon, tau, met,mcorrections):

	muvec= ROOT.TVector3()
        muvec.SetPtEtaPhi(muon.pt*mcorrections,muon.eta,0)
        tauvec= ROOT.TVector3()
        tauvec.SetPtEtaPhi(tau.pt,tau.eta,0)
        etvec= ROOT.TVector3()
        etvec.SetPtEtaPhi(met.pt,0,met.phi)
        zetaU=ROOT.TVector3()
	zetaU= muvec*(1/muvec.Mag()) +tauvec*(1/tauvec.Mag())
	pall=(muvec + tauvec + etvec).Dot(zetaU)
        pvis=(muvec+ tauvec).Dot(zetaU)
        pzeta=pall-0.85*(pvis)
	if pzeta > -25: return True
	else: return False


    def ZInvarMass(self, muon, tau, mcorrections):
	muonV = ROOT.Math.PtEtaPhiMVector(muon.pt*mcorrections, muon.eta, muon.phi, muon.mass)
        tauV =  ROOT.Math.PtEtaPhiMVector(tau.pt, tau.eta, tau.phi, tau.mass)
	ZV = muonV + tauV
	return ZV

    def MCTC_Cut(self, muon, tau):
	if np.signbit(muon.charge) != np.signbit(tau.charge): return True
	else: return False

    def TightID_TightIso_ScaleFactor(self,eta,pt,yeartag):
        ptn=0
        etan=0
        if yeartag == 2018: lookup_table = {
	(0.9, 30): 0.974,
        (1.2, 30): 0.961,
        (2.1, 30): 1.008,
        (2.4, 30): 0.976,

	(0.9, 40): 0.980,
        (1.2, 40): 0.972,
        (2.1, 40): 1.008,
        (2.4, 40): 0.999,

        (0.9, 50): 0.980,
        (1.2, 50): 0.974,
        (2.1, 50): 1.006,
        (2.4, 50): 1.010,

        (0.9, 60): 0.980,
        (1.2, 60): 0.974,
        (2.1, 60): 1.005,
        (2.4, 60): 1.013,

        (0.9, 120): 0.977,
        (1.2, 120): 0.970,
        (2.1, 120): 1.002,
        (2.4, 120): 1.013,

        (0.9, 200): 0.974,
        (1.2, 200): 0.963,
        (2.1, 200): 1.004,
        (2.4, 200): 1.025

	}
        if abs(eta)< 0.9 : etan = 0.9
        if 0.9 <= abs(eta)< 1.2 : etan = 1.2
        if 1.2 <= abs(eta)< 2.1 : etan = 2.1
        if abs(eta)> 2.1 : etan = 2.4

        if pt <30: ptn= 30
        if 30<= pt <40: ptn= 40
        if 40<= pt <50: ptn= 50
        if 50<= pt <60: ptn= 60
        if 60<= pt <120: ptn= 120
        if pt >120: ptn= 200

        scale_factor= lookup_table.get((etan,ptn))
        return scale_factor

    def MediumID_MediumIso_ScaleFactor(self,eta,pt,yeartag):
        ptn=0
        etan=0
        if yeartag == 2018: lookup_table = {
        (0.9, 30): 0.970,
        (1.2, 30): 0.955,
        (2.1, 30): 1.007,
        (2.4, 30): 0.984,

        (0.9, 40): 0.975,
        (1.2, 40): 0.966,
        (2.1, 40): 1.005,
        (2.4, 40): 1.006,

        (0.9, 50): 0.976,
        (1.2, 50): 0.968,
        (2.1, 50): 1.004,
        (2.4, 50): 1.016,

        (0.9, 60): 0.976,
        (1.2, 60): 0.968,
        (2.1, 60): 1.004,
        (2.4, 60): 1.016,

        (0.9, 120): 0.974,
        (1.2, 120): 0.963,
        (2.1, 120): 0.999,
        (2.4, 120): 1.017,

        (0.9, 200): 0.972,
        (1.2, 200): 0.955,
        (2.1, 200): 0.999,
        (2.4, 200): 1.026
        }

	if yeartag == 2017: lookup_table = {
        (0.9, 30): 0.967,
        (1.2, 30): 0.949,
        (2.1, 30): 1.008,
        (2.4, 30): 0.864,

        (0.9, 40): 0.971,
        (1.2, 40): 0.946,
        (2.1, 40): 1.015,
        (2.4, 40): 0.947,

        (0.9, 50): 0.970,
        (1.2, 50): 0.946,
        (2.1, 50): 1.014,
        (2.4, 50): 0.990,

        (0.9, 60): 0.969,
        (1.2, 60): 0.944,
        (2.1, 60): 1.011,
        (2.4, 60): 1.004,

        (0.9, 120): 0.966,
        (1.2, 120): 0.937,
        (2.1, 120): 1.008,
        (2.4, 120): 1.010,

        (0.9, 200): 0.965,
        (1.2, 200): 0.930,
        (2.1, 200): 0.998,
        (2.4, 200): 1.039
        }

	if yeartag == 2016: lookup_table = {
        (0.9, 30): 0.955,
        (1.2, 30): 0.949,
        (2.1, 30): 0.984,
        (2.4, 30): 0.928,

        (0.9, 40): 0.966,
        (1.2, 40): 0.962,
        (2.1, 40): 0.995,
        (2.4, 40): 0.961,

        (0.9, 50): 0.972,
        (1.2, 50): 0.966,
        (2.1, 50): 1.000,
        (2.4, 50): 0.970,

        (0.9, 60): 0.973,
        (1.2, 60): 0.967,
        (2.1, 60): 1.000,
        (2.4, 60): 0.977,

        (0.9, 120): 0.973,
        (1.2, 120): 0.964,
        (2.1, 120): 0.998,
        (2.4, 120): 0.978,

        (0.9, 200): 0.967,
        (1.2, 200): 0.957,
        (2.1, 200): 0.996,
        (2.4, 200): 0.978
        }




        if abs(eta)< 0.9 : etan = 0.9
        if 0.9 <= abs(eta)< 1.2 : etan = 1.2
        if 1.2 <= abs(eta)< 2.1 : etan = 2.1
        if abs(eta)> 2.1 : etan = 2.4

        if pt <30: ptn= 30
        if 30<= pt <40: ptn= 40
        if 40<= pt <50: ptn= 50
        if 50<= pt <60: ptn= 60
        if 60<= pt <120: ptn= 120
        if pt >120: ptn= 200

        scale_factor= lookup_table.get((etan,ptn))
        return scale_factor


    def Highpt_Trigger_Eff(self,eta,pt,yeartag):
	ptn=0
	etan=0
	if yeartag == 2018: lookup_table = {
        (0.9, 55): 0.980,
        (1.2, 55): 0.982,
        (2.1, 55): 0.983,
        (2.4, 55): 0.987,

	(0.9, 60): 0.982,
        (1.2, 60): 0.966,
        (2.1, 60): 1.012,
        (2.4, 60): 0.952,

	(0.9, 120): 0.976,
        (1.2, 120): 0.967,
        (2.1, 120): 1.002,
        (2.4, 120): 0.999,

	(0.9, 200): 0.980,
        (1.2, 200): 0.964,
        (2.1, 200): 1.003,
        (2.4, 200): 1.003,

	(0.9, 300): 0.978,
        (1.2, 300): 0.960,
        (2.1, 300): 1.009,
        (2.4, 300): 1.024,

	(0.9, 500): 0.973,
        (1.2, 500): 0.985,
        (2.1, 500): 0.993,
        (2.4, 500): 0.993,

	(0.9, 501): 0.957,
        (1.2, 501): 0.988,
        (2.1, 501): 1.062,
        (2.4, 501): 1.070
	}

	if yeartag == 2017: lookup_table = {
        (0.9, 55): 0.968,
        (1.2, 55): 0.929,
        (2.1, 55): 0.980,
        (2.4, 55): 0.869,

        (0.9, 60): 0.968,
        (1.2, 60): 0.953,
        (2.1, 60): 0.981,
        (2.4, 60): 0.864,

        (0.9, 120): 0.969,
        (1.2, 120): 0.945,
        (2.1, 120): 0.985,
        (2.4, 120): 0.926,

        (0.9, 200): 0.967,
        (1.2, 200): 0.946,
        (2.1, 200): 0.991,
        (2.4, 200): 0.992,

        (0.9, 300): 0.964,
	    (1.2, 300): 0.948,
        (2.1, 300): 0.986,
        (2.4, 300): 1.023,

        (0.9, 500): 0.963,
        (1.2, 500): 0.939,
        (2.1, 500): 1.002,
        (2.4, 500): 0.905,

        (0.9, 501): 0.951,
        (1.2, 501): 0.896,
        (2.1, 501): 0.922,
        (2.4, 501): 0.424
	}

	if yeartag == 2016: lookup_table = {
        (0.9, 55): 0.980,
        (1.2, 55): 0.938,
        (2.1, 55): 0.996,
        (2.4, 55): 0.934,

        (0.9, 60): 0.983,
        (1.2, 60): 0.958,
        (2.1, 60): 1.005,
        (2.4, 60): 0.954,

        (0.9, 120): 0.980,
        (1.2, 120): 0.956,
        (2.1, 120): 0.990,
        (2.4, 120): 0.948,

        (0.9, 200): 0.980,
        (1.2, 200): 0.948,
        (2.1, 200): 0.988,
        (2.4, 200): 0.935,

        (0.9, 300): 0.980,
	(1.2, 300): 0.927,
        (2.1, 300): 0.977,
        (2.4, 300): 0.888,

        (0.9, 500): 0.962,
        (1.2, 500): 0.911,
        (2.1, 500): 1.010,
        (2.4, 500): 0.953,

        (0.9, 501): 0.983,
        (1.2, 501): 0.924,
        (2.1, 501): 1.020,
        (2.4, 501): 1.134
	}


	if abs(eta)< 0.9 : etan = 0.9
        if 0.9 <= abs(eta)< 1.2 : etan = 1.2
        if 1.2 <= abs(eta)< 2.1 : etan = 2.1
        if abs(eta)> 2.1 : etan = 2.4

        if pt <55: ptn= 55
        if 55<= pt <60: ptn= 60
        if 60<= pt <120: ptn= 120
        if 120<= pt <200: ptn= 200
        if 200<= pt <300: ptn= 300
	if 300<= pt <500: ptn =500
        if pt >500: ptn= 501

	scale_factor= lookup_table.get((etan,ptn))
        return scale_factor
