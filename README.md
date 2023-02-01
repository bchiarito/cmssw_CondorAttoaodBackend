### interactive testing
```
python nano_postproc.py infiles.txt . --drop my_ana_drop.txt --filter="one_hpid_photon" --add_recophi HPID -n=10 --data
python metadata_create.py report.txt NANOAOD_TwoProng_986_Skim.root datasetname 0
```

### setup instructions

```
cmsrel CMSSW_10_6_20
cd CMSSW_10_6_20/src
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
cd PhysicsTools/NanoAODTools/python
git clone https://github.com/bchiarito/cmssw_temp_attoframework_backend.git fmk_atto
cd ../../..
scram b
```
