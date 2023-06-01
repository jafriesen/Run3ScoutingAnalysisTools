# Run3ScoutingAnalysisTools

#### Setup
```
cmsrel CMSSW_12_4_2
cd CMSSW_12_4_2/src
cmsenv
git cms-init
git clone  https://github.com/jafriesen/Run3ScoutingAnalysisTools.git
scram b -j 8
```

#### Run on the scouting dataset
Run a basic example on a file of the Run 3 Scouting dataset:
```
voms-proxy-init --voms cms --valid 168:00
cmsRun tree.py
```
Run a basic example on a file of the Run 2 Scouting dataset:
```
voms-proxy-init --voms cms --valid 168:00
cmsRun treeRun2.py
```
