
Use gfal-ls to create input lists:
```
gfal-ls 'davs://xrootd.cmsaf.mit.edu:1094/store/user/jfriesen/ScoutingCaloMuon/scoutingTreeRun2_23May2023/230523_224825/0000'>list0000.txt
gfal-ls 'davs://xrootd.cmsaf.mit.edu:1094/store/user/jfriesen/ScoutingCaloMuon/scoutingTreeRun2_23May2023/230523_224825/0001'>list0001.txt
```
Use submitFillHistogram.py to create Lxy histograms:
```
mkdir /eos/user/j/jfriesen/fillHistogram/histoFullRun3
python3 submitFillHistogram.py -i 'root://cmsxrootd.fnal.gov//store/user/jfriesen/ScoutingPFRun3/scoutingTreeRun3_23May2023/230523_182013' -o h_lxy -e /eos/user/j/jfriesen/fillHistogram/histoFullRun3 -t histoFullRun3_26May2023 -n 400 -s /afs/cern.ch/user/j/jfriesen/CMSSW_12_4_2/src/Run3ScoutingAnalysisTools/FillHistogram/fillHistogram.py -r run3 -c full
mkdir /eos/user/j/jfriesen/fillHistogram/histoFullRun2
python3 submitFillHistogram.py -i 'root://cmsxrootd.fnal.gov//store/user/jfriesen/ScoutingCaloMuon/scoutingTreeRun2_23May2023/230523_224825' -o h_lxy -e /eos/user/j/jfriesen/fillHistogram/histoFullRun2 -t histoFullRun2_26May2023 -n 400 -s /afs/cern.ch/user/j/jfriesen/CMSSW_12_4_2/src/Run3ScoutingAnalysisTools/FillHistogram/fillHistogram.py -r run2 -c full
```
Use submitFillHistogram.py to create mass histograms:
```
mkdir /eos/user/j/jfriesen/fillHistogram/massRun2lxy02to1
python3 submitFillHistogram.py -i 'root://cmsxrootd.fnal.gov//store/user/jfriesen/ScoutingCaloMuon/scoutingTreeRun2_23May2023/230523_224825' -o h_mass -e /eos/user/j/jfriesen/fillHistogram/massRun2lxy02to1 -t massRun2lxy02to1_2June2023 -n 400 -s /afs/cern.ch/user/j/jfriesen/CMSSW_12_4_2/src/Run3ScoutingAnalysisTools/FillHistogram/fillMassHistogram.py -r run2 -c vtx --min 0.2 --max 1.0
mkdir /eos/user/j/jfriesen/fillHistogram/massRun3lxy02to1
python3 submitFillHistogram.py -i 'root://cmsxrootd.fnal.gov//store/user/jfriesen/ScoutingPFRun3/scoutingTreeRun3_23May2023/230523_182013' -o h_mass -e /eos/user/j/jfriesen/fillHistogram/massRun3lxy02to1 -t massRun3lxy02to1_2June2023 -n 400 -s /afs/cern.ch/user/j/jfriesen/CMSSW_12_4_2/src/Run3ScoutingAnalysisTools/FillHistogram/fillMassHistogram.py -r run3 -c vtx --min 0.2 --max 1.0
```
