
./Make.sh ZTT_XSection.cc

# ./ZTT_XSection.exe output/SingleMu.root root://131.225.204.161:1094//store/user/cmsdas/2019/long_exercises/ZTauTau/SingleMuon.root
# ./ZTT_XSection.exe output/SingleEl.root root://131.225.204.161:1094//store/user/cmsdas/2019/long_exercises/ZTauTau/SingleElectron.root

./ZTT_XSection.exe output/DYJetsToTauTau.root root://131.225.204.161:1094//store/user/cmsdas/2019/long_exercises/ZTauTau/DYJetsToLL_M-50_Inc.root
./ZTT_XSection.exe output/DYJetsToLL.root root://131.225.204.161:1094//store/user/cmsdas/2019/long_exercises/ZTauTau/DYJetsToLL_M-50_Inc.root

# ./ZTT_XSection.exe output/TTJets.root root://131.225.204.161:1094//store/user/cmsdas/2019/long_exercises/ZTauTau/TTbar.root
# ./ZTT_XSection.exe output/WJetsToLNu.root root://131.225.204.161:1094//store/user/cmsdas/2019/long_exercises/ZTauTau/WJetsToLNu_Inc.root
# ./ZTT_XSection.exe output/WW.root root://131.225.204.161:1094//store/user/cmsdas/2019/long_exercises/ZTauTau/WW.root
# ./ZTT_XSection.exe output/WZ.root root://131.225.204.161:1094//store/user/cmsdas/2019/long_exercises/ZTauTau/WZ.root
# ./ZTT_XSection.exe output/ZZ.root root://131.225.204.161:1094//store/user/cmsdas/2019/long_exercises/ZTauTau/ZZ.root

cd output
python xs_calculator_prefit.py etau
python xs_calculator_prefit.py mutau
cd ..