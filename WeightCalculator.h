#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TSystem.h"
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"
#include <TLorentzVector.h>
using namespace std;

float XSection(std::string OutName) {
    
    
//    https://docs.google.com/spreadsheets/d/1rWM3AlFKO8IJVaeoQkWZYWwSvicQ1QCXYSzH74QyZqE/edit?alt=json#gid=398123591
    
    if (OutName.compare("WJets") == 0) return 61526;
    else if (OutName.compare("DYJetsToLL") == 0) return 5765.4; // 3 * 1921=57654   //6025.2
    else if (OutName.compare("TTJets") == 0) return 831.76 ;
    
    //Di-boson
    else if (OutName.compare("WW") == 0) return 115.0;
    else if (OutName.compare("WZ") == 0) return 47.13;
    else if (OutName.compare("ZZ") == 0) return 16.523;
    
    else return 0;
}


float weightCalc(TH1F *Histo,std::string outputName) {
    
//    cout<< "outputName is "<<outputName << "  and histoname is " <<Histo->GetName()<<  " Histo->GetBinContent(1)="<<Histo->GetBinContent(1)<< " XSection(wjet)=" <<XSection("WJets")<<"\n";
    
        float luminosity=35900;
    
    size_t isSingleMu = outputName.find("SingleMu");
    size_t isSingleEle = outputName.find("SingleEle");
    size_t isWjet = outputName.find("WJets");
    size_t isZJet = outputName.find("DYJetsToLL");
    size_t isTTbar = outputName.find("TTbar");
    size_t isWW = outputName.find("WW");
    size_t isWZ = outputName.find("WZ");
    size_t isZZ = outputName.find("ZZ");
    
    
    if (isSingleMu != string::npos || isSingleEle!= string::npos)   return 1;
    else if ( isWjet != string::npos ) return (luminosity * XSection("WJets")*1.0) / Histo->GetBinContent(2) ;
    else if ( isZJet != string::npos )   return luminosity * XSection("DYJetsToLL")*1.0 / Histo->GetBinContent(2) ;
    else if ( isTTbar != string::npos )   return luminosity * XSection("TTJets")*1.0 / Histo->GetBinContent(2) ;
    else if ( isWW != string::npos )   return luminosity * XSection("WW")*1.0 / Histo->GetBinContent(2) ;
    else if ( isWZ != string::npos )   return luminosity * XSection("WZ")*1.0 / Histo->GetBinContent(2) ;
    else if ( isZZ != string::npos )   return luminosity * XSection("ZZ")*1.0 / Histo->GetBinContent(2) ;
    else    return 0;
    
}



float TMass_F(float pt3lep, float px3lep, float py3lep, float met, float metPhi) {
    return sqrt(pow(pt3lep + met, 2) - pow(px3lep + met * cos(metPhi), 2) - pow(py3lep + met * sin(metPhi), 2));
}











