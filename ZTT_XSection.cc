////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Compiling the code:   ./Make.sh ZTT_XSection.cc
//   Running the code:     ./ZTT_XSection.exe OutPut.root   Input.root
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TreeReader.h"
#include "WeightCalculator.h"
#include <string> //remove space for successful compilation
#include <ostream> //remove space for successful compilation
#include <ctime>

void initBranch(TTree *tree, std::string name, void *add);
void registerBranches(TTree *tree);


int main(int argc, char** argv) {
  using namespace std;

  std::string out = *(argv + 1);

  cout << "\n\n\n OUTPUT NAME IS:    " << out << endl;     //PRINTING THE OUTPUT FILE NAME
  TFile *fout = TFile::Open(out.c_str(), "RECREATE");
  std::string input = *(argv + 2);
  cout << "\n\n\n INPUT NAME IS:    " << input << endl;     //PRINTING THE INPUT FILE NAME
  TFile * myFile = TFile::Open(input.c_str());
  TH1F * HistoTot = (TH1F*) myFile->Get("hcount");

  bool isDYJetsToLL = false;
  bool isDYJetsToTauTau = false;

  if (out=="output/DYJetsToLL.root"){
    isDYJetsToLL = true;
  } else if (out=="output/DYJetsToTauTau.root"){
    isDYJetsToTauTau = true;
  }


  //global setting for Luminosity
  float LumiWeight = 1;
  if (HistoTot) LumiWeight = weightCalc(HistoTot, input);
  cout << "LumiWeight is " << LumiWeight << "\n";

  //global setting for Pileup Reweighting
  TFile * PUData = new TFile("data_pu/MyDataPileupHistogram2016.root");
  TH1F * HistoPUData = (TH1F *) PUData->Get("pileup");
  HistoPUData->Scale(1.0/HistoPUData->Integral());
  TFile * PUMC= new TFile("data_pu/mcMoriondPU.root");
  TH1F * HistoPUMC= (TH1F *) PUMC->Get("pileup");
  HistoPUMC->Scale(1.0/HistoPUMC->Integral());

  //global setting for reading ntuple trees
  TTree *Run_Tree = (TTree*) myFile->Get("EventTree");
  cout.setf(ios::fixed, ios::floatfield);
  registerBranches(Run_Tree);
  
  //global declare output
  //add the histrograms of muon and tau visible mass (both for opposite sign and same sign pair )
  TH1F *    mutau_visibleMassOS = new TH1F ("mutau_visibleMassOS","mutau_visibleMassOS", 30, 0, 300);
  TH1F *    mutau_visibleMassSS = new TH1F ("mutau_visibleMassSS","mutau_visibleMassSS", 30, 0, 300);
  //add the histrograms of electron and tau visible mass (both for opposite sign and same sign pair )
  TH1F *    etau_visibleMassOS = new TH1F ("etau_visibleMassOS","etau_visibleMassOS", 30, 0, 300);
  TH1F *    etau_visibleMassSS = new TH1F ("etau_visibleMassSS","etau_visibleMassSS", 30, 0, 300);


  //global setting for lepton mass
  float muMass  = 0.10565837;
  float eleMass = 0.000511;


  //loop over events
  Int_t nentries_wtn = (Int_t) Run_Tree->GetEntries();
  cout<<"nentries_wtn is " << nentries_wtn << "\n";
  for ( Int_t i = 0; i < nentries_wtn; i++) {

    Run_Tree->GetEntry(i);

    if (i % 1000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);

    fflush(stdout);

    TLorentzVector Lep4Momentum, Tau4Momentum;

    //0.1 event trigger
    bool PassMuTrigger = (HLTEleMuX >> 19 & 1) == 1; // if (name.find("HLT_IsoMu24_v") != string::npos) bitEleMuX = 19;
    bool PassElTrigger = (HLTEleMuX >> 0  & 1) == 1; // if (name.find("HLT_Ele25_eta2p1_WPTight_Gsf_v")  string::npos) bitEleMuX = 0; 
    
    if ((! PassMuTrigger) && (! PassElTrigger)) continue;

    //0.2 event PU weight
    float PUWeight = 1;
    if (!isData){
      int puNUmmc=int(puTrue->at(0)*10);
      int puNUmdata=int(puTrue->at(0)*10);
      float PUMC_=HistoPUMC->GetBinContent(puNUmmc+1);
      float PUData_=HistoPUData->GetBinContent(puNUmdata+1);
      PUWeight= PUData_/PUMC_;
    }

    // cout<<"event PU weight is "<< PUWeight <<endl;
    // cout<<"nEle is "<<nEle<<", nMu is "<<nMu<<", nTau is "<<nMu<<endl;
    // continue;
    ///////////////////////////////////////////////
    //Important Analysis Loop Will Happen Here!!!//
    ///////////////////////////////////////////////

    //////////////////////
    // 1.preselection
    //////////////////////

    // 1.0 gen particle

    int numGenTau=0;
    int numGenEle=0;

    for (int igen=0;igen < nMC; igen++){
      if ( fabs(mcPID->at(igen)) ==15 && mcMomPID->at(igen)==23 ) numGenTau++;
      if ( fabs(mcPID->at(igen)) ==11 && mcMomPID->at(igen)==23 ) numGenEle++;
    }

    if (isDYJetsToLL && numGenTau>0) continue;
    if (isDYJetsToTauTau && numGenTau<1) continue;


    // 1.1 muon
    std::vector<int>  idxPassMu,idxFailMu;
    if (nMu>0){
      for  (int imu=0 ; imu < nMu; imu++){

        float IsoMu = muPFChIso->at(imu)/muPt->at(imu);
        if ( (muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu) )  > 0.0){
          IsoMu = ( muPFChIso->at(imu)/muPt->at(imu) + muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu))/muPt->at(imu);
        }

        bool selected = true;
        selected = selected && (muPt->at(imu) > 25);
        selected = selected && (fabs(muEta->at(imu)) < 2.1);
        selected = selected && (muIDbit->at(imu) >> 2 & 2); // 2 is tight;
        if (selected && (IsoMu < 0.15)) {
          idxPassMu.push_back(imu);
        }
        if (selected && (IsoMu > 0.15)) {
          idxFailMu.push_back(imu);
        }
      } // End of muon loop
    } 


    // 1.2 electron 
    std::vector<int>  idxPassEle, idxFailEle;
    if (nEle>0){
      for  (int iele=0 ; iele < nEle; iele++){

        float IsoEle=elePFChIso->at(iele)/elePt->at(iele);
        if ( (elePFNeuIso->at(iele) + elePFPhoIso->at(iele) - 0.5* elePFPUIso->at(iele))  > 0.0){
          IsoEle= (elePFChIso->at(iele)/elePt->at(iele) + elePFNeuIso->at(iele) + elePFPhoIso->at(iele) - 0.5* elePFPUIso->at(iele))/elePt->at(iele);
        }
                    
        bool eleMVAId= false;
        if (fabs (eleSCEta->at(iele)) <= 0.8 && eleIDMVA->at(iele) >   0.941  ) eleMVAId= true;
        else if (fabs (eleSCEta->at(iele)) >  0.8 &&fabs (eleSCEta->at(iele)) <=  1.5 && eleIDMVA->at(iele) >   0.899  ) eleMVAId= true;
        else if ( fabs (eleSCEta->at(iele)) >=  1.5 && eleIDMVA->at(iele) >  0.758  ) eleMVAId= true;
        else eleMVAId= false;

        bool selected = true;
        selected = selected && (elePt->at(iele) > 30);
        selected = selected && (fabs(eleEta->at(iele)) < 2.1);
        selected = selected && (eleMVAId);
        selected = selected ;
        if (selected && (IsoEle<0.1) ) {
          idxPassEle.push_back(iele);
        }
        if (selected && (IsoEle>0.1) ) {
          idxFailEle.push_back(iele);
        }
      } // End of Electron Loop
    } 


    // 1.3 taus
    std::vector<int>  idxPassTau;
    if(nTau > 0){
      for  (int itau=0 ; itau < nTau; itau++){
        TLorentzVector tauP4;
        tauP4.SetPtEtaPhiM(tauPt->at(itau),tauEta->at(itau),tauPhi->at(itau),tauMass->at(itau));

        // find if tau is overlap with good electrons
        bool eleOverlay = false;
        if (idxPassEle.size()>0){
          for (int i=0; i<idxPassEle.size(); i++){
            TLorentzVector lepP4;
            lepP4.SetPtEtaPhiM(elePt->at(idxPassEle[i]),eleEta->at(idxPassEle[i]),elePhi->at(idxPassEle[i]),eleMass);
            if (tauP4.DeltaR(lepP4) < 0.5){
              eleOverlay = true;
            }
          }
        }

        // find if tau is overlap with good muons
        bool muOverlay = false;
        if (idxPassMu.size()>0){
          for (int i=0; i<idxPassMu.size(); i++){
            TLorentzVector lepP4;
            lepP4.SetPtEtaPhiM(muPt->at(idxPassMu[i]),muEta->at(idxPassMu[i]),muPhi->at(idxPassMu[i]),muMass);
            if (tauP4.DeltaR(lepP4) < 0.5){
              muOverlay = true;
            }
          }
        }


        bool tauIdIso = false;
        if ( taupfTausDiscriminationByDecayModeFinding->at(itau) > 0.5 
            && tauByTightMuonRejection3->at(itau) > 0
            && tauByMVA6TightElectronRejection->at(itau) > 0
            && tauByTightIsolationMVArun2v1DBoldDMwLT->at(itau) > 0
        ) tauIdIso = true;

        bool selected = true;
        selected = selected && (tauPt->at(itau) > 30);
        selected = selected && (fabs(tauEta->at(itau)) < 2.3);
        selected = selected && (!eleOverlay);
        selected = selected && (!muOverlay);
        selected = selected && (tauIdIso);
        if (selected) {
          idxPassTau.push_back(itau);
        }
      } // End of tau loop
    }

    // 1.4 jet and btag
    int nPassJet = 0;
    int nPassBJet = 0;
    if (nJet>0){
      for (int ijet= 0 ; ijet < nJet ; ijet++){
        TLorentzVector jetP4;
        jetP4.SetPtEtaPhiE(jetPt->at(ijet),jetEta->at(ijet),jetPhi->at(ijet),jetEn->at(ijet));

        // find if tau is overlap with good electrons
        bool eleOverlay = false;
        if (idxPassEle.size()>0){
          for (int i=0; i<idxPassEle.size(); i++){
            TLorentzVector lepP4;
            lepP4.SetPtEtaPhiM(elePt->at(idxPassEle[i]),eleEta->at(idxPassEle[i]),elePhi->at(idxPassEle[i]),eleMass);
            if (jetP4.DeltaR(lepP4) < 0.5){
              eleOverlay = true;
            }
          }
        }

        // find if tau is overlap with good muons
        bool muOverlay  = false;
        if (idxPassMu.size()>0){
          for (int i=0; i<idxPassMu.size(); i++){
            TLorentzVector lepP4;
            lepP4.SetPtEtaPhiM(muPt->at(idxPassMu[i]),muEta->at(idxPassMu[i]),muPhi->at(idxPassMu[i]),muMass);
            if (jetP4.DeltaR(lepP4) < 0.5){
              muOverlay = true;
            }
          }
        }

        // find if tau is overlap with good muons
        bool tauOverlay = false;
        if (idxPassTau.size()>0){
          for (int i=0; i<idxPassTau.size(); i++){
            TLorentzVector lepP4;
            lepP4.SetPtEtaPhiM(tauPt->at(idxPassTau[i]),tauEta->at(idxPassTau[i]),tauPhi->at(idxPassTau[i]),tauMass->at(idxPassTau[i]));
            if (jetP4.DeltaR(lepP4) < 0.5){
              tauOverlay = true;
            }
          }
        }

        bool selected = true;
        selected = selected && (jetPt->at(ijet) > 20);
        selected = selected && (fabs(jetEta->at(ijet)) < 2.5);
        selected = selected && jetPFLooseId->at(ijet)==1;
        selected = selected && (!eleOverlay);
        selected = selected && (!muOverlay);
        selected = selected && (!tauOverlay);
        

        if (selected){
          nPassJet ++;
          if (jetCSV2BJetTags->at(ijet) > 0.8484){
            nPassBJet ++;
          }
        }      
      } // End of jet loop
    }


    //////////////////////
    // 2.analysis selection
    //////////////////////

    //cout<<idxPassEle.size()<<idxPassMu.size()<<idxPassTau.size()<<idxFailMu.size()<<idxPassEle.size()<<nPassBJet<<endl;

    // 2.1 mutau channel
    if (   idxPassEle.size()==0 
        && idxPassMu.size() ==1
        && idxPassTau.size()>=1
        // && idxFailEle.size()==0 
        // && idxFailMu.size() ==0 
        && nPassBJet==0 && PassMuTrigger){


      Lep4Momentum.SetPtEtaPhiM(muPt->at(idxPassMu[0]), muEta->at(idxPassMu[0]), muPhi->at(idxPassMu[0]), muMass);
      Tau4Momentum.SetPtEtaPhiM(tauPt->at(idxPassTau[0]), tauEta->at(idxPassTau[0]), tauPhi->at(idxPassTau[0]), tauMass->at(idxPassTau[0]));
      TLorentzVector Dilepton4Momentum = Lep4Momentum + Tau4Momentum;

      float LepMetTranverseMass = TMass_F(Lep4Momentum.Pt(), Lep4Momentum.Pt()*cos(Lep4Momentum.Phi()),Lep4Momentum.Pt()*sin(Lep4Momentum.Phi()) ,  pfMET, pfMETPhi);
      bool wveto = LepMetTranverseMass<40;
      bool  OS = muCharge->at(idxPassMu[0]) * tauCharge->at(idxPassTau[0]) < 0;
      bool  SS = muCharge->at(idxPassMu[0]) * tauCharge->at(idxPassTau[0]) > 0;

      if (wveto && OS){
        //Check if there is an OS and TMass(mu.MET) < 40 and then fill the weighted histogram as below:
        mutau_visibleMassOS->SetDefaultSumw2();
        mutau_visibleMassOS->Fill(Dilepton4Momentum.M(),LumiWeight*PUWeight);
      } else if (wveto && SS){
        //Check if there is a SS and TMass(mu.MET) < 40 and then fill the weighted histogram as below:
        mutau_visibleMassSS->SetDefaultSumw2();
        mutau_visibleMassSS->Fill(Dilepton4Momentum.M(),LumiWeight*PUWeight);
      }
    }

    // 2.2 etau channel
    else if (  idxPassEle.size()==1 
            && idxPassMu.size() ==0 
            && idxPassTau.size()>=1 
            // && idxFailEle.size()==0 
            // && idxFailMu.size() ==0 
            && nPassBJet==0 && PassElTrigger){



      Lep4Momentum.SetPtEtaPhiM(elePt->at(idxPassEle[0]), eleEta->at(idxPassEle[0]), elePhi->at(idxPassEle[0]), eleMass);
      Tau4Momentum.SetPtEtaPhiM(tauPt->at(idxPassTau[0]), tauEta->at(idxPassTau[0]), tauPhi->at(idxPassTau[0]), tauMass->at(idxPassTau[0]));
      
      if (isDYJetsToLL){
        Tau4Momentum *= 1.1;
        if (abs(Tau4Momentum.Eta())<1.44) {
          PUWeight *= 1.4;
        } else if (abs(Tau4Momentum.Eta())>1.44 && abs(Tau4Momentum.Eta())<2.5){
          PUWeight *= 1.9;
        }
        
      }

      TLorentzVector Dilepton4Momentum = Lep4Momentum + Tau4Momentum;

      float LepMetTranverseMass = TMass_F(Lep4Momentum.Pt(), Lep4Momentum.Pt()*cos(Lep4Momentum.Phi()),Lep4Momentum.Pt()*sin(Lep4Momentum.Phi()) ,  pfMET, pfMETPhi);
      bool wveto = LepMetTranverseMass<40;
      bool  OS = eleCharge->at(idxPassEle[0]) * tauCharge->at(idxPassTau[0]) < 0;
      bool  SS = eleCharge->at(idxPassEle[0]) * tauCharge->at(idxPassTau[0]) > 0;

      if (wveto && OS){
        //Check if there is an OS and TMass(mu.MET) < 40 and then fill the weighted histogram as below:
        etau_visibleMassOS->SetDefaultSumw2();
        etau_visibleMassOS->Fill(Dilepton4Momentum.M(),LumiWeight*PUWeight);
      } else if (wveto && SS){
        //Check if there is a SS and TMass(mu.MET) < 40 and then fill the weighted histogram as below:
        etau_visibleMassSS->SetDefaultSumw2();
        etau_visibleMassSS->Fill(Dilepton4Momentum.M(),LumiWeight*PUWeight);
      }
    }


  } //End Processing all entries
   

  //end of analysis code, close and write histograms/file
  fout->cd();
  mutau_visibleMassOS->Write();
  mutau_visibleMassSS->Write();
  etau_visibleMassOS->Write();
  etau_visibleMassSS->Write();
  fout->Close();
    
}


void initBranch(TTree *tree, std::string name, void *add){
  const char *bname = name.c_str();
  tree->SetBranchStatus(bname, 1);
  tree->SetBranchAddress(bname, add);
}

void registerBranches(TTree *tree) {
  tree->SetBranchStatus("*", 0);


  //########################################   General Info
  initBranch(tree, "isData", &isData);
  initBranch(tree,"run", &run);
  initBranch(tree,"lumis", &lumis);
  initBranch(tree,"event", &event);
  initBranch(tree,"genWeight",&genWeight);
  initBranch(tree,"HLTEleMuX", &HLTEleMuX);
  initBranch(tree,"puTrue", &puTrue);
  initBranch(tree,"nVtx",&nVtx);
        
  //########################################   MC Info
  initBranch(tree,"nMC", &nMC);
  initBranch(tree,"mcPID", &mcPID);
  initBranch(tree,"mcStatus", &mcStatus);
  initBranch(tree,"mcPt", &mcPt );
  initBranch(tree,"mcEta", &mcEta );
  initBranch(tree,"mcPhi", &mcPhi );
  initBranch(tree,"mcE", &mcE );
  initBranch(tree,"mcMass", &mcMass );
  initBranch(tree,"mcMomPID", &mcMomPID );
  initBranch(tree,"mcGMomPID", &mcGMomPID );
        
  //########################################   Tau Info
  initBranch(tree,"nTau", &nTau);
  initBranch(tree,"tauPt"  ,&tauPt);
  initBranch(tree,"tauEta"  ,&tauEta);
  initBranch(tree,"tauPhi"  ,&tauPhi);
  initBranch(tree,"tauMass"  ,&tauMass);
  initBranch(tree,"tauCharge"  ,&tauCharge);

  initBranch(tree,"taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding);

  initBranch(tree,"tauByTightMuonRejection3", &tauByTightMuonRejection3);
  initBranch(tree,"tauByLooseMuonRejection3", &tauByLooseMuonRejection3);

  initBranch(tree,"tauByMVA6TightElectronRejection"  ,&tauByMVA6TightElectronRejection);
  initBranch(tree,"tauByMVA6MediumElectronRejection"  ,&tauByMVA6MediumElectronRejection);
  initBranch(tree,"tauByMVA6LooseElectronRejection", &tauByMVA6LooseElectronRejection);

  initBranch(tree,"tauDxy",&tauDxy);
  initBranch(tree,"tauDecayMode",&tauDecayMode);

  initBranch(tree,"tauByLooseIsolationMVArun2v1DBoldDMwLT",&tauByLooseIsolationMVArun2v1DBoldDMwLT);
  initBranch(tree,"tauByVLooseIsolationMVArun2v1DBoldDMwLT",&tauByVLooseIsolationMVArun2v1DBoldDMwLT);
  initBranch(tree,"tauByTightIsolationMVArun2v1DBoldDMwLT",&tauByTightIsolationMVArun2v1DBoldDMwLT);
        
  //########################################   Mu Info
  initBranch(tree,"nMu", &nMu);
  initBranch(tree,"muPt"  ,&muPt);
  initBranch(tree,"muEta"  ,&muEta);
  initBranch(tree,"muPhi"  ,&muPhi);
  initBranch(tree,"muIsoTrk", &muIsoTrk);
  initBranch(tree,"muCharge",&muCharge);
  initBranch(tree,"muIDbit",&muIDbit);//NEW
  initBranch(tree,"muPFChIso", &muPFChIso);
  initBranch(tree,"muPFPhoIso", &muPFPhoIso);
  initBranch(tree,"muPFNeuIso", &muPFNeuIso);
  initBranch(tree,"muPFPUIso", &muPFPUIso);
  initBranch(tree,"muD0",&muD0);
  initBranch(tree,"muDz",&muDz);
        
  //########################################   Ele Info
  initBranch(tree,"nEle", &nEle);
  initBranch(tree,"elePt"  ,&elePt);
  initBranch(tree,"eleEta"  ,&eleEta);
  initBranch(tree,"elePhi"  ,&elePhi);
  initBranch(tree,"elePFChIso", &elePFChIso);
  initBranch(tree,"eleIDMVA", &eleIDMVA);//NEW
  initBranch(tree,"eleCharge",&eleCharge);
  initBranch(tree,"eleSCEta",&eleSCEta);
  initBranch(tree,"elePFChIso", &elePFChIso);
  initBranch(tree,"elePFPhoIso", &elePFPhoIso);
  initBranch(tree,"elePFNeuIso", &elePFNeuIso);
  initBranch(tree,"elePFPUIso", &elePFPUIso);
  initBranch(tree,"eleD0",&eleD0);
  initBranch(tree,"eleDz",&eleDz);
  initBranch(tree,"eleMissHits", &eleMissHits);
  initBranch(tree,"eleConvVeto", &eleConvVeto);
  initBranch(tree,"eleSCEta", &eleSCEta );
        
  //########################################   Jet Info
  initBranch(tree,"nJet",&nJet);
  initBranch(tree,"jetPt",&jetPt);
  initBranch(tree,"jetEta",&jetEta);
  initBranch(tree,"jetPhi",&jetPhi);
  initBranch(tree,"jetEn",&jetEn);
  initBranch(tree,"jetCSV2BJetTags",&jetCSV2BJetTags);
  initBranch(tree,"jetPFLooseId",&jetPFLooseId);
  initBranch(tree,"jetPUID",&jetPUID);
  initBranch(tree,"jetRawPt",&jetRawPt);
  initBranch(tree,"jetJECUnc",&jetJECUnc);
  initBranch(tree,"jetRawEn",&jetRawEn);
  initBranch(tree,"jetHadFlvr",&jetHadFlvr);
        
  //########################################   MET Info
  initBranch(tree,"pfMET",&pfMET);
  initBranch(tree,"pfMETPhi",&pfMETPhi);
  initBranch(tree,"metFilters",&metFilters);
  initBranch(tree,"genHT",&genHT);
}
