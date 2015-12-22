#define cutoptim_cxx
#include "cutoptim.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TString.h>

#include <iostream>
#include <map>

using namespace std;

void cutoptim::Loop()
{
//   In a ROOT session, you can do:
//      root> .L cutoptim.C
//      root> cutoptim t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *f = new TFile("histos.root","RECREATE");

   const int nvar = 18;
   const char* varname[18] = {"isGoodMuon", "highPurity", "TrkMuArb", "TMOneStaTight", "nPixValHits",
      "nMuValHits", "nTrkHits", "normChi2_inner", "normChi2_global", "nPixWMea",
      "nTrkWMea", "StationsMatched", "dxy", "dxyErr", "dz",
      "dzErr", "ptErr_inner", "ptErr_global"};
   const int nbins[18] = {2, 2, 2, 2, 11,
      56, 36, 100, 100, 6,
      19, 7, 100, 100, 100,
      100, 100, 100};
   const double lowedge[18] = {0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0};
   const double hiedge[18] = {2, 2, 2, 2, 11,
      56, 36, 10, 50, 6,
      19, 7, 2, 0.25, 20, 
      5, 0.15, 0.15};
   map<TString, TH1F*> hists_sig, hists_bkg;

   for (int i=0; i<nvar; i++) {
      hists_sig[Form("%s",varname[i])] = new TH1F(Form("h%s_sig",varname[i]),Form("h%s_sig",varname[i]),nbins[i],lowedge[i],hiedge[i]);
      hists_bkg[Form("%s",varname[i])] = new TH1F(Form("h%s_bkg",varname[i]),Form("h%s_bkg",varname[i]),nbins[i],lowedge[i],hiedge[i]);
   }


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // loop on reco dimuons
      for (int i=0; i<Reco_QQ_size; i++) {
         if (Reco_QQ_sign[i]!=0) continue;

         TLorentzVector *tlvmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(i);
         TLorentzVector *tlvmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(i);

         if (!IsAccept(tlvmupl->Pt(), tlvmupl->Eta()) || !IsAccept(tlvmumi->Pt(), tlvmumi->Eta())) continue;

         TLorentzVector *tlvqq = (TLorentzVector*) Reco_QQ_4mom->At(i);
         double mass = tlvqq->M();

         bool issig = (mass>3 && mass<3.2);
         bool isbkg = ((mass>2.7 && mass<2.9) || (mass>3.3 && mass<3.5));

         // fill histos
         for (int ihist=0; ihist<nvar; ihist++) {
            for (int j=1; j<=hists_sig[varname[ihist]]->GetNbinsX(); j++) {
               if (Cut(varname[ihist],i,hists_sig[varname[ihist]]->GetBinLowEdge(j))) {
                  if (issig) hists_sig[varname[ihist]]->Fill(hists_sig[varname[ihist]]->GetBinCenter(j));
                  if (isbkg) hists_bkg[varname[ihist]]->Fill(hists_bkg[varname[ihist]]->GetBinCenter(j));
               }
               // else if (j==1) cout << "mmm... this event fails the minimum cut on " << varname[ihist] <<  " (" << hists_sig[varname[ihist]]->GetBinLowEdge(1) << ")" << endl;
            }
         }
      }
   }

   f->Write();
   f->Close();
}
