void plot(const char* filename_sig, const char* filename_bkg) {
   TFile *fsig = new TFile(filename_sig);
   TFile *fbkg = new TFile(filename_bkg);
   const char* varname[18] = {"isGoodMuon", "highPurity", "TrkMuArb", "TMOneStaTight", "nPixValHits",
      "nMuValHits", "nTrkHits", "normChi2_inner", "normChi2_global", "nPixWMea",
      "nTrkWMea", "StationsMatched", "dxy", "dxyErr", "dz",
      "dzErr", "ptErr_inner", "ptErr_global"};

   const double scale_sig = 0.04;

   for (int i=0; i<18; i++) {
      TCanvas *c1 = new TCanvas();
      TH1F *hsig = (TH1F*) fsig->Get(Form("h%s_sig",varname[i]));
      TH1F *hbkg = (TH1F*) fbkg->Get(Form("h%s_bkg",varname[i]));
      TH1F *hbkg2 = (TH1F*) fsig->Get(Form("h%s_bkg",varname[i]));

      double nbkg = hbkg->GetMaximum();
      double nbkg2 = hbkg2->GetMaximum();
      double nsig = scale_sig*(hsig->GetMaximum() - 0.5*nbkg2);

      int nbins = hsig->GetNbinsX();
      double binmin = hsig->GetXaxis()->GetXmin();
      double binmax = hsig->GetXaxis()->GetXmax();

      TH1F *hsigeff = new TH1F("hsigeff",Form("hsigeff;%s;Signal efficiency",varname[i]),nbins,binmin,binmax);
      hsigeff->SetLineColor(kRed);
      TH1F *hbkgeff = new TH1F("hbkgeff",Form("hbkgeff;%s;1 - Background efficiency",varname[i]),nbins,binmin,binmax);
      hbkgeff->SetLineColor(kBlue);
      TH1F *hsignif = new TH1F("hsignif",Form("hsignif;%s;S / #sqrt{S+B}",varname[i]),nbins,binmin,binmax);
      hsignif->SetLineColor(kGreen);

      for (int j=1; j<nbins+1; j++) {
         double nbkgj = hbkg->GetBinContent(j);
         double nbkgj2 = hbkg2->GetBinContent(j);
         double nsigj = scale_sig*(hsig->GetBinContent(j) - 0.5*nbkgj2);
         // cout << nsigj << " " << nbkgj << " " << nsigj/sqrt(nsigj+nbkgj) << endl;
         hsigeff->SetBinContent(j,nsigj/nsig);
         hbkgeff->SetBinContent(j,1.-nbkgj/nbkg);
         hsignif->SetBinContent(j,nsigj+nbkgj>0 ? nsigj/sqrt(nsigj+nbkgj) : 0);
      }

      hsigeff->GetYaxis()->SetRangeUser(0,1.05);
      hsigeff->Draw("");
      hbkgeff->Draw("same");
      c1->Update();

      // scale hsignif to the pad coordinates
      Float_t rightmax = 1.1*hsignif->GetMaximum();
      Float_t rightmin = 0.98*hsignif->GetMaximum();
      Float_t scale = gPad->GetUymax()/rightmax;
      hsignif->Scale(scale);
      hsignif->Draw("same");
      // draw an axis on the right side
      TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
            gPad->GetUxmax(), gPad->GetUymax(),rightmin,rightmax,510,"+L");
      axis->SetLineColor(kGreen);
      axis->SetLabelColor(kGreen);
      axis->Draw();

      TLegend *tleg = new TLegend(0.5,0.2,0.9,0.6);
      tleg->SetBorderSize(0);
      tleg->AddEntry(hsigeff,"Signal efficiency","lp");
      tleg->AddEntry(hbkgeff,"1 - Background efficiency","lp");
      tleg->AddEntry(hsignif,"S / #sqrt{S+B}","lp");
      tleg->Draw();

      c1->Update();
      c1->SaveAs(Form("%s.pdf",varname[i]));
   }
}
