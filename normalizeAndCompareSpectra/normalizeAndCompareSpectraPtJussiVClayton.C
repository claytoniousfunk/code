#include "/home/clayton/Analysis/code/myProcesses/hiforest/plugin/eventMap_hiForest.h"
#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TLatex.h>
#include <TCut.h>
#include "TDatime.h"
#include <vector>
#include "TCanvas.h"
#include <array>

TH1D *pythia_pt, *data_pt, *r;

double func_temp(double *x, double *par){
    double pt = x[0];
    //int bin = r->FindBin(pt);
    //double fitval = par[0] + par[1]*r->GetBinContent(bin) + par[2]*TMath::Exp(-1*par[3]*r->GetBinContent(bin));
    //double fitval = par[0] + par[1]*pt + par[2]*TMath::Exp(par[3]*pt);
    double fitval = par[0]*TMath::Exp(par[1]*pt) + par[2]*TMath::Exp(par[3]*pt);
    return fitval;
}



void normalizeAndCompareSpectraPtJussiVClayton(){

	TFile* f_clayton =TFile::Open("/home/clayton/Analysis/data/ppskim/pythia_skim_newWeights_muptcut_10.root"); 
	TFile* f_jussi =TFile::Open("/home/clayton/Analysis/data/jussiPPSpectra/ppMC2017_RecoJets_Pythia8_pfJets_wtaAxis_jetWeight_JECv4_processed_2019-11-21.root"); 


	
	
  // pt
    TH1D *clayton_pt, *jussi_pt;
	f_clayton->GetObject("h_jetpt",clayton_pt);
	f_jussi->GetObject("anyJet/anyJetPt_C0",jussi_pt);
    
    
	
  // normalize

  clayton_pt->Scale(1.0/clayton_pt->Integral());
  jussi_pt->Scale(1.0/jussi_pt->Integral());

  double pt_axis[16] = {50.0,60.0,70.0,80.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,250.0,300.0,350.0,400.0,500.0};
  
  TH1D *clayton_pt_rebin = (TH1D*) clayton_pt->Rebin(15,"clayton_pt_rebin",pt_axis);
  TH1D *jussi_pt_rebin = (TH1D*) jussi_pt->Rebin(15,"jussi_pt_rebin",pt_axis);
  TH1D *clayton_pt_rebin_xnorm = (TH1D*) clayton_pt_rebin->Clone("clayton_pt_rebin_xnorm");
  TH1D *jussi_pt_rebin_xnorm = (TH1D*) jussi_pt_rebin->Clone("jussi_pt_rebin_xnorm");
  
  int Nc = clayton_pt_rebin_xnorm->GetSize();
  for(int i=0;i<Nc;i++){
      double xc = clayton_pt_rebin_xnorm->GetBinWidth(i);
      double yc = clayton_pt_rebin_xnorm->GetBinContent(i);
      double zc = clayton_pt_rebin_xnorm->GetBinError(i);
      if(xc!=0){
          clayton_pt_rebin_xnorm->SetBinContent(i,yc/xc);
          clayton_pt_rebin_xnorm->SetBinError(i,zc/xc);
      }
  }
  TF1 *fit_fxn = new TF1("fit_fxn","[0]*exp([1]*x) + [2]*exp([3]*x)",50.0,500.0);
  
  double a0 = 5.65463e-01;
  double e0 = 6.76497e-03;
  double a1 = -2.92413e-03;
  double e1 = 4.93041e-05;
  double a2 = 2.96635e+00;
  double e2 = 1.94034e-02;
  double a3 = -2.59675e-02;
  double e3 = 1.98715e-04;

  fit_fxn->SetParameters(a0,a1,a2,a3);

  TH1D *clayton_pt_rebin_xnorm_fitscaled = (TH1D*) clayton_pt_rebin_xnorm->Clone("clayton_pt_rebin_xnorm_fitscaled");
  for(int k=0;k<Nc;k++){
        double pteval = clayton_pt_rebin_xnorm_fitscaled->GetXaxis()->GetBinCenter(k);
        double fitscale = fit_fxn->Eval(pteval);
        double oldval = clayton_pt_rebin_xnorm_fitscaled->GetBinContent(k);
        double olderror = clayton_pt_rebin_xnorm_fitscaled->GetBinError(k);
        double newval = oldval/fitscale;
        double newerror = newval*TMath::Sqrt(TMath::Power(olderror/oldval,2) + TMath::Power(e0/a0,2) + TMath::Power(e1/a1,2) + TMath::Power(e2/a2,2) + TMath::Power(e3/a3,2));
        clayton_pt_rebin_xnorm_fitscaled->SetBinContent(k,newval);
        clayton_pt_rebin_xnorm_fitscaled->SetBinError(k,newerror);
  }





  int Nj = jussi_pt_rebin_xnorm->GetSize();
  for(int j=0;j<Nj;j++){
      double xj = jussi_pt_rebin_xnorm->GetBinWidth(j);
      double yj = jussi_pt_rebin_xnorm->GetBinContent(j);
      double zj = jussi_pt_rebin_xnorm->GetBinError(j);
      if(xj!=0){
          jussi_pt_rebin_xnorm->SetBinContent(j,yj/xj);
          jussi_pt_rebin_xnorm->SetBinError(j,zj/xj);
      }
  }
  

//////////////////////////////////////////////////////////////////////////////////    spectrum plot    /////////////////////////////////////////////////////////////////////
    TCanvas *c1 = new TCanvas("c1","c1",500,500);

// Define the Canvas
    
      TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
      
      c2->cd();
      
      //c2->SetLogy();
      TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.3, 1., 1.);
      //pad1->SetGridx();         // Vertical grid
      pad2->SetLeftMargin(0.15);
      pad2->SetLogy();
      pad2->Draw();             // Draw the upper pad: pad1
      pad2->cd();
      

      // jet pT 
    
    
      Double_t msize = 1.0;

      //clayton_pt_rebin_xnorm_fitscaled->SetMarkerStyle(24);
      //clayton_pt_rebin_xnorm_fitscaled->SetMarkerColor(kBlue);
      //clayton_pt_rebin_xnorm_fitscaled->SetMarkerSize(msize);
      //clayton_pt_rebin_xnorm_fitscaled->SetFillColorAlpha(kBlue,0.7);
      //clayton_pt_rebin_xnorm_fitscaled->SetStats(0);
      //clayton_pt_rebin_xnorm_fitscaled->SetTitle("");
      clayton_pt_rebin_xnorm->SetMarkerStyle(24);
      clayton_pt_rebin_xnorm->SetMarkerColor(kBlue);
      clayton_pt_rebin_xnorm->SetMarkerSize(msize);
      clayton_pt_rebin_xnorm->SetFillColorAlpha(kBlue,0.7);
      clayton_pt_rebin_xnorm->SetStats(0);
      clayton_pt_rebin_xnorm->SetTitle("");


      jussi_pt_rebin_xnorm->SetMarkerStyle(24);
      jussi_pt_rebin_xnorm->SetMarkerColor(kRed);
      jussi_pt_rebin_xnorm->SetMarkerSize(msize);
      jussi_pt_rebin_xnorm->SetFillColorAlpha(kRed,0.7);
      jussi_pt_rebin_xnorm->SetStats(0);
      jussi_pt_rebin_xnorm->SetTitle("");   


      //clayton_pt_rebin_xnorm_fitscaled->Draw("e2");
      clayton_pt_rebin_xnorm->Draw("e2");
      jussi_pt_rebin_xnorm->Draw("e2 same");
      
      

      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      //legend->AddEntry(clayton_pt_rebin_xnorm_fitscaled,"Clayton","p");
      legend->AddEntry(clayton_pt_rebin_xnorm,"Clayton","p");
      legend->AddEntry(jussi_pt_rebin_xnorm,"Jussi","p");
      legend->SetBorderSize(0);
      legend->Draw();
      

     
      //clayton_pt_rebin_xnorm_fitscaled->GetYaxis()->SetTitle("1/N^{jet} dN^{jet}/dp_{T}");
      //clayton_pt_rebin_xnorm_fitscaled->GetYaxis()->SetTitleSize(0.04);


      clayton_pt_rebin_xnorm->GetYaxis()->SetTitle("1/N^{jet} dN^{jet}/dp_{T}");
      clayton_pt_rebin_xnorm->GetYaxis()->SetTitleSize(0.04);
      clayton_pt_rebin_xnorm->SetMinimum(1e-7);
      clayton_pt_rebin_xnorm->SetMaximum(1e-1);
      
	


///////////////////////////////////////////////////////////////////////////  ratio plot  ////////////////////////////////////////////////////////////////////////////


    c2->cd();
    TPad *pad1 = new TPad("pad1","pad1",0.,0.0,1.,0.3);
    pad1->SetBottomMargin(0.3);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    //r = (TH1D*)clayton_pt_rebin_xnorm_fitscaled->Clone("r");
    r = (TH1D*)clayton_pt_rebin_xnorm->Clone("r");
    r->Divide(jussi_pt_rebin_xnorm);
    r->SetMarkerStyle(24);
    r->SetMarkerColor(kBlack);
    r->SetMarkerSize(msize);
    r->SetFillColorAlpha(kBlack,0.7);
    r->SetStats(0);
    r->SetMinimum(0.0);
    r->SetMaximum(5.0);
    r->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
    r->GetXaxis()->SetTitleSize(0.09);
    r->GetYaxis()->SetTitle("Clayton / Jussi");
    r->GetYaxis()->SetTitleSize(0.08);
    r->SetTitle("");
    r->Draw("e2");

    c2->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/pythiaDataPtCompareJussiVClayton.pdf");





   

}
