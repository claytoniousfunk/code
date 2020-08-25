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

void normalizeAndCompareSpectra(){

	TFile* f_pythia =TFile::Open("/home/clayton/Analysis/data/ppskim/clayton_pp_mc_skim_muptcut_10.root"); 
	TFile* f_data =TFile::Open("/home/clayton/Analysis/code/makeData/V2/rootFiles/makeDataV2_pp_data_muptcut_10.root"); 


	
	TH1D *hp, *hd;
  // pt
  /*
	f_pythia->GetObject("h_jetpt",hp);
	f_data->GetObject("h_jetpt",hd);
  */
  // eta
  f_pythia->GetObject("h_jeteta",hp);
	f_data->GetObject("h_jeteta",hd);

  // normalize

  hp->Scale(1.0/hp->Integral());
  hd->Scale(1.0/hd->Integral());

  double pt_axis[16] = {50.0,60.0,70.0,80.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,250.0,300.0,350.0,400.0,500.0};
  double eta_axis[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
  
  



//////////////////////////////////////////////////////////////////////////////////    spectrum plot    /////////////////////////////////////////////////////////////////////
    

// Define the Canvas
    
      TCanvas *c2 = new TCanvas("c2", "canvas", 800, 600);
      c2->SetLeftMargin(0.6);
      c2->cd();
      //c2->SetLogy();
      TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.3, 1., 1.);
      //pad1->SetGridx();         // Vertical grid
      
      pad2->Draw();             // Draw the upper pad: pad1
      pad2->cd();
      

      // jet pT 
      TH1D *hh1 = (TH1D*)hd->Clone("hh1"); // pythia
      TH1D *hh2a = (TH1D*)hp->Clone("hh2a"); // data
      TH1D *h1 = (TH1D*)hh1->Rebin(30,"",eta_axis);
      TH1D *h2a = (TH1D*)hh2a->Rebin(30,"",eta_axis);
      
      
      Double_t msize = 1.2;

      //h1 -> SetMinimum(5*pow(10,-5));
      //h1 -> SetMaximum(1*pow(10,0));

      h1 -> SetMarkerStyle(8);
      h1 -> SetMarkerColor(kBlue);
      h1 -> SetMarkerSize(msize);
      h1 -> SetFillColorAlpha(kBlue,0.7);
      h1 -> SetStats(0);
      

      h2a -> SetMarkerStyle(8);
      h2a -> SetMarkerColor(kRed);
      h2a -> SetMarkerSize(msize);
      h2a -> SetFillColorAlpha(kRed,0.7);
      h2a -> SetStats(0);
    

      h1->Draw("e2");
      h2a->Draw("e2 same");
      

      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      legend->AddEntry(h1,"Pythia","fp");
      legend->AddEntry(h2a,"Data","fp");
      legend->Draw();
      

      h1->GetXaxis()->SetTitle("Jet #eta");
      h1->GetYaxis()->SetTitle("(N^{jet})^{-1} dN^{jet}/d#eta");


      h1->SetTitle("");
      c2->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/pythiaDataPtSpectra.pdf");
	


///////////////////////////////////////////////////////////////////////////  ratio plot  ////////////////////////////////////////////////////////////////////////////


    c2->cd();
    TPad *pad1 = new TPad("pad1","pad1",0.,0.0,1.,0.3);
    //pad1->SetGridx();
    pad1->Draw();
    pad1->cd();
    // pt
    /*
    TH1D *rr = (TH1D*)hp->Clone("rr");
    rr->Divide(hd);
    TH1D *r = (TH1D*)rr->Rebin(15,"",pt_axis);
    */
    // eta
    TH1D *rr = (TH1D*)h1->Clone("rr");
    rr->Divide(h2a);
    TH1D *r = (TH1D*)rr->Rebin(30,"",eta_axis);
    
    
    // normalize r by binwidth
  /*
    Double_t Nbin = r->GetSize();

    for (int i=0;i<Nbin;i++){

      Double_t xp = r->GetBinWidth(i);
      Double_t yp = r->GetBinContent(i);
      Double_t zp = r->GetBinError(i);
      Double_t rp = 0.0;
      Double_t ep = 0.0;
      
      if(xp!=0){

        rp = yp/xp;    
        ep = zp/xp;
        r->SetBinContent(i,rp); 
        r->SetBinError(i,ep);
 
               
      }
      
    }
    */
    r->SetMarkerStyle(8);
    r->SetMarkerColor(kBlack);
    r->SetMarkerSize(1.2);
    r->SetFillColorAlpha(kBlack,0.7);
    r->SetStats(0);
    r->Draw("e2");
    //r->GetXaxis()->SetTitle("jet p_{T} (GeV/c)");
    r->GetXaxis()->SetTitle("jet #eta");
    //r->GetYaxis()->SetTitle("Pythia / data");
    r->SetTitle("Pythia-data #eta ratio");
    r->SetMinimum(0.8);
    r->SetMaximum(1.2);
    //c1->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/pythiaDataPtRatio.pdf");
    //c1->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/pythiaDataEtaRatio.pdf");














}