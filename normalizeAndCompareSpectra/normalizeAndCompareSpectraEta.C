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

void normalizeAndCompareSpectraEta(){

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
    
      TCanvas *c2 = new TCanvas("c2", "canvas", 600, 600);
      c2->SetLeftMargin(0.6);
      c2->cd();
      //c2->SetLogy();
      TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.3, 1., 1.);
      //pad1->SetGridx();         // Vertical grid
      pad2->SetLeftMargin(0.15);
      
      pad2->Draw();             // Draw the upper pad: pad1
      pad2->cd();
      

      // jet pT 
      TH1D *hh1 = (TH1D*)hd->Clone("hh1"); // pythia
      TH1D *hh2 = (TH1D*)hp->Clone("hh2"); // data
      TH1D *h1 = (TH1D*)hh1->Rebin(30,"",eta_axis);
      TH1D *h2 = (TH1D*)hh2->Rebin(30,"",eta_axis);

       // normalize h1 by bin width
  
    Double_t Nbin1 = h1->GetSize();

    for (int i=0;i<Nbin1;i++){

      Double_t x1 = h1->GetBinWidth(i);
      Double_t y1 = h1->GetBinContent(i);
      Double_t z1 = h1->GetBinError(i);
      Double_t r1 = 0.0;
      Double_t e1 = 0.0;
      
      if(x1!=0){

        r1 = y1/x1;    
        e1 = z1/x1;
        h1->SetBinContent(i,r1); 
        h1->SetBinError(i,e1);
 
               
      }
      
    }
      
  // normalize h2 by bin width
  
    Double_t Nbin2 = h2->GetSize();

    for (int i=0;i<Nbin2;i++){

      Double_t x2 = h2->GetBinWidth(i);
      Double_t y2 = h2->GetBinContent(i);
      Double_t z2 = h2->GetBinError(i);
      Double_t r2 = 0.0;
      Double_t e2 = 0.0;
      
      if(x2!=0){

        r2 = y2/x2;    
        e2 = z2/x2;
        h2->SetBinContent(i,r2); 
        h2->SetBinError(i,e2);
 
               
      }
      
    }

      
      
      Double_t msize = 1.2;

      //h1 -> SetMinimum(5*pow(10,-5));
      //h1 -> SetMaximum(1*pow(10,0));

      h1 -> SetMarkerStyle(8);
      h1 -> SetMarkerColor(kBlue);
      h1 -> SetMarkerSize(msize);
      h1 -> SetFillColorAlpha(kBlue,0.7);
      h1 -> SetStats(0);
      

      h2 -> SetMarkerStyle(8);
      h2 -> SetMarkerColor(kRed);
      h2 -> SetMarkerSize(msize);
      h2 -> SetFillColorAlpha(kRed,0.7);
      h2 -> SetStats(0);
    

      h1->Draw("e2");
      h2->Draw("e2 same");
      

      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      legend->AddEntry(h1,"Data","fp");
      legend->AddEntry(h2,"Pythia","fp");
      legend->Draw();
      

     
      h1->GetYaxis()->SetTitle("(N^{jet})^{-1} dN^{jet}/d#eta");
      h1->GetYaxis()->SetTitleSize(0.04);


      h1->SetTitle("");
      
	


///////////////////////////////////////////////////////////////////////////  ratio plot  ////////////////////////////////////////////////////////////////////////////


    c2->cd();
    TPad *pad1 = new TPad("pad1","pad1",0.,0.0,1.,0.3);
    //pad1->SetGridx();
    pad1->SetLeftMargin(0.15);
    pad1->SetBottomMargin(0.3);
    pad1->Draw();
    pad1->cd();
    // pt
    /*
    TH1D *rr = (TH1D*)hp->Clone("rr");
    rr->Divide(hd);
    TH1D *r = (TH1D*)rr->Rebin(15,"",pt_axis);
    */
    // eta
    TH1D *r = (TH1D*)h2->Clone("rr");
    r->Divide(h1);
    
    
    
    r->SetMarkerStyle(8);
    r->SetMarkerColor(kBlack);
    r->SetMarkerSize(1.2);
    r->SetFillColorAlpha(kBlack,0.7);
    r->SetStats(0);
    r->Draw("e2");
    r->GetXaxis()->SetTitle("jet #eta");
    r->GetXaxis()->SetTitleSize(0.09);
    r->GetYaxis()->SetTitle("Pythia / data");
    r->GetYaxis()->SetTitleSize(0.07);
    r->SetTitle("");
    r->SetMinimum(0.8);
    r->SetMaximum(1.2);
    //c1->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/pythiaDataPtRatio.pdf");
    //c1->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/pythiaDataEtaRatio.pdf");


c2->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/pythiaDataEtaCompare.pdf");











}