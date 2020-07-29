
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
void ptProjectionPlotPthatRegions(Double_t etacut1=-1.5, Double_t etacut2=1.5){

/// LOAD DATA ///
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_ref.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_reco.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_pthat_50.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_recojets_pthat_30_jetptcut_50_muptcut_2-5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV10_refjets_pthat_30_jetptcut_50_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV10_refjets_pthat_50_jetptcut_50_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/etaPtHistoV10_refjets_pthat_80_jetptcut_50_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_30_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_50_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_refjets_pthat_80_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_recojets_pthat_30_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_recojets_pthat_50_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV1_recojets_pthat_80_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeDataV3_refjets_pthat_50_muptcut_5.root");
//TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V4/rootFiles/makeDataV4_refjets_pthatregions.root");
//  TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_refjets_pthat_30_muptcut_5.root");
TFile *f = TFile::Open("/home/clayton/Analysis/code/makeData/V8/rootFiles/makeDataV8_recojets_pthat_30_muptcut_5.root");

TH2D* h2_1;
TH2D* h2_2;
TH2D* h2_3;


f->GetObject("h2_pthatReg1",h2_1);
f->GetObject("h2_pthatReg2",h2_2);
f->GetObject("h2_pthatReg3",h2_3);

TH1D* hx_test=h2_1->ProjectionX();
TAxis *xaxis = hx_test->GetXaxis();
int firstxbin = xaxis->FindBin(etacut1);
cout << firstxbin << endl;
int lastxbin = xaxis->FindBin(etacut2);
cout << lastxbin << endl;

TH1D* h_jetpt_1 = h2_1->ProjectionY("h_jetpt_1",firstxbin,lastxbin);
TH1D* h_jetpt_2 = h2_2->ProjectionY("h_jetpt_2",firstxbin,lastxbin);
TH1D* h_jetpt_3 = h2_3->ProjectionY("h_jetpt_3",firstxbin,lastxbin);

TH1D *h1 = (TH1D*)h_jetpt_1->Clone("h1"); 
TH1D *h2 = (TH1D*)h_jetpt_2->Clone("h2"); 
TH1D *h3 = (TH1D*)h_jetpt_3->Clone("h3"); 


// generate dN/dp_T plot

    double Nbin = h1->GetSize();

    for (int i=0;i<Nbin;i++){

      double x = h1->GetBinWidth(i);

      double y1 = h1->GetBinContent(i);
      double y2 = h2->GetBinContent(i);
      double y3 = h3->GetBinContent(i);
      
      double z1 = h1->GetBinError(i);
      double z2 = h2->GetBinError(i);
      double z3 = h3->GetBinError(i);

      double r1 = 0;
      double r2 = 0;
      double r3 = 0;
      
      double e1 = 0;
      double e2 = 0;
      double e3 = 0;
     
      if(x!=0){

        r1 = y1/x;
        r2 = y2/x;
        r3 = y3/x;

        e1 = z1/x;
        e2 = z2/x;
        e3 = z3/x;

        h1->SetBinContent(i,r1);
        h2->SetBinContent(i,r2);
        h3->SetBinContent(i,r3);

        h1->SetBinError(i,e1);
        h2->SetBinError(i,e2);
        h3->SetBinError(i,e3);
       
       
      }
      
    }

    // Define the Canvas
      TCanvas *c1 = new TCanvas("c1", "canvas", 1200, 800);
      c1->cd();
      c1->SetLogy();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      //pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1

      // jet pT 
      
            
      double msize = 1.2;

      h1 -> SetMinimum(1*pow(10,-4));
      h1 -> SetMaximum(1*pow(10,4));

      h1 -> SetMarkerStyle(8);
      h1 -> SetMarkerColor(kAzure+9);
      h1 -> SetMarkerSize(msize);
      h1 -> SetFillColorAlpha(kAzure+10,0.5);
      h1 -> SetStats(1); 
      
      h1->Draw("e2");

      auto legend1 = new TLegend(0.7,0.4,0.9,0.5);
      legend1->AddEntry(h1,"inclusive","fp");
      legend1->Draw();
     
      h1->SetTitle("30 < #hat{p}_{T} < 50 GeV/c");
      h1->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      h1->GetYaxis()->SetTitle("#frac{dN^{jet}}{dp_{T}}");
      c1->SaveAs("/home/clayton/Analysis/code/ptProjectionPlotPthatRegions/figures/pthatreg1.png");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Define the Canvas
      TCanvas *c2 = new TCanvas("c2", "canvas", 1200, 800);
      c2->cd();
      c2->SetLogy();
      TPad *pad2 = new TPad("pad2", "pad2", 0., 0., 1, 1);
      //pad1->SetGridx();         // Vertical grid
      pad2->Draw();             // Draw the upper pad: pad1

      // jet pT 
    

      h2 -> SetMinimum(1*pow(10,-5));
      h2 -> SetMaximum(1*pow(10,3));

      h2 -> SetMarkerStyle(22);
      h2 -> SetMarkerColor(kRed);
      h2 -> SetMarkerSize(msize);
      h2 -> SetFillColorAlpha(kRed,0.5);
      h2 -> SetStats(1); 
      
      h2->Draw("e2");

      auto legend2 = new TLegend(0.7,0.4,0.9,0.5);
      legend2->AddEntry(h2,"inclusive","fp");
      legend2->Draw();
     
      h2->SetTitle("50 < #hat{p}_{T} < 80 GeV/c");
      h2->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      h2->GetYaxis()->SetTitle("#frac{dN^{jet}}{dp_{T}}");
      c2->SaveAs("/home/clayton/Analysis/code/ptProjectionPlotPthatRegions/figures/pthatreg2.png");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Define the Canvas
      TCanvas *c3 = new TCanvas("c3", "canvas", 1200, 800);
      c3->cd();
      c3->SetLogy();
      TPad *pad3 = new TPad("pad3", "pad3", 0., 0., 1, 1);
      //pad1->SetGridx();         // Vertical grid
      pad3->Draw();             // Draw the upper pad: pad1

      // jet pT 
    

      h3 -> SetMinimum(1*pow(10,-3));
      h3 -> SetMaximum(1*pow(10,3));

      h3 -> SetMarkerStyle(23);
      h3 -> SetMarkerColor(kGreen-2);
      h3 -> SetMarkerSize(msize);
      h3 -> SetFillColorAlpha(kGreen-2,0.5);
      h3 -> SetStats(1); 
      
      h3->Draw("e2");

      auto legend3 = new TLegend(0.7,0.4,0.9,0.5);
      legend3->AddEntry(h3,"inclusive","fp");
      legend3->Draw();
     
      h3->SetTitle("#hat{p}_{T} > 80 GeV/c");
      h3->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      h3->GetYaxis()->SetTitle("#frac{dN^{jet}}{dp_{T}}");
      c3->SaveAs("/home/clayton/Analysis/code/ptProjectionPlotPthatRegions/figures/pthatreg3.png");





} // END PROGRAM
