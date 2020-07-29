





#include "myProcesses/hiforest/plugin/eventMap_hiForest.h"
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
void etaProjectionPlotV1(Double_t ptcut1=100, Double_t ptcut2=600){

/// LOAD DATA ///
TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_10.root");

TH2D* h2;
TH2D* h2_d;
TH2D* h2_u;
TH2D* h2_s;
TH2D* h2_c;
TH2D* h2_b;
TH2D* h2_dbar;
TH2D* h2_ubar;
TH2D* h2_sbar;
TH2D* h2_cbar;
TH2D* h2_bbar;
TH2D* h2_hq;
TH2D* h2_lq;
TH2D* h2_qall;
TH2D* h2_dall;
TH2D* h2_uall;
TH2D* h2_sall;
TH2D* h2_call;
TH2D* h2_ball;
TH2D* h2_g;
TH2D* h2_ud;
TH2D* h2_ubardbar;
TH2D* h2_other;

f->GetObject("h2",h2);
f->GetObject("h2_d",h2_d);
f->GetObject("h2_u",h2_u);
f->GetObject("h2_s",h2_s);
f->GetObject("h2_c",h2_c);
f->GetObject("h2_b",h2_b);
f->GetObject("h2_dbar",h2_dbar);
f->GetObject("h2_ubar",h2_ubar);
f->GetObject("h2_sbar",h2_sbar);
f->GetObject("h2_cbar",h2_cbar);
f->GetObject("h2_bbar",h2_bbar);
f->GetObject("h2_hq",h2_hq);
f->GetObject("h2_lq",h2_lq);
f->GetObject("h2_qall",h2_qall);
f->GetObject("h2_dall",h2_dall);
f->GetObject("h2_uall",h2_uall);
f->GetObject("h2_sall",h2_sall);
f->GetObject("h2_call",h2_call);
f->GetObject("h2_ball",h2_ball);
f->GetObject("h2_g",h2_g);
f->GetObject("h2_ud",h2_ud);
f->GetObject("h2_ubardbar",h2_ubardbar);
f->GetObject("h2_other",h2_other);

TH1D* hy_test=h2->ProjectionY();
TAxis *yaxis = hy_test->GetXaxis();
Int_t firstybin = yaxis->FindBin(ptcut1);
cout << firstybin << endl;
Int_t lastybin = yaxis->FindBin(ptcut2);
cout << lastybin << endl;

TH1D* h_jeteta = h2->ProjectionX("All jets",firstybin,lastybin);
TH1D* h_jeteta_g = h2_g->ProjectionX("gluon jets",firstybin,lastybin);
TH1D* h_jeteta_lq = h2_lq->ProjectionX("light quark jets",firstybin,lastybin);
TH1D* h_jeteta_sall = h2_sall->ProjectionX("s & #bar{s}",firstybin,lastybin);
TH1D* h_jeteta_call = h2_call->ProjectionX("c & #bar{c}",firstybin,lastybin);
TH1D* h_jeteta_ball = h2_ball->ProjectionX("b & #bar{b}",firstybin,lastybin);


/*
// Define the Canvas

      TCanvas *c1 = new TCanvas("c", "canvas", 800, 800);
      c1->cd();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
      TH1D *h0 = (TH1D*)h_jeteta->Clone("h0"); // all jets
      TH1D *h1 = (TH1D*)h_jeteta_g->Clone("h1"); // gluons
      TH1D *h2a = (TH1D*)h_jeteta_lq->Clone("h2a"); // light quarks
      TH1D *h3 = (TH1D*)h_jeteta_ball->Clone("h3"); // b and bbar
      TH1D *h4 = (TH1D*)h_jeteta_sall->Clone("h4"); // s and sbar
      TH1D *h5 = (TH1D*)h_jeteta_call->Clone("h5"); // c and cbar

      Double_t msize = 1.2;


    

      h1 -> SetMarkerStyle(8);
      h1 -> SetMarkerColor(kAzure+9);
      h1 -> SetMarkerSize(msize);
      h1 -> SetFillColorAlpha(kAzure+10,0.5);
      h1 -> SetStats(0);
      h1 -> Divide(h0);
      
      h2a -> SetMarkerStyle(21);
      h2a -> SetMarkerColor(kTeal+9);
      h2a -> SetMarkerSize(msize);
      h2a -> SetFillColorAlpha(kTeal+10,0.5);
      h2a -> SetStats(0);
      h2a -> Divide(h0);
      
      h3 -> SetMarkerStyle(22);
      h3 -> SetMarkerColor(kRed-4);
      h3 -> SetMarkerSize(msize);
      h3 -> SetFillColorAlpha(kRed-7,0.5);
      h3 -> SetStats(0);
      h3 -> Divide(h0);
     
      h4 -> SetMarkerStyle(33);
      h4 -> SetMarkerColor(kMagenta+3);
      h4 -> SetMarkerSize(msize);
      h4 -> SetFillColorAlpha(kMagenta+1,0.5);
      h4 -> SetStats(0);
      h4 -> Divide(h0);
      
      h5 -> SetMarkerStyle(34);
      h5 -> SetMarkerColor(kGreen-5);
      h5 -> SetMarkerSize(msize);
      h5 -> SetFillColorAlpha(kGreen-8,0.5);
      h5 -> SetStats(0);
      h5 -> Divide(h0);
     // h5 -> GetXaxis()->SetLimits(-1.5,1.5);
     
     
      fraction_stack->Add(h5);
      fraction_stack->Add(h4);
      fraction_stack->Add(h3);
      fraction_stack->Add(h2a);
      fraction_stack->Add(h1);


      fraction_stack->Draw("e1 hist");
      auto legend = new TLegend(0.91,0.4,0.99,0.6);

      
      legend->AddEntry(h1,"gluons","fp");
      legend->AddEntry(h2a,"light quarks","fp");
      legend->AddEntry(h3,"b & #bar{b}","fp");
      legend->AddEntry(h4,"s & #bar{s}","fp");
      legend->AddEntry(h5,"c & #bar{c}","fp");
      legend->Draw();
      fraction_stack->GetXaxis()->SetTitle("Jet #eta");
      fraction_stack->GetYaxis()->SetTitle("Jet particle fraction");
      //fraction_stack->GetXaxis()->SetLimits(-1.5,1.5);
      //fraction_stack->GetXaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);
      c1->SaveAs("etaProjectionPlotV1-1.pdf");

     
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    // fifth plot: eta profile of gluons, light quarks (together), and heavy quarks (separate)
    
    // Define the Canvas
      TCanvas *c2 = new TCanvas("c2", "canvas", 1100, 800);
      c2->cd();
      c2->SetLogy();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1

      // jet pT 
      TH1D *h1 = (TH1D*)h_jeteta_g->Clone("h1"); // gluons
      TH1D *h2a = (TH1D*)h_jeteta_lq->Clone("h2a"); // light quarks
      TH1D *h3 = (TH1D*)h_jeteta_ball->Clone("h3"); // b and bbar
      TH1D *h4 = (TH1D*)h_jeteta_sall->Clone("h4"); // s and sbar
      TH1D *h5 = (TH1D*)h_jeteta_call->Clone("h5"); // c and cbar
      
      Double_t msize = 1.2;

      h5 -> SetMinimum(1*pow(10,-2));
      h5 -> SetMaximum(1*pow(10,2));

      h1 -> SetMarkerStyle(8);
      h1 -> SetMarkerColor(kAzure+9);
      h1 -> SetMarkerSize(msize);
      h1 -> SetFillColorAlpha(kAzure+10,0.5);
      h1 -> SetStats(0);
      

      h2a -> SetMarkerStyle(21);
      h2a -> SetMarkerColor(kTeal+9);
      h2a -> SetMarkerSize(msize);
      h2a -> SetFillColorAlpha(kTeal+10,0.5);
      h2a -> SetStats(0);
    

      h3 -> SetMarkerStyle(22);
      h3 -> SetMarkerColor(kRed-4);
      h3 -> SetMarkerSize(msize);
      h3 -> SetFillColorAlpha(kRed-7,0.5);
      h3 -> SetStats(0);
      

      h4 -> SetMarkerStyle(33);
      h4 -> SetMarkerColor(kMagenta+3);
      h4 -> SetMarkerSize(msize);
      h4 -> SetFillColorAlpha(kMagenta+1,0.5);
      h4 -> SetStats(0);
      

      h5 -> SetMarkerStyle(34);
      h5 -> SetMarkerColor(kGreen-5);
      h5 -> SetMarkerSize(msize);
      h5 -> SetFillColorAlpha(kGreen-8,0.5);
      h5 -> SetStats(0);
      
      h5->Draw("e2");
      h4->Draw("same e2");
      h3->Draw("same e2");
      h2a->Draw("same e2");
      h1->Draw("same e2");

      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      legend->AddEntry(h1,"gluons","fp");
      legend->AddEntry(h2a,"light quarks","fp");
      legend->AddEntry(h3,"b & #bar{b}","fp");
      legend->AddEntry(h4,"s & #bar{s}","fp");
      legend->AddEntry(h5,"c & #bar{c}","fp");
      legend->Draw();
      
      //fraction_stack->GetYaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);
      h5->GetXaxis()->SetTitle("Jet #eta");
      h5->GetYaxis()->SetTitle("#frac{dN^{jet}}{d#eta}");
      h5->SetTitle("");

      c2->SaveAs("etaProjectionPlotV1-2.pdf");









} // END PROGRAM