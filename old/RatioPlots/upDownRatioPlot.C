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

//void jetParticleFraction(TString inf){
void upDownRatioPlot(){
	/// LOAD DATA ///
    TFile* f =TFile::Open("/home/clayton/Analysis/code/jetPtFractionV5.root");
    TH1D* h_jetpt_u;
    TH1D* h_jetpt_d;
    f->GetObject("h_jetpt_u",h_jetpt_u);
    f->GetObject("h_jetpt_d",h_jetpt_d);
    
    TH1D *h1 = (TH1D*)h_jetpt_u->Clone("h1"); // u quarks
    TH1D *h2 = (TH1D*)h_jetpt_d->Clone("h2"); // d quarks
    THStack *hstack = new THStack("hstack","");
   

    // Define the Canvas
    TCanvas *c = new TCanvas("c", "canvas", 800, 800);

    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    h1->SetStats(0);          // No statistics on upper plot
    h1->SetMinimum(0);
    hstack->Add(h2);
    hstack->Add(h1);
    hstack->Draw("hist");        // Draw h2 on top of h1
    auto legend = new TLegend(0.2,0.75,0.3,0.85);
    legend->AddEntry(h1,"u");
    legend->AddEntry(h2,"d");
    legend->Draw();


   	// Do not draw the Y axis label on the upper plot and redraw a small
    // axis instead, in order to avoid the first label (0) to be clipped.
    h1->GetYaxis()->SetLabelSize(0.);
    TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
    axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axis->SetLabelSize(15);
    axis->Draw();

    // lower plot will be in pad
    c->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.1, 1, 0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridx(); // vertical grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad

    // Define the ratio plot
    TH1F *h3 = (TH1F*)h1->Clone("h3");
    h3->Divide(h2);
    h3->SetLineColor(kBlack);
    Double_t rmin = h3->GetMinimum();
    Double_t rmax = h3->GetMaximum();
    h3->SetMinimum(0.9*rmin);  // Define Y ..
    h3->SetMaximum(1.1*rmax); // .. range
    h3->SetStats(0);      // No statistics on lower plot
    h3->SetMarkerStyle(21);
    h3->SetMarkerSize(0.5);
    h3->Draw("e1p");       // Draw the ratio plot

    // h1 settings
    h1->SetFillColor(kGreen-3);
    h1->SetLineWidth(1);

   // Y axis h1 plot settings
   h1->GetYaxis()->SetTitleSize(20);
   h1->GetYaxis()->SetTitleFont(43);
   h1->GetYaxis()->SetTitleOffset(1.55);

   // h2 settings
   h2->SetFillColor(kBlue);
   h2->SetLineWidth(1);

   // Ratio plot (h3) settings
   h3->SetTitle(""); // Remove the ratio title

   // Y axis ratio plot settings
   h3->GetYaxis()->SetTitle("ratio u/d ");
   h3->GetYaxis()->SetNdivisions(505);
   h3->GetYaxis()->SetTitleSize(18);
   h3->GetYaxis()->SetTitleFont(43);
   h3->GetYaxis()->SetTitleOffset(1.55);
   h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h3->GetYaxis()->SetLabelSize(15);
  
   
   // X axis ratio plot settings
   h3->GetXaxis()->SetTitle("Jet pT (GeV/c)");
   h3->GetXaxis()->SetTitleSize(18);
   h3->GetXaxis()->SetTitleFont(43);
   h3->GetXaxis()->SetTitleOffset(3.);
   h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h3->GetXaxis()->SetLabelSize(12);
}