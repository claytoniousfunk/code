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
void jetPtFractionPlot(){
	/// LOAD DATA ///
    TFile* f =TFile::Open("/home/clayton/Analysis/code/jetPtFractionV5.root");
    TH1D* h_jetpt;
    TH1D* h_jetpt_g;
    TH1D* h_jetpt_d;
    TH1D* h_jetpt_dbar;
    TH1D* h_jetpt_u;
    TH1D* h_jetpt_ubar;
    TH1D* h_jetpt_hq;
    TH1D* h_jetpt_other;
    f->GetObject("h_jetpt",h_jetpt);
    f->GetObject("h_jetpt_g",h_jetpt_g);
    f->GetObject("h_jetpt_d",h_jetpt_d);
    f->GetObject("h_jetpt_dbar",h_jetpt_dbar);
    f->GetObject("h_jetpt_u",h_jetpt_u);
    f->GetObject("h_jetpt_ubar",h_jetpt_ubar);
    f->GetObject("h_jetpt_hq",h_jetpt_hq);
    f->GetObject("h_jetpt_other",h_jetpt_other);
    
    // Define the Canvas
        TCanvas *c1 = new TCanvas("c1", "canvas", 800, 800);
        c1->cd();
        TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
        pad1->SetGridx();         // Vertical grid
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad
        //pad1->SetLogy();
       THStack *fraction_stack1 = new THStack("fraction_stack1","");
      // jet pT stack
       TH1D *h1 = (TH1D*)h_jetpt_g->Clone("h1"); // gluons
       TH1D *h2 = (TH1D*)h_jetpt_d->Clone("h2"); // d quarks
       TH1D *h3 = (TH1D*)h_jetpt_dbar->Clone("h3"); // d antiquarks
       TH1D *h4 = (TH1D*)h_jetpt_u->Clone("h4"); // u quarks
       TH1D *h5 = (TH1D*)h_jetpt_ubar->Clone("h5"); // u antiquarks
       TH1D *h6 = (TH1D*)h_jetpt_hq->Clone("h6"); // heavy quarks
       TH1D *h7 = (TH1D*)h_jetpt_other->Clone("h7"); // unidentified particles

       h1 -> SetFillColor(kCyan+1);
       h2 -> SetFillColor(kBlue);
       h3 -> SetFillColor(kOrange+8);
       h4 -> SetFillColor(kGreen-3);
       h5 -> SetFillColor(kYellow-7);
       h6 -> SetFillColor(kRed);
     h7 -> SetFillColor(kBlue-10);
       h1 -> SetStats(0);
       h1 -> Divide(h_jetpt);
       h2 -> Divide(h_jetpt);
       h3 -> Divide(h_jetpt);
       h4 -> Divide(h_jetpt);
       h5 -> Divide(h_jetpt);
       h6 -> Divide(h_jetpt);
     h7 -> Divide(h_jetpt);
  
       fraction_stack1->Add(h1);
       fraction_stack1->Add(h2);
       fraction_stack1->Add(h3);
       fraction_stack1->Add(h4);
       fraction_stack1->Add(h5);
       fraction_stack1->Add(h6);
     //fraction_stack1->Add(h7);
       fraction_stack1->Draw("ehist");
       auto legend = new TLegend(0.91,0.4,0.99,0.6);
       TLatex latex;
       legend->AddEntry(h6,"heavy","f");
       legend->AddEntry(h5,"#bar{u}","f");
       legend->AddEntry(h4,"u","f");
       legend->AddEntry(h3,"#bar{d}","f");
       legend->AddEntry(h2,"d","f");
       legend->AddEntry(h1,"g","f");
       
       
       
       
     //legend->AddEntry(h7,"unidentified");
       legend->Draw();

       fraction_stack1->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
       //fraction_stack->GetXaxis()->SetTitleSize(18);
       //fraction_stack->GetXaxis()->SetTitleFont(43);
       //fraction_stack->GetXaxis()->SetTitleOffset(3.);
       //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       //fraction_stack->GetXaxis()->SetLabelSize(12);

       c1->SaveAs("jetPtFractionPlot1.pdf");


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


       // Define the Canvas
        TCanvas *c2 = new TCanvas("c2", "canvas", 800, 800);
        c1->cd();
        TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
        pad1->SetGridx();         // Vertical grid
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad
        //pad1->SetLogy();
       THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
       TH1D *h1 = (TH1D*)h_jetpt_g->Clone("h1"); // gluons
       TH1D *h2 = (TH1D*)h_jetpt_d->Clone("h2"); // d quarks
       TH1D *h3 = (TH1D*)h_jetpt_dbar->Clone("h3"); // d antiquarks
       TH1D *h4 = (TH1D*)h_jetpt_u->Clone("h4"); // u quarks
       TH1D *h5 = (TH1D*)h_jetpt_ubar->Clone("h5"); // u antiquarks
       TH1D *h6 = (TH1D*)h_jetpt_hq->Clone("h6"); // heavy quarks
       TH1D *h7 = (TH1D*)h_jetpt_other->Clone("h7"); // unidentified particles

       h1 -> SetFillColor(kCyan+1);
       h2 -> SetFillColor(kBlue);
       h3 -> SetFillColor(kOrange+8);
       h4 -> SetFillColor(kGreen-3);
       h5 -> SetFillColor(kYellow-7);
       h6 -> SetFillColor(kRed);
     h7 -> SetFillColor(kBlue-10);
       h1 -> SetStats(0);
       h1 -> Divide(h_jetpt);
       h2 -> Divide(h_jetpt);
       h3 -> Divide(h_jetpt);
       h4 -> Divide(h_jetpt);
       h5 -> Divide(h_jetpt);
       h6 -> Divide(h_jetpt);
     h7 -> Divide(h_jetpt);
  
       fraction_stack->Add(h1);
       fraction_stack->Add(h2);
       fraction_stack->Add(h3);
       fraction_stack->Add(h4);
       fraction_stack->Add(h5);
       fraction_stack->Add(h6);
     //fraction_stack->Add(h7);
       fraction_stack->Draw("ehist");
       auto legend = new TLegend(0.91,0.4,0.99,0.6);
       TLatex latex;
       legend->AddEntry(h6,"heavy","f");
       legend->AddEntry(h5,"#bar{u}","f");
       legend->AddEntry(h4,"u","f");
       legend->AddEntry(h3,"#bar{d}","f");
       legend->AddEntry(h2,"d","f");
       legend->AddEntry(h1,"g","f");
       
       
       
       
     //legend->AddEntry(h7,"unidentified");
       legend->Draw();

       fraction_stack->GetXaxis()->SetTitle("Jet pT (GeV/c)");
       //fraction_stack->GetXaxis()->SetTitleSize(18);
       //fraction_stack->GetXaxis()->SetTitleFont(43);
       //fraction_stack->GetXaxis()->SetTitleOffset(3.);
       //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       //fraction_stack->GetXaxis()->SetLabelSize(12);

       c1->SaveAs("jetPtFractionPlot.pdf");


}