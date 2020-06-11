



// V2: new ratio plots recommended by Olga

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
void jetPtFractionPlotV2(){
	/// LOAD DATA ///
    TFile* f =TFile::Open("/home/clayton/Analysis/code/jetPtFractionV5.root");
    TH1D* h_jetpt;
    TH1D* h_jetpt_g;
    TH1D* h_jetpt_d;
    TH1D* h_jetpt_dbar;
    TH1D* h_jetpt_dall;
    TH1D* h_jetpt_u;
    TH1D* h_jetpt_ubar;
    TH1D* h_jetpt_uall;
    TH1D* h_jetpt_s;
    TH1D* h_jetpt_sbar;
    TH1D* h_jetpt_sall;
    TH1D* h_jetpt_c;
    TH1D* h_jetpt_cbar;
    TH1D* h_jetpt_call;
    TH1D* h_jetpt_b;
    TH1D* h_jetpt_bbar;
    TH1D* h_jetpt_ball;
    TH1D* h_jetpt_hq;
    TH1D* h_jetpt_lq;
    TH1D* h_jetpt_qall;
    TH1D* h_jetpt_other;
    TH1D* h_jetpt_ud;
    TH1D* h_jetpt_ubardbar;
    f->GetObject("h_jetpt",h_jetpt);
    f->GetObject("h_jetpt_g",h_jetpt_g);
    f->GetObject("h_jetpt_d",h_jetpt_d);
    f->GetObject("h_jetpt_dbar",h_jetpt_dbar);
    f->GetObject("h_jetpt_dall",h_jetpt_dall);
    f->GetObject("h_jetpt_u",h_jetpt_u);
    f->GetObject("h_jetpt_ubar",h_jetpt_ubar);
    f->GetObject("h_jetpt_uall",h_jetpt_uall);
    f->GetObject("h_jetpt_s",h_jetpt_s);
    f->GetObject("h_jetpt_sbar",h_jetpt_sbar);
    f->GetObject("h_jetpt_sall",h_jetpt_sall);
    f->GetObject("h_jetpt_c",h_jetpt_c);
    f->GetObject("h_jetpt_cbar",h_jetpt_cbar);
    f->GetObject("h_jetpt_call",h_jetpt_call);
    f->GetObject("h_jetpt_b",h_jetpt_b);
    f->GetObject("h_jetpt_bbar",h_jetpt_bbar);
    f->GetObject("h_jetpt_ball",h_jetpt_ball);
    f->GetObject("h_jetpt_hq",h_jetpt_hq);
    f->GetObject("h_jetpt_lq",h_jetpt_lq);
    f->GetObject("h_jetpt_qall",h_jetpt_qall);
    f->GetObject("h_jetpt_other",h_jetpt_other);
    f->GetObject("h_jetpt_ud",h_jetpt_ud);
    f->GetObject("h_jetpt_ubardbar",h_jetpt_ubardbar);




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    // first plot: heavy flavors only with quarks and antiquarks grouped together
    
    // Define the Canvas
      TCanvas *c1 = new TCanvas("c1", "canvas", 800, 800);
      c1->cd();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
       THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
      TH1D *h1 = (TH1D*)h_jetpt_sall->Clone("h1"); // s quarks/antiquarks
      TH1D *h2 = (TH1D*)h_jetpt_call->Clone("h2"); // c quarks/antiquarks
      TH1D *h3 = (TH1D*)h_jetpt_bbar->Clone("h3"); // b quarks/antiquarks   
      h1 -> SetFillColor(kCyan+1);
      h2 -> SetFillColor(kAzure-3);
      h3 -> SetFillColor(kTeal+5);
      h1 -> SetStats(0);
      h1 -> Divide(h_jetpt);
      h2 -> Divide(h_jetpt);
      h3 -> Divide(h_jetpt);      
      fraction_stack->Add(h1);
      fraction_stack->Add(h2);
      fraction_stack->Add(h3);
      fraction_stack->Draw("ehist");
      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      legend->AddEntry(h3,"b & #bar{b}","f");
      legend->AddEntry(h2,"c & #bar{c}","f");
      legend->AddEntry(h1,"s & #bar{s}","f");
      legend->Draw();
      fraction_stack->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      //fraction_stack->GetXaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);
      c1->SaveAs("jetPtFractionPlot1.pdf");
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    // second plot: light flavors only with quarks and antiquarks combined
    
    // Define the Canvas
      TCanvas *c2 = new TCanvas("c2", "canvas", 800, 800);
      c2->cd();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
      TH1D *h1 = (TH1D*)h_jetpt_dall->Clone("h1"); // s quarks/antiquarks
      TH1D *h2 = (TH1D*)h_jetpt_uall->Clone("h2"); // c quarks/antiquarks
   
      h1 -> SetFillColor(kRed-8);
      h2 -> SetFillColor(kGreen-7);
      h1 -> SetStats(0);
      h1 -> Divide(h_jetpt);
      h2 -> Divide(h_jetpt);   
      fraction_stack->Add(h1);
      fraction_stack->Add(h2);
      fraction_stack->Draw("ehist");
      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      legend->AddEntry(h2,"d & #bar{d}","f");
      legend->AddEntry(h1,"u & #bar{u}","f");
      legend->Draw();
      fraction_stack->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      //fraction_stack->GetXaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);
      c2->SaveAs("jetPtFractionPlot2.pdf");

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    // third plot: gluons and all light quarks
    
    // Define the Canvas
      TCanvas *c3 = new TCanvas("c3", "canvas", 800, 800);
      c3->cd();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
      TH1D *h1 = (TH1D*)h_jetpt_lq->Clone("h1"); // all quarks
      TH1D *h2 = (TH1D*)h_jetpt_g->Clone("h2"); // gluons
   
      h1 -> SetFillColor(kOrange+1);
      h2 -> SetFillColor(kCyan+1);
      h1 -> SetStats(0);
      h1 -> Divide(h_jetpt);
      h2 -> Divide(h_jetpt);   
      fraction_stack->Add(h1);
      fraction_stack->Add(h2);
      fraction_stack->Draw("ehist");
      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      legend->AddEntry(h2,"light quarks","f");
      legend->AddEntry(h1,"gluons","f");
      legend->Draw();
      fraction_stack->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      //fraction_stack->GetXaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);
      c3->SaveAs("jetPtFractionPlot3.pdf");

*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // fourth plot: gluons, u & d together, ubar & dbar together, heavy's together
    
    // Define the Canvas
      TCanvas *c4 = new TCanvas("c3", "canvas", 800, 800);
      c4->cd();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
      TH1D *h1 = (TH1D*)h_jetpt_ud->Clone("h1"); // all quarks
      TH1D *h2 = (TH1D*)h_jetpt_ubardbar->Clone("h2"); // gluons
      TH1D *h3 = (TH1D*)h_jetpt_hq->Clone("h3"); // gluons
      TH1D *h4 = (TH1D*)h_jetpt_g->Clone("h4"); // gluons
   
      h1 -> SetFillColor(kOrange+4);
      h2 -> SetFillColor(kCyan+6);
      h3 -> SetFillColor(kBlue);
      h4 -> SetFillColor(kRed);
      h1 -> SetStats(0);
      h1 -> Divide(h_jetpt);
      h2 -> Divide(h_jetpt);  
      h3 -> Divide(h_jetpt);
      h4 -> Divide(h_jetpt); 
      fraction_stack->Add(h1);
      fraction_stack->Add(h2);
      fraction_stack->Add(h3);
      fraction_stack->Add(h4);
      fraction_stack->Draw("ehist");
      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      legend->AddEntry(h4,"gluons","f");
      legend->AddEntry(h3,"heavy quarks","f");
      legend->AddEntry(h2,"#bar{u} & #bar{d}","f");
      legend->AddEntry(h1,"u & d","f");
      legend->Draw();
      fraction_stack->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      //fraction_stack->GetXaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);
      c4->SaveAs("jetPtFractionPlot3.pdf");







}