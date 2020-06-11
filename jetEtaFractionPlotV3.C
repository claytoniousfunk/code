// 
// V3: based off of JetPtFractionPlotV3.C

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
void jetEtaFractionPlotV3(){
	/// LOAD DATA ///
    TFile* f =TFile::Open("/home/clayton/Analysis/code/jetEtaFractionV6.root");
    TH1D* h_jeteta;
    TH1D* h_jeteta_g;
    TH1D* h_jeteta_d;
    TH1D* h_jeteta_dbar;
    TH1D* h_jeteta_dall;
    TH1D* h_jeteta_u;
    TH1D* h_jeteta_ubar;
    TH1D* h_jeteta_uall;
    TH1D* h_jeteta_s;
    TH1D* h_jeteta_sbar;
    TH1D* h_jeteta_sall;
    TH1D* h_jeteta_c;
    TH1D* h_jeteta_cbar;
    TH1D* h_jeteta_call;
    TH1D* h_jeteta_b;
    TH1D* h_jeteta_bbar;
    TH1D* h_jeteta_ball;
    TH1D* h_jeteta_hq;
    TH1D* h_jeteta_lq;
    TH1D* h_jeteta_qall;
    TH1D* h_jeteta_other;
    TH1D* h_jeteta_ud;
    TH1D* h_jeteta_ubardbar;
    f->GetObject("h_jeteta",h_jeteta);
    f->GetObject("h_jeteta_g",h_jeteta_g);
    f->GetObject("h_jeteta_d",h_jeteta_d);
    f->GetObject("h_jeteta_dbar",h_jeteta_dbar);
    f->GetObject("h_jeteta_dall",h_jeteta_dall);
    f->GetObject("h_jeteta_u",h_jeteta_u);
    f->GetObject("h_jeteta_ubar",h_jeteta_ubar);
    f->GetObject("h_jeteta_uall",h_jeteta_uall);
    f->GetObject("h_jeteta_s",h_jeteta_s);
    f->GetObject("h_jeteta_sbar",h_jeteta_sbar);
    f->GetObject("h_jeteta_sall",h_jeteta_sall);
    f->GetObject("h_jeteta_c",h_jeteta_c);
    f->GetObject("h_jeteta_cbar",h_jeteta_cbar);
    f->GetObject("h_jeteta_call",h_jeteta_call);
    f->GetObject("h_jeteta_b",h_jeteta_b);
    f->GetObject("h_jeteta_bbar",h_jeteta_bbar);
    f->GetObject("h_jeteta_ball",h_jeteta_ball);
    f->GetObject("h_jeteta_hq",h_jeteta_hq);
    f->GetObject("h_jeteta_lq",h_jeteta_lq);
    f->GetObject("h_jeteta_qall",h_jeteta_qall);
    f->GetObject("h_jeteta_other",h_jeteta_other);
    f->GetObject("h_jeteta_ud",h_jeteta_ud);
    f->GetObject("h_jeteta_ubardbar",h_jeteta_ubardbar);


    







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
      c1->SaveAs("jetEtaFractionPlotV3-1.pdf");
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
      c2->SaveAs("jetEtaFractionPlotV3-2.pdf");

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // third plot: gluons, all light quarks (together), and all heavys (separate)
    
    // Define the Canvas
      TCanvas *c3 = new TCanvas("c3", "canvas", 800, 800);
      c3->cd();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
      TH1D *h0 = (TH1D*)h_jeteta->Clone("h0"); // all jets
      TH1D *h1 = (TH1D*)h_jeteta_g->Clone("h1"); // gluons
      TH1D *h2 = (TH1D*)h_jeteta_lq->Clone("h2"); // light quarks
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
      
      h2 -> SetMarkerStyle(21);
      h2 -> SetMarkerColor(kTeal+9);
      h2 -> SetMarkerSize(msize);
      h2 -> SetFillColorAlpha(kTeal+10,0.5);
      h2 -> SetStats(0);
      h2 -> Divide(h0);
      
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
      fraction_stack->Add(h2);
      fraction_stack->Add(h1);


      fraction_stack->Draw("e1 hist");
      auto legend = new TLegend(0.91,0.4,0.99,0.6);

      
      legend->AddEntry(h1,"gluons","fp");
      legend->AddEntry(h2,"light quarks","fp");
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
      c3->SaveAs("jetEtaFractionPlotV3-3.pdf");


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    // fourth plot: gluons, u & d together, ubar & dbar together, heavy's together
    
    // Define the Canvas
      TCanvas *c4 = new TCanvas("c3", "canvas", 800, 800);
      c4->cd();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
      TH1D *h1 = (TH1D*)h_jetpt_g->Clone("h1"); // gluons
      TH1D *h2 = (TH1D*)h_jetpt_ud->Clone("h2"); // u and d
      TH1D *h3 = (TH1D*)h_jetpt_ubardbar->Clone("h3"); // ubar and dbar
      TH1D *h4 = (TH1D*)h_jetpt_hq->Clone("h4"); // heavy quarks
      
   
      h1 -> SetFillColor(kMagenta-10);
      h2 -> SetFillColor(kAzure+9);
      h3 -> SetFillColor(kGreen-9);
      h4 -> SetFillColor(kGreen+1);
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
      legend->AddEntry(h4,"heavy quarks","f");
      legend->AddEntry(h3,"#bar{u} & #bar{d}","f");
      legend->AddEntry(h2,"u & d","f");
      legend->AddEntry(h1,"gluons","f");
      legend->Draw();
      
      //fraction_stack->GetYaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);
      fraction_stack->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      fraction_stack->GetYaxis()->SetTitle("#frac{dN^{jet}}{dp_{T}}");

      c4->SaveAs("jetEtaFractionPlotV3-4.pdf");

*/



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    // fifth plot: eta profile of gluons, light quarks (together), and heavy quarks (separate)
    
    // Define the Canvas
      TCanvas *c5 = new TCanvas("c5", "canvas", 1100, 800);
      c5->cd();
      c5->SetLogy();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1

      // jet pT 
      TH1D *h1 = (TH1D*)h_jeteta_g->Clone("h1"); // gluons
      TH1D *h2 = (TH1D*)h_jeteta_lq->Clone("h2"); // light quarks
      TH1D *h3 = (TH1D*)h_jeteta_ball->Clone("h3"); // b and bbar
      TH1D *h4 = (TH1D*)h_jeteta_sall->Clone("h4"); // s and sbar
      TH1D *h5 = (TH1D*)h_jeteta_call->Clone("h5"); // c and cbar
      
      Double_t msize = 1.2;

      h5 -> SetMinimum(2*pow(10,4));
      h5 -> SetMaximum(3*pow(10,6));

      h1 -> SetMarkerStyle(8);
      h1 -> SetMarkerColor(kAzure+9);
      h1 -> SetMarkerSize(msize);
      h1 -> SetFillColorAlpha(kAzure+10,0.5);
      h1 -> SetStats(0);
      

      h2 -> SetMarkerStyle(21);
      h2 -> SetMarkerColor(kTeal+9);
      h2 -> SetMarkerSize(msize);
      h2 -> SetFillColorAlpha(kTeal+10,0.5);
      h2 -> SetStats(0);
    

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
      h2->Draw("same e2");
      h1->Draw("same e2");

      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      legend->AddEntry(h1,"gluons","fp");
      legend->AddEntry(h2,"light quarks","fp");
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

      c5->SaveAs("jetEtaFractionPlotV3-5.pdf");
*/








}