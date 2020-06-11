

//



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
void muRelPtProjectionPlotV1(Float_t ptcut1=50.0, Float_t ptcut2=600.0){

/// LOAD DATA ///
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_5.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_12-5.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_16.root");
TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_20.root");
TH2D* h2_relPtJetPt;
TH2D* h2_relPtJetPt_g;
TH2D* h2_relPtJetPt_lq;
TH2D* h2_relPtJetPt_sall;
TH2D* h2_relPtJetPt_ball;
TH2D* h2_relPtJetPt_call;


f->GetObject("h2_relPtJetPt",h2_relPtJetPt);
f->GetObject("h2_relPtJetPt_g",h2_relPtJetPt_g);
f->GetObject("h2_relPtJetPt_lq",h2_relPtJetPt_lq);
f->GetObject("h2_relPtJetPt_sall",h2_relPtJetPt_sall);
f->GetObject("h2_relPtJetPt_ball",h2_relPtJetPt_ball);
f->GetObject("h2_relPtJetPt_call",h2_relPtJetPt_call);

TH1D* hy_test=h2_relPtJetPt->ProjectionY();
TAxis *yaxis = hy_test->GetXaxis();
Int_t firstybin = yaxis->FindBin(ptcut1);
cout << firstybin << endl;
Int_t lastybin = yaxis->FindBin(ptcut2);
cout << lastybin << endl;

TH1D* h0 = h2_relPtJetPt->ProjectionX("All jets",firstybin,lastybin);
TH1D* h1 = h2_relPtJetPt_g->ProjectionX("gluon jets",firstybin,lastybin);
TH1D* h2a = h2_relPtJetPt_lq->ProjectionX("light quark jets",firstybin,lastybin);
TH1D* h3 = h2_relPtJetPt_ball->ProjectionX("s & #bar{s}",firstybin,lastybin);
TH1D* h4 = h2_relPtJetPt_sall->ProjectionX("c & #bar{c}",firstybin,lastybin);
TH1D* h5 = h2_relPtJetPt_call->ProjectionX("b & #bar{b}",firstybin,lastybin);



// Define the Canvas

      TCanvas *c1 = new TCanvas("c", "canvas", 1000, 800);
      c1->cd();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
      /*
      TH1D *h0 = (TH1D*)h_jeteta->Clone("h0"); // all jets
      TH1D *h1 = (TH1D*)h_jeteta_g->Clone("h1"); // gluons
      TH1D *h2a = (TH1D*)h_jeteta_lq->Clone("h2a"); // light quarks
      TH1D *h3 = (TH1D*)h_jeteta_ball->Clone("h3"); // b and bbar
      TH1D *h4 = (TH1D*)h_jeteta_sall->Clone("h4"); // s and sbar
      TH1D *h5 = (TH1D*)h_jeteta_call->Clone("h5"); // c and cbar
      */

      Double_t msize = 1.2;


    

      h1 -> SetMarkerStyle(8);
      h1 -> SetMarkerColor(kAzure+9);
      h1 -> SetMarkerSize(msize);
      h1 -> SetFillColorAlpha(kAzure+10,0.5);
      h1 -> SetStats(0);
      //h1 -> Divide(h0);
      
      h2a -> SetMarkerStyle(21);
      h2a -> SetMarkerColor(kTeal+9);
      h2a -> SetMarkerSize(msize);
      h2a -> SetFillColorAlpha(kTeal+10,0.5);
      h2a -> SetStats(0);
      //h2a -> Divide(h0);
      
      h3 -> SetMarkerStyle(22);
      h3 -> SetMarkerColor(kRed-4);
      h3 -> SetMarkerSize(msize);
      h3 -> SetFillColorAlpha(kRed-7,0.5);
      h3 -> SetStats(0);
      //h3 -> Divide(h0);
     
      h4 -> SetMarkerStyle(33);
      h4 -> SetMarkerColor(kMagenta+3);
      h4 -> SetMarkerSize(msize);
      h4 -> SetFillColorAlpha(kMagenta+1,0.5);
      h4 -> SetStats(0);
      //h4 -> Divide(h0);
      
      h5 -> SetMarkerStyle(34);
      h5 -> SetMarkerColor(kGreen-5);
      h5 -> SetMarkerSize(msize);
      h5 -> SetFillColorAlpha(kGreen-8,0.5);
      h5 -> SetStats(0);
      //h5 -> Divide(h0);
     // h5 -> GetXaxis()->SetLimits(-1.5,1.5);
     
     
      fraction_stack->Add(h1);
      fraction_stack->Add(h2a);
      fraction_stack->Add(h4);
      fraction_stack->Add(h5);
      fraction_stack->Add(h3);

      fraction_stack->Draw("e1 hist");
      //h1->Draw("e2");
      //h2a->Draw("e2 same");
      //h3->Draw("e2 same");
      //h4->Draw("e2 same");
      //h5->Draw("e2 same");



      auto legend = new TLegend(0.905,0.4,0.995,0.7);

      
      legend->AddEntry(h3,"b","fp");
      legend->AddEntry(h5,"c","fp");
      legend->AddEntry(h4,"s","fp");
      legend->AddEntry(h2a,"u,d","fp");
      legend->AddEntry(h1,"g","fp");
      legend->Draw();

      //h1->GetXaxis()->SetTitle("Muon rel-p_{T} (GeV/c)");
      //h1->GetYaxis()->SetRangeUser(0.,10000.);

      //h1->SetTitle("p_{T,cut}^{#mu} = 10 GeV/c, 50 < p_{T}^{jet} < 100 GeV/c");
      //h1->SetTitle("p_{T,cut}^{#mu} = 10 GeV/c, 100 < p_{T}^{jet} < 150 GeV/c");
      //h1->SetTitle("p_{T,cut}^{#mu} = 10 GeV/c, 150 < p_{T}^{jet} < 200 GeV/c");
      //h1->SetTitle("p_{T,cut}^{#mu} = 10 GeV/c, 200 < p_{T}^{jet} < 250 GeV/c");
      //h1->SetTitle("p_{T,cut}^{#mu} = 10 GeV/c, 250 < p_{T}^{jet} < 300 GeV/c");
      //h1->SetTitle("p_{T,cut}^{#mu} = 10 GeV/c, 50 < p_{T}^{jet} < 600 GeV/c");

      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 5 GeV/c, 50 < p_{T}^{jet} < 600 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 12.5 GeV/c, 50 < p_{T}^{jet} < 600 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 16 GeV/c, 50 < p_{T}^{jet} < 600 GeV/c");
      fraction_stack->SetTitle("p_{T,cut}^{#mu} = 20 GeV/c, 50 < p_{T}^{jet} < 600 GeV/c");

      fraction_stack->GetXaxis()->SetTitle("p_{T,rel}^{#mu} (GeV/c)");
      //fraction_stack->GetYaxis()->SetTitle("Jet particle fraction");
      //fraction_stack->GetXaxis()->SetLimits(-1.5,1.5);
      //fraction_stack->GetXaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);

      //c1->SaveAs("muRelPtProjectionPlotV1-1_jetpt_50-100.pdf");
      //c1->SaveAs("muRelPtProjectionPlotV1-1_jetpt_100-150.pdf");
      //c1->SaveAs("muRelPtProjectionPlotV1-1_jetpt_150-200.pdf");
      //c1->SaveAs("muRelPtProjectionPlotV1-1_jetpt_200-250.pdf");
      //c1->SaveAs("muRelPtProjectionPlotV1-1_jetpt_250-300.pdf");
      //c1->SaveAs("muRelPtProjectionPlotV1-1_weights.pdf");

      //c1->SaveAs("muRelPtProjectionPlotV1-2_muptcut_5.pdf");
      //c1->SaveAs("muRelPtProjectionPlotV1-2_muptcut_12-5.pdf");
      //c1->SaveAs("muRelPtProjectionPlotV1-2_muptcut_16.pdf");
      c1->SaveAs("muRelPtProjectionPlotV1-2_muptcut_20.pdf");
      
      

     

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

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

      h5 -> SetMinimum(1*pow(10,1));
      h5 -> SetMaximum(1*pow(10,4));

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

      c2->SaveAs("etaProjectionPlotV2-2.pdf");

*/







} // END PROGRAM