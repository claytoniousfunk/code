

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
void muRelPtProjectionPlotV1(Float_t ptcut1=50.0, Float_t ptcut2=500.0){

/// LOAD DATA ///

//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_5.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_12-5.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_16.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_20.root");

TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_recojets_pthat_30_jetptcut_50_muptcut_2-5.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_recojets_pthat_30_jetptcut_50_muptcut_5.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_recojets_pthat_30_jetptcut_50_muptcut_7-5.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_recojets_pthat_30_jetptcut_50_muptcut_10.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_recojets_pthat_30_jetptcut_50_muptcut_12-5.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_recojets_pthat_30_jetptcut_50_muptcut_15.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_recojets_pthat_30_jetptcut_50_muptcut_17-5.root");
//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV9_recojets_pthat_30_jetptcut_50_muptcut_20.root");

TH2D* h2;
TH2D* h2_relPtJetPt;
TH2D* h2_relPtJetPt_g;
TH2D* h2_relPtJetPt_lq;
TH2D* h2_relPtJetPt_sall;
TH2D* h2_relPtJetPt_ball;
TH2D* h2_relPtJetPt_call;

f->GetObject("h2",h2);
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

TH1D* h_all = h2->ProjectionX("All jets",firstybin,lastybin);
TH1D* h0 = h2_relPtJetPt->ProjectionX("All muon jets",firstybin,lastybin);
TH1D* h1 = h2_relPtJetPt_g->ProjectionX("gluon jets",firstybin,lastybin);
TH1D* h2a = h2_relPtJetPt_lq->ProjectionX("light quark jets",firstybin,lastybin);
TH1D* h3 = h2_relPtJetPt_ball->ProjectionX("s & #bar{s}",firstybin,lastybin);
TH1D* h4 = h2_relPtJetPt_sall->ProjectionX("c & #bar{c}",firstybin,lastybin);
TH1D* h5 = h2_relPtJetPt_call->ProjectionX("b & #bar{b}",firstybin,lastybin);

cout << "Total number of jets = " << h_all->Integral() << endl;

// generate dN/dp_T plot

      // INCLUDE BIN ERROR!

    Double_t Nbin = h1->GetSize();

    for (int i=0;i<Nbin;i++){

      Double_t x = h1->GetBinWidth(i);

      Double_t y1 = h1->GetBinContent(i);
      Double_t y2a = h2a->GetBinContent(i);
      Double_t y3 = h3->GetBinContent(i);
      Double_t y4 = h4->GetBinContent(i);
      Double_t y5 = h5->GetBinContent(i);

      Double_t z1 = h1->GetBinError(i);
      Double_t z2a = h2a->GetBinError(i);
      Double_t z3 = h3->GetBinError(i);
      Double_t z4 = h4->GetBinError(i);
      Double_t z5 = h5->GetBinError(i);

      Double_t r1 = 0;
      Double_t r2a = 0;
      Double_t r3 = 0;
      Double_t r4 = 0;
      Double_t r5 = 0;

      Double_t e1 = 0;
      Double_t e2a = 0;
      Double_t e3 = 0;
      Double_t e4 = 0;
      Double_t e5 = 0;
      
      if(x!=0){
        r1 = y1/x;
        r2a = y2a/x;
        r3 = y3/x;
        r4 = y4/x;
        r5 = y5/x;

        e1 = z1/x;
        e2a = z2a/x;
        e3 = z3/x;
        e4 = z4/x;
        e5 = z5/x;

        h1->SetBinContent(i,r1);
        h2a->SetBinContent(i,r2a);
        h3->SetBinContent(i,r3);
        h4->SetBinContent(i,r4);
        h5->SetBinContent(i,r5);

        h1->SetBinError(i,e1);
        h2a->SetBinError(i,e2a);
        h3->SetBinError(i,e3);
        h4->SetBinError(i,e4);
        h5->SetBinError(i,e5);
       
      }
      
    }



// Define the Canvas

      TCanvas *c1 = new TCanvas("c", "canvas", 1200, 800);

      c1->cd();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      //pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      pad1->cd();
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

      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 5 GeV/c, 50 < p_{T}^{jet} < 500 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 12.5 GeV/c, 50 < p_{T}^{jet} < 500 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 16 GeV/c, 50 < p_{T}^{jet} < 500 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 2.5 GeV/c, 50 < p_{T}^{jet} < 500 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 5 GeV/c, 50 < p_{T}^{jet} < 500 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 7.5 GeV/c, 50 < p_{T}^{jet} < 500 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 10 GeV/c, 50 < p_{T}^{jet} < 500 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 12.5 GeV/c, 50 < p_{T}^{jet} < 500 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 15 GeV/c, 50 < p_{T}^{jet} < 500 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 17.5 GeV/c, 50 < p_{T}^{jet} < 500 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 20 GeV/c, 50 < p_{T}^{jet} < 500 GeV/c");

      fraction_stack->GetXaxis()->SetTitle("p_{T,rel}^{#mu} (GeV/c)");
      fraction_stack->GetYaxis()->SetTitle("#frac{dN}{dp_{T,rel}^{#mu}}");
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

      TPad *pad2 = new TPad("pad2","pad2",0.5,0.2,0.8,0.85);
      pad2->Draw();
      pad2->cd();

      TString d1;
      d1 = "#int #frac{dN^{g}}{dp_{T,rel}^{#mu}} dp_{T,rel}^{#mu} = ";
      d1 += h1->Integral();
      TLatex *l1 = new TLatex(0.0,1.0,d1.Data());
      l1->SetTextSize(0.05);
      l1->DrawLatex(0.2,0.1,d1.Data());

      TString d2;
      d2 = "#int #frac{dN^{u,d}}{dp_{T,rel}^{#mu}} dp_{T,rel}^{#mu} = ";
      d2 += h2a->Integral();
      TLatex *l2 = new TLatex(0.0,1.0,d2.Data());
      l2->SetTextSize(0.05);
      l2->DrawLatex(0.2,0.25,d2.Data());

      TString d3;
      d3 = "#int #frac{dN^{s}}{dp_{T,rel}^{#mu}} dp_{T,rel}^{#mu} = ";
      d3 += h4->Integral();
      TLatex *l3 = new TLatex(0.0,1.0,d3.Data());
      l3->SetTextSize(0.05);
      l3->DrawLatex(0.2,0.4,d3.Data());

      TString d4;
      d4 = "#int #frac{dN^{c}}{dp_{T,rel}^{#mu}} dp_{T,rel}^{#mu} = ";
      d4 += h5->Integral();
      TLatex *l4 = new TLatex(0.0,1.0,d4.Data());
      l4->SetTextSize(0.05);
      l4->DrawLatex(0.2,0.55,d4.Data());

      TString d5;
      d5 = "#int #frac{dN^{b}}{dp_{T,rel}^{#mu}} dp_{T,rel}^{#mu} = ";
      d5 += h3->Integral();
      TLatex *l5 = new TLatex(0.0,1.0,d5.Data());
      l5->SetTextSize(0.05);
      l5->DrawLatex(0.2,0.7,d5.Data());

      TString d6;
      d6 = "#int #frac{dN^{tot}}{dp_{T,rel}^{#mu}} dp_{T,rel}^{#mu} = ";
      d6 += h1->Integral() + h2a->Integral()+ h3->Integral() + h4->Integral()+h5->Integral();
      TLatex *l6 = new TLatex(0.0,1.0,d6.Data());
      l6->SetTextSize(0.05);
      l6->DrawLatex(0.2,0.85,d6.Data());
















      //c1->SaveAs("muRelPtProjectionPlotV1-2.pdf");
      
      

     

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