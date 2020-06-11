
// V2: muon jets

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
void ptProjectionPlotV2(Double_t etacut1=-1.5, Double_t etacut2=1.5){

/// LOAD DATA ///
TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7.root");

TH2D* h2_MJ;
TH2D* h2_d_MJ;
TH2D* h2_u_MJ;
TH2D* h2_s_MJ;
TH2D* h2_c_MJ;
TH2D* h2_b_MJ;
TH2D* h2_dbar_MJ;
TH2D* h2_ubar_MJ;
TH2D* h2_sbar_MJ;
TH2D* h2_cbar_MJ;
TH2D* h2_bbar_MJ;
TH2D* h2_hq_MJ;
TH2D* h2_lq_MJ;
TH2D* h2_qall_MJ;
TH2D* h2_dall_MJ;
TH2D* h2_uall_MJ;
TH2D* h2_sall_MJ;
TH2D* h2_call_MJ;
TH2D* h2_ball_MJ;
TH2D* h2_g_MJ;
TH2D* h2_ud_MJ;
TH2D* h2_ubardbar_MJ;
TH2D* h2_other_MJ;

f->GetObject("h2_MJ",h2_MJ);
f->GetObject("h2_d_MJ",h2_d_MJ);
f->GetObject("h2_u_MJ",h2_u_MJ);
f->GetObject("h2_s_MJ",h2_s_MJ);
f->GetObject("h2_c_MJ",h2_c_MJ);
f->GetObject("h2_b_MJ",h2_b_MJ);
f->GetObject("h2_dbar_MJ",h2_dbar_MJ);
f->GetObject("h2_ubar_MJ",h2_ubar_MJ);
f->GetObject("h2_sbar_MJ",h2_sbar_MJ);
f->GetObject("h2_cbar_MJ",h2_cbar_MJ);
f->GetObject("h2_bbar_MJ",h2_bbar_MJ);
f->GetObject("h2_hq_MJ",h2_hq_MJ);
f->GetObject("h2_lq_MJ",h2_lq_MJ);
f->GetObject("h2_qall_MJ",h2_qall_MJ);
f->GetObject("h2_dall_MJ",h2_dall_MJ);
f->GetObject("h2_uall_MJ",h2_uall_MJ);
f->GetObject("h2_sall_MJ",h2_sall_MJ);
f->GetObject("h2_call_MJ",h2_call_MJ);
f->GetObject("h2_ball_MJ",h2_ball_MJ);
f->GetObject("h2_g_MJ",h2_g_MJ);
f->GetObject("h2_ud_MJ",h2_ud_MJ);
f->GetObject("h2_ubardbar_MJ",h2_ubardbar_MJ);
f->GetObject("h2_other_MJ",h2_other_MJ);



TH1D* hx_test=h2_MJ->ProjectionX();
TAxis *xaxis = hx_test->GetXaxis();
Int_t firstxbin = xaxis->FindBin(etacut1);
cout << firstxbin << endl;
Int_t lastxbin = xaxis->FindBin(etacut2);
cout << lastxbin << endl;

TH1D* h_jetpt = h2_MJ->ProjectionY("All jets",firstxbin,lastxbin);
TH1D* h_jetpt_g = h2_g_MJ->ProjectionY("gluon jets",firstxbin,lastxbin);
TH1D* h_jetpt_lq = h2_lq_MJ->ProjectionY("light quark jets",firstxbin,lastxbin);
TH1D* h_jetpt_sall = h2_sall_MJ->ProjectionY("s & #bar{s}",firstxbin,lastxbin);
TH1D* h_jetpt_call = h2_call_MJ->ProjectionY("c & #bar{c}",firstxbin,lastxbin);
TH1D* h_jetpt_ball = h2_ball_MJ->ProjectionY("b & #bar{b}",firstxbin,lastxbin);
TH1D* h_jetpt_other = h2_other_MJ->ProjectionY("ghost",firstxbin,lastxbin);



// Define the Canvas

      TCanvas *c1 = new TCanvas("c", "canvas", 1000, 800);
      c1->cd();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      THStack *fraction_stack = new THStack("fraction_stack","");
      // jet pT stack
      TH1D *h0 = (TH1D*)h_jetpt->Clone("h0"); // all jets
      TH1D *h1 = (TH1D*)h_jetpt_g->Clone("h1"); // gluons
      TH1D *h2a = (TH1D*)h_jetpt_lq->Clone("h2a"); // light quarks
      TH1D *h3 = (TH1D*)h_jetpt_ball->Clone("h3"); // b and bbar
      TH1D *h4 = (TH1D*)h_jetpt_sall->Clone("h4"); // s and sbar
      TH1D *h5 = (TH1D*)h_jetpt_call->Clone("h5"); // c and cbar
      TH1D *h6 = (TH1D*)h_jetpt_other->Clone("h6"); // ghost

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

      h6 -> SetMarkerStyle(29);
      h6 -> SetMarkerColor(kGray);
      h6 -> SetMarkerSize(msize);
      h6 -> SetFillColorAlpha(kBlack,0.5);
      h6 -> SetStats(0);
      h6 -> Divide(h0);
     
     
           
      fraction_stack->Add(h5);
      fraction_stack->Add(h4);
      fraction_stack->Add(h3);
      fraction_stack->Add(h2a);
      fraction_stack->Add(h1);
      fraction_stack->Add(h6); 



      fraction_stack->Draw("e1 hist");
      auto legend = new TLegend(0.905,0.4,0.995,0.7);

      legend->AddEntry(h6,"?","fp");
      legend->AddEntry(h1,"g","fp");
      legend->AddEntry(h2a,"u, d","fp");
      legend->AddEntry(h3,"b","fp");
      legend->AddEntry(h4,"s","fp");
      legend->AddEntry(h5,"c","fp");
      legend->Draw();
      fraction_stack->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      fraction_stack->GetYaxis()->SetTitle("Jet fraction");
      fraction_stack->SetTitle("");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 2.5 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 5 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 10 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 15 GeV/c");
      //fraction_stack->GetXaxis()->SetLimits(-1.5,1.5);
      //fraction_stack->GetXaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);
      c1->SaveAs("ptProjectionPlotV2-1.pdf");


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
      TH1D *h1 = (TH1D*)h_jetpt_g->Clone("h1"); // gluons
      TH1D *h2a = (TH1D*)h_jetpt_lq->Clone("h2a"); // light quarks
      TH1D *h3 = (TH1D*)h_jetpt_ball->Clone("h3"); // b and bbar
      TH1D *h4 = (TH1D*)h_jetpt_sall->Clone("h4"); // s and sbar
      TH1D *h5 = (TH1D*)h_jetpt_call->Clone("h5"); // c and cbar
      
      Double_t msize = 1.2;

      h5 -> SetMinimum(1*pow(10,0));
      h5 -> SetMaximum(1*pow(10,3));

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
    
      
      h5->Draw("e2");
      h4->Draw("same e2");
      h3->Draw("same e2");
      h2a->Draw("same e2");
      h1->Draw("same e2");

      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      legend->AddEntry(h1,"g","fp");
      legend->AddEntry(h2a,"u,d","fp");
      legend->AddEntry(h3,"b","fp");
      legend->AddEntry(h4,"s","fp");
      legend->AddEntry(h5,"c","fp");
      legend->Draw();
      
      //fraction_stack->GetYaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);
      h5->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      h5->GetYaxis()->SetTitle("#frac{dN^{jet}}{dp_{T}}");
      h5->SetTitle("");

      c2->SaveAs("ptProjectionPlotV2-2.pdf");
*/







} // END PROGRAM