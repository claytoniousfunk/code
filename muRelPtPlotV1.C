

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

void muRelPtPlotV1(Double_t ptcut1=0, Double_t ptcut2 = 5.0){
	/// LOAD DATA ///
	//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7.root");
	//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_2-5.root");
	//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_5.root");
	//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_10.root");
	TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_muptcut_15.root");
	//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_jetptcut_50.root");
	//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_jetptcut_100.root");
	//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_jetptcut_200.root");
	//TFile* f =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_jetptcut_300.root");
	
	TH1D *h0, *h1, *h2a, *h3, *h4, *h5;
	f->GetObject("h_muRelPt",h0);
	f->GetObject("h_muRelPt_g",h1);
	f->GetObject("h_muRelPt_lq",h2a);
	f->GetObject("h_muRelPt_ball",h3);
	f->GetObject("h_muRelPt_sall",h4);
	f->GetObject("h_muRelPt_call",h5);

	TAxis *xaxis = h0->GetXaxis();
	Int_t firstxbin = xaxis->FindBin(ptcut1);
	cout << "First x bin = " << firstxbin << endl;
	Int_t lastxbin = xaxis->FindBin(ptcut2);
	cout << "Last x bin = " << lastxbin << endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Define the Canvas

      TCanvas *c1 = new TCanvas("c", "canvas", 1000, 800);
      c1->cd();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      THStack *fraction_stack = new THStack("fraction_stack","");

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
      fraction_stack->GetXaxis()->SetTitle("Muon rel-p_{T} (GeV/c)");
      //fraction_stack->GetYaxis()->SetTitle("Jet fraction");
      //fraction_stack->GetXaxis()->SetLimits(-1.5,1.5);
      //fraction_stack->GetXaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);
      
      //fraction_stack->SetTitle("");
      //c1->SaveAs("muRelPtPlotV1-1.pdf");

      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 2.5 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 5 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 10 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{#mu} = 15 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{jet} = 50 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{jet} = 100 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{jet} = 200 GeV/c");
      //fraction_stack->SetTitle("p_{T,cut}^{jet} = 300 GeV/c");

      //h1->SetTitle("p_{T,cut}^{#mu} = 2.5 GeV/c");
      //h1->SetTitle("p_{T,cut}^{#mu} = 5 GeV/c");
      //h1->SetTitle("p_{T,cut}^{#mu} = 10 GeV/c");
      h1->SetTitle("p_{T,cut}^{#mu} = 15 GeV/c");


      //c1->SaveAs("muRelPtPlotV1-1_muptcut_2-5.pdf");
      //c1->SaveAs("muRelPtPlotV1-1_muptcut_5.pdf");
      //c1->SaveAs("muRelPtPlotV1-1_muptcut_10.pdf");
      //c1->SaveAs("muRelPtPlotV1-1_muptcut_15.pdf");
      //c1->SaveAs("muRelPtPlotV1-1_jetptcut_50.pdf");
      //c1->SaveAs("muRelPtPlotV1-1_jetptcut_100.pdf");
      //c1->SaveAs("muRelPtPlotV1-1_jetptcut_200.pdf");
      //c1->SaveAs("muRelPtPlotV1-1_jetptcut_300.pdf");

      //c1->SaveAs("muRelPtPlotV1-2_muptcut_2-5.pdf");
      //c1->SaveAs("muRelPtPlotV1-2_muptcut_5.pdf");
      //c1->SaveAs("muRelPtPlotV1-2_muptcut_10.pdf");
      //c1->SaveAs("muRelPtPlotV1-2_muptcut_15.pdf");
      //c1->SaveAs("muRelPtPlotV1-2_jetptcut_50.pdf");
      //c1->SaveAs("muRelPtPlotV1-2_jetptcut_100.pdf");
      //c1->SaveAs("muRelPtPlotV1-2_jetptcut_200.pdf");
      //c1->SaveAs("muRelPtPlotV1-2_jetptcut_300.pdf");

      //c1->SaveAs("muRelPtPlotV1-3_muptcut_2-5.pdf");
      //c1->SaveAs("muRelPtPlotV1-3_muptcut_5.pdf");
      //c1->SaveAs("muRelPtPlotV1-3_muptcut_10.pdf");
      c1->SaveAs("muRelPtPlotV1-3_muptcut_15.pdf");
      //c1->SaveAs("muRelPtPlotV1-3_jetptcut_50.pdf");
      //c1->SaveAs("muRelPtPlotV1-3_jetptcut_100.pdf");
      //c1->SaveAs("muRelPtPlotV1-3_jetptcut_200.pdf");
      //c1->SaveAs("muRelPtPlotV1-3_jetptcut_300.pdf");




      ///////////////////////////////// purity calculation /////////////////////////////////
      /*
      Double_t N=0;
      Double_t x_g=0;
      Double_t x_lq=0;
      Double_t x_sall=0;
      Double_t x_ball=0;
      Double_t x_call=0;

      N = ptcut2-ptcut1; // width of fraction plot
      x_g = (h1->Integral("width"))/N;
      x_lq = (h2a->Integral("width"))/N;
      x_ball = (h3->Integral("width"))/N;
      x_sall = (h4->Integral("width"))/N;
      x_call = (h5->Integral("width"))/N;

      cout<<"Gluon purity = "<< x_g<<endl;
      cout<<"Light quark purity= "<<x_lq<<endl;
      cout<<"b puritiy = "<<x_ball<<endl;
      cout<<"s purity = "<<x_sall<<endl;
      cout<<"c purity = "<<x_call<<endl;
	  */






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}  // end program