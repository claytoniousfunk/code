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

void normalizeAndCompareSpectra(){

	TFile* fc =TFile::Open("/home/clayton/Analysis/code/etaPtHistoV7_pthat_50.root"); // clayton's file
	TFile* fx =TFile::Open("/home/clayton/Analysis/xiaoSpectra/bjtc_djetMC_format2_std_v2_pthat50.root"); // xiao's file


	TH2D *h2c;
	TH1D *hx, *hc;
	fc->GetObject("h2",h2c);
	hc = h2c->ProjectionY();
	fx->GetObject("jetQASets/incl_RecoLevel_pt_C0",hx);


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////// DEFINE BINNING ///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	/////////////////////////////////////////////////////////////////// create pT edges array /////////////////////////////////////////////////////////////////////////
	
	const Int_t Nbins = 1000;
	Double_t ptEdges[Nbins+1];
	Double_t deltaPt = 10.0; // initial width of bin (pT) // 10 GeV up to 100, 20 GeV up to 250, 250 out = 50 GeV, 400-600 = one bin
	Double_t ptmin=120.0;
	Double_t ptmax=500.0;
	Double_t jtpt_check=ptmin;
	ptEdges[0]=ptmin;
	int numUsedPtBins = 0;
	for(int yy=1; jtpt_check<ptmax; yy++){

		if(jtpt_check<100){
			deltaPt=10.;	
			ptEdges[yy]=ptEdges[yy-1]+deltaPt;
			numUsedPtBins=numUsedPtBins+1;
		}
		else if (jtpt_check<250){
			deltaPt=20.;
			ptEdges[yy]=ptEdges[yy-1]+deltaPt;
			numUsedPtBins=numUsedPtBins+1;
		}
		else if (jtpt_check<400){
			deltaPt=50;
			ptEdges[yy]=ptEdges[yy-1]+deltaPt;
			numUsedPtBins=numUsedPtBins+1;
			if(ptEdges[yy]>=400){
				ptEdges[yy]=400;
				ptEdges[yy+1]=ptmax;
				numUsedPtBins=numUsedPtBins+1;
				break;
			}
		}
		
		
		
		jtpt_check=jtpt_check+deltaPt;
		
	}
	

	Double_t usedPtEdges[numUsedPtBins+1];
	for(int ii=0; ii<numUsedPtBins+1; ii++){
		usedPtEdges[ii]=ptEdges[ii];
		cout << usedPtEdges[ii] << endl;
	}

	// REBIN!
	//hc-> Rebin(numUsedPtBins,"",usedPtEdges);
	//hx-> Rebin(numUsedPtBins,"",usedPtEdges);


	// INCLUDE BIN ERROR!
	

	// normalize hc

    Double_t NbinC = hc->GetSize();

    for (int i=0;i<NbinC;i++){

      Double_t xc = hc->GetBinWidth(i);

      Double_t yc = hc->GetBinContent(i);
     
      Double_t zc = hc->GetBinError(i);
       
      Double_t rc = 0;

      Double_t ec = 0;

      
      if(xc!=0){

        rc = yc/xc;    
        ec = zc/xc;
        hc->SetBinContent(i,rc); 
        hc->SetBinError(i,ec);
 
               
      }
      
    }
    
    

// normalize hx

   Double_t NbinX = hx->GetSize();

    for (int j=0;j<NbinX;j++){

      Double_t xx = hx->GetBinWidth(j);

      Double_t yx = hx->GetBinContent(j);
     
      Double_t zx = hx->GetBinError(j);
       
      Double_t rx = 0;

      Double_t ex = 0;

      
      if(xx!=0){

        rx = yx/xx;    
        ex = zx/xx;
        hx->SetBinContent(j,rx); 
        hx->SetBinError(j,ex);
 
               
      }
      
    }
    


// Define the Canvas
      TCanvas *c2 = new TCanvas("c2", "canvas", 1200, 800);
      c2->cd();
      c2->SetLogy();
      TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1, 1);
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1

      // jet pT 
      TH1D *h1 = (TH1D*)hc->Clone("h1"); // clayton
      TH1D *h2a = (TH1D*)hx->Clone("h2a"); // xiao
      
      
      Double_t msize = 1.2;

      h1 -> SetMinimum(1*pow(10,-5));
      h1 -> SetMaximum(1*pow(10,4));

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
    

      h1->Draw("e2");
      h2a->Draw("e2 same");
      

      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      legend->AddEntry(h1,"Clayton","fp");
      legend->AddEntry(h2a,"Xiao","fp");
      legend->Draw();
      
      
      //fraction_stack->GetYaxis()->SetTitleSize(18);
      //fraction_stack->GetXaxis()->SetTitleFont(43);
      //fraction_stack->GetXaxis()->SetTitleOffset(3.);
      //fraction_stack->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      //fraction_stack->GetXaxis()->SetLabelSize(12);
      h1->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      h1->GetYaxis()->SetTitle("#frac{dN^{jet}}{dp_{T}}");
      
      

      //h5->SetTitle("");
      //c2->SaveAs("ptProjectionPlotV1-2.pdf");

      //h5->SetTitle("Ref jets");
      //c2->SaveAs("ptProjectionPlotV1-2_ref.pdf");
      
      //h5->SetTitle("Reco jets");
      //c2->SaveAs("ptProjectionPlotV1-2_reco.pdf");

      h1->SetTitle("");
      c2->SaveAs("normalizeAndCompareSpectra_claytonVsXiao.pdf");













}