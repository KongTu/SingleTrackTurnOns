#include "RiceStyle.h"

using namespace std;

void L1plusHLTSingleTrackTurnOns(){



	TH1D* hist = makeHist("hist","", "p_{T} (GeV)", "Efficiency", 100,0,100);

	TFile * file[3];
	// MinBias L1 + HLT SingleTrack for PbPb:
	file[0] = new TFile("../rootfiles/PbPb_SingleTrackTurnOns_v2.root");

	TH1D* allEvents1;
	TH1D* allEvents2;
	TH1D* triggerPath[10];
	TGraphAsymmErrors * gr[5];

	vector<TH1D*> spectra = loadingHistogram(file[0], "ana_PbPb", 10);
	int binwidth = 50;

	TCanvas* c2 = makeMultiCanvas("c2","",3,2);

	for( int i = 0; i < 5; i++ ){

		c2->cd(i+1);
		gPad->SetTicks();
	 	gPad->SetLeftMargin(0.14);
 		gPad->SetBottomMargin(0.14);
		hist->Draw();
		gr[i] = new TGraphAsymmErrors();

		spectra[i]->Rebin(binwidth);
		spectra[i+5]->Rebin(binwidth);

		gr[i]->Divide(spectra[i+5], spectra[i], "cp");

		if( i == 4 ){
			gr[i]->SetLineColor(6);
			gr[i]->SetMarkerColor(6);
		}
		else{
			gr[i]->SetLineColor(i+1);
			gr[i]->SetMarkerColor(i+1);
		}
		gr[i]->SetMarkerStyle(20);
		gr[i]->Draw("Psame");

	}
	
	c2->cd(6);

	TLegend *w1 = new TLegend(0.10,0.20,0.90,0.80);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->AddEntry(gr[0],"L1 + HLT SingleTrack 12","P");
    w1->AddEntry(gr[1],"L1 + HLT SingleTrack 18","P");
    w1->AddEntry(gr[2],"L1 + HLT SingleTrack 24","P");
    w1->AddEntry(gr[3],"L1 + HLT SingleTrack 34","P");
    w1->AddEntry(gr[4],"L1 + HLT SingleTrack 45","P");
    w1->Draw("");

	// c2->SaveAs("../files/triggerEfficiency.png");
 // 	c2->SaveAs("../files/triggerEfficiency.pdf");



}
