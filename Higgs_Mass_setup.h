#include "TFile.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooCrystalBall.h"
// #include "RooMyPDF_DSCB.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "RooVoigtian.h"
#include <fstream>
#include "RooLandau.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TRandom.h"
#include "TRandom3.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

// using namespace TMVA;
using namespace RooFit;


TH1F *h_blank = new TH1F("h_blank", "h_blank", 10, -0.5, 9.5);

TH1F* h_num_eventi;
float num_event;
Long64_t nentries;

TH2F *hMuScaleFac;
TH2F *hMuScaleFacUnc;

TH1F* h_mu_B;
TH1F* h_mu_E;
TH1F* h_el_B;
TH1F* h_el_E;
TAxis* x_pTaxis;

ULong64_t Run, Event, LumiSect;

Double_t pT1_FromMuonBestTrack, pT2_FromMuonBestTrack;
Double_t eta1_FromMuonBestTrack, eta2_FromMuonBestTrack;
Double_t phi1_FromMuonBestTrack, phi2_FromMuonBestTrack;
Double_t muonPV_x1, muonPV_x2;
Double_t muonPV_y1, muonPV_y2;
Double_t muonPV_z1, muonPV_z2;
Double_t m1, pT1, eta1, phi1, pterr1old, Iso1;
Double_t pT_FSR1, eta_FSR1, phi_FSR1, m_FSR1;
Int_t Id1;
Double_t d0BS1, d0BS2;
Double_t d0PV1, d0PV2;
Double_t dzPV1, dzPV2;
Double_t Sip1, Sip2;
Double_t genLep_pt1,genLep_pt2;
Double_t pterr1, pterr1_single, pterr1_VX, pterr1_VX_BS;
Double_t genLep_eta1, genLep_eta2;
Double_t genLep_phi1, genLep_phi2;
Double_t m2, pT2, eta2, phi2, pterr2old, Iso2;
Double_t pT_FSR2, eta_FSR2, phi_FSR2, m_FSR2;
Int_t Id2;
Double_t pterr2, pterr2_single, pterr2_VX, pterr2_VX_BS;
// Double_t weight;
Float_t weight;
Float_t k_ggZZ;
Float_t k_qqZZ_ewk;
Float_t k_qqZZ_qcd_M;
Int_t lep1_ecalDriven;
Int_t lep2_ecalDriven;
Double_t massZ, massZErr;
Double_t massZ_FSR, massZErr_FSR;
Double_t genzm;
Double_t GENmass2l;
float GENmass4l;
Int_t nZXCRFailedLeptons;


Float_t genWeight, pileupWeight, prefiringWeight, dataMCWeight,pileupWeightUp, pileupWeightDn;

Double_t massZ_single_BS, massZ_single_BS_FSR, massZ_single_BS_BS, massZ_single_BS_BS_FSR;
Double_t massZ_vtx, massZ_vtx_FSR, massZ_vtx_BS, massZ_vtx_BS_FSR;
Double_t massZErr_single_BS, massZErr_single_BS_FSR, massZErr_vtx, massZErr_vtx_FSR, massZErr_vtx_BS, massZErr_vtx_BS_FSR;
Double_t single_BS_pT1, single_BS_pT2, single_BS_eta1, single_BS_eta2, single_BS_m1, single_BS_m2, single_BS_phi1, single_BS_phi2;
Double_t vtx_pT1, vtx_pT2, vtx_eta1, vtx_eta2, vtx_m1, vtx_m2, vtx_phi1, vtx_phi2;
Double_t vtx_BS_pT1, vtx_BS_pT2, vtx_BS_eta1, vtx_BS_eta2, vtx_BS_m1, vtx_BS_m2, vtx_BS_phi1, vtx_BS_phi2;
Double_t single_BS_pT_FSR1, single_BS_pT_FSR2, single_BS_eta_FSR1, single_BS_eta_FSR2, single_BS_m_FSR1, single_BS_m_FSR2, single_BS_phi_FSR1, single_BS_phi_FSR2;
Double_t vtx_pT_FSR1, vtx_pT_FSR2, vtx_eta_FSR1, vtx_eta_FSR2, vtx_m_FSR1, vtx_m_FSR2, vtx_phi_FSR1, vtx_phi_FSR2;
Double_t vtx_BS_pT_FSR1, vtx_BS_pT_FSR2, vtx_BS_eta_FSR1, vtx_BS_eta_FSR2, vtx_BS_m_FSR1, vtx_BS_m_FSR2, vtx_BS_phi_FSR1, vtx_BS_phi_FSR2;
Double_t d0BS_vtx_BS1, d0BS_vtx_BS2;
Double_t pT1_genFromReco, pT2_genFromReco;
Double_t Tracker1, Tracker2;
Int_t Tight1, Tight2;
float massZ2;

std::vector<int> *lep_id = 0; 
std::vector<float> *lep_tightId = 0;
std::vector<float> *lep_RelIso = 0;
std::vector<float> *GENlep_hasFSR = 0;
std::vector<float> *GENlep_status = 0;
std::vector<float> *GENlep_MomId = 0;
std::vector<float> *GENlep_MomMomId = 0;
std::vector<float> *GENlep_id = 0;
std::vector<float> *GENlep_mass = 0;
std::vector<float> *GENlep_pt = 0; 
std::vector<float> *pho_pt = 0;
std::vector<float> *fsrPhotons_pt = 0;
std::vector<float> *lep_pt = 0;
std::vector<float> *lep_pt_UnS = 0;
std::vector<float> *lep_pt_genFromReco = 0;
std::vector<float> *lep_pt_FromMuonBestTrack = 0;
std::vector<float> *lep_trackerLayersWithMeasurement = 0;
std::vector<float> *lepFSR_pt = 0;
std::vector<float> *vtxLep_pt = 0;
std::vector<float> *vtxLep_BS_pt_NoRoch = 0;
std::vector<float> *vtxLep_BS_pt = 0;
std::vector<float> *vtxLepFSR_pt = 0;
std::vector<float> *vtxLepFSR_BS_pt = 0;
std::vector<float> *GENlep_eta = 0;
std::vector<float> *lep_eta = 0;
std::vector<float> *pho_eta = 0;
std::vector<float> *fsrPhotons_eta = 0;
std::vector<float> *lepFSR_eta = 0;
std::vector<float> *vtxLep_eta = 0;
std::vector<float> *vtxLep_BS_eta = 0;
std::vector<float> *vtxLepFSR_eta = 0;
std::vector<float> *vtxLepFSR_BS_eta = 0;
std::vector<float> *GENlep_phi = 0;
std::vector<float> *lep_phi = 0;
std::vector<float> *pho_phi = 0;
std::vector<float> *fsrPhotons_phi = 0;
std::vector<float> *lepFSR_phi = 0;
std::vector<float> *vtxLep_phi = 0;
std::vector<float> *vtxLep_BS_phi = 0;
std::vector<float> *vtxLepFSR_phi = 0;
std::vector<float> *vtxLepFSR_BS_phi = 0;
std::vector<float> *lep_mass = 0;
std::vector<float> *lepFSR_mass = 0;
std::vector<float> *vtxLep_mass = 0;
std::vector<float> *vtxLep_BS_mass = 0;
std::vector<float> *vtxLepFSR_mass = 0;
std::vector<float> *vtxLepFSR_BS_mass = 0;
std::vector<int> *lep_numberOfValidPixelHits = 0;
std::vector<float> *vtxLep_BS_d0 = 0;
std::vector<float> *lep_d0BS = 0;
std::vector<int> *lep_ecalDriven = 0;

std::vector<float> *jet_pt = 0;
std::vector<float> *jet_eta = 0;
std::vector<float> *jet_phi = 0;
std::vector<float> *jet_mass = 0;
std::vector<float> *jet_csv_cTag_vsL = 0;
std::vector<float> *jet_csv_cTag_vsB = 0;
std::vector<float> *jet_partonFlavour = 0;
std::vector<float> *jet_QGTagger = 0;

Int_t njets_pt30_eta4p7;

int lep_Hindex[4];


std::vector<float> *lep_dataMC = 0;
std::vector<float> *lep_dataMCErr = 0;
std::vector<float> *lep_normalizedChi2 = 0;
std::vector<float> *lep_numberOfValidMuonHits = 0;
std::vector<float> *lep_numberOfMatchedStations = 0;
std::vector<float> *lep_Sip = 0;
std::vector<float> *lep_d0PV = 0;
std::vector<float> *lep_dzPV = 0;
	
std::vector<float> *commonPV_x = 0;
std::vector<float> *commonPV_y = 0;
std::vector<float> *commonPV_z = 0;
std::vector<double> *commonBS_x = 0;
std::vector<double> *commonBS_y = 0;
std::vector<double> *commonBS_z = 0;
	
float PV_x, PV_y, PV_z;
float BS_x, BS_y, BS_z;
float BS_xErr, BS_yErr, BS_zErr;
float BeamWidth_x, BeamWidth_y;
	
float pTL1, pTL2, pTL3, pTL4;
float etaL1, etaL2, etaL3, etaL4;
float phiL1, phiL2, phiL3, phiL4;
float mL1, mL2, mL3, mL4;	

float met;

float GENMH;

int finalState;
bool passedFullSelection;
Int_t is2P2F, is3P1F, isMCzz;
Float_t eventWeightFR;
Float_t eventWeightFR_up;
Float_t eventWeightFR_down;
Float_t fr3;
Float_t fr2;
float mass4l, mass4lREFIT, mass4mu, mass2e2mu, mass4l_vtx, mass4l_vtxFSR, mass4lREFIT_vtx, mass4l_vtx_BS, mass4l_vtxFSR_BS, mass4lREFIT_vtx_BS;
float mass4lErr, mass4lErrREFIT, mass4lErr_vtx, mass4lErr_vtx_BS, mass4lErrREFIT_vtx, mass4lErrREFIT_vtx_BS;
float D_bkg_kin, D_bkg_kin_vtx_BS;
float massZ1, massZ1REFIT;
float GENmassZ1;

std::vector<float> *lep_pterr = 0;
std::vector<float> *vtxLep_ptError = 0;
std::vector<float> *vtxLep_BS_ptError = 0;


int nVtx, nInt;
int NVtx, NInt;



TString directory;
TString execute;
TString filename;
TString save_nome;

TCanvas* c1;


	std::vector<std::vector< Double_t >> prova_2;
	std::vector<Double_t> prova_4mu;
	std::vector<Double_t> prova_4e;
	std::vector<Double_t> prova_2e2mu;
	std::vector<Double_t> prova_2mu2e;


std::vector<Double_t> pT_bins;
std::vector<Double_t> pT_bins_err;
std::vector<Double_t> pT_bins_center;
std::vector<Double_t> eta_bins;
std::vector<TString> eta_bins_name;
std::vector<Double_t> eta_bins_count;
std::vector<Double_t> phi_bins;
// std::vector<Double_t> pT_bins_graph;
// std::vector<Double_t> pT_err_bins_graph;

std::vector<Double_t> MC_mean_vs_gamma_DSCB;
std::vector<Double_t> MC_meanErr_vs_gamma_DSCB;
std::vector<Double_t> DATA_mean_vs_gamma_DSCB;
std::vector<Double_t> DATA_meanErr_vs_gamma_DSCB;

std::vector<Double_t> MC_mean_vs_gamma_CB;
std::vector<Double_t> MC_meanErr_vs_gamma_CB;
std::vector<Double_t> DATA_mean_vs_gamma_CB;
std::vector<Double_t> DATA_meanErr_vs_gamma_CB;

std::vector<Double_t> MC_mean_vs_pt_DSCB;
std::vector<Double_t> MC_meanErr_vs_pt_DSCB;
std::vector<Double_t> DATA_mean_vs_pt_DSCB;
std::vector<Double_t> DATA_meanErr_vs_pt_DSCB;

std::vector<Double_t> MC_mean_vs_pt_CB;
std::vector<Double_t> MC_meanErr_vs_pt_CB;
std::vector<Double_t> DATA_mean_vs_pt_CB;
std::vector<Double_t> DATA_meanErr_vs_pt_CB;

std::vector<Double_t> scale_DATA_MC_DSCB;
std::vector<Double_t> scaleErr_DATA_MC_DSCB;
std::vector<Double_t> scale_DATA_MC_CB;
std::vector<Double_t> scaleErr_DATA_MC_CB;

std::vector<Double_t> spreadScale_DSCB;
std::vector<Double_t> spreadScaleErr_DSCB;
std::vector<Double_t> spreadScale_CB;
std::vector<Double_t> spreadErrScale_CB;

std::vector<Double_t> spreadMean_DSCB;
std::vector<Double_t> spreadMeanErr_DSCB;
std::vector<Double_t> spreadMean_CB;
std::vector<Double_t> spreadErrMean_CB;


float MC_mean, DATA_mean;
float MC_meanErr, DATA_meanErr;
float spreadScale_DSCB_min = 999; 
float spreadScale_DSCB_max = 0;
float spreadMean_DSCB_min = 999; 
float spreadMean_DSCB_max = 0;

float spreadScale_CB_min = 999; 
float spreadScale_CB_max = 0;
float spreadMean_CB_min = 999; 
float spreadMean_CB_max = 0;


TH2F *hist_lut;


bool DATA;

TString fs;
TString mode;
TString histo_name;
TString fit_param_latex;
TString canvas_name;

int chi_dof;

TLine* lineRef = new TLine(60,0,120,0.);
TLine* lineRef_H = new TLine(105,0,140,0.);

TFile* _file0;
TTree *tree;                                             

TLorentzVector lep_1, lep_2;
TLorentzVector lep_1err, lep_2err;
TLorentzVector ZPrime; 


TLine* scale_up;
TLine* scale_down;

std::vector<TString> mass_type;

std::vector<Double_t> gamma_bins;
std::vector<Double_t> gammaErr_bins;		

std::vector<int> massPoint;
std::vector<TString> massPointName;
std::vector<TString> ProdMode;
std::vector<TString> ProdMode_File;
std::vector<TString> year;
std::vector<TString> category;
std::vector<TString> decayMode;
std::vector<float> luminosity;

std::vector<int> numberEvent_ggH;
std::vector<int> numberEvent_VBF;
std::vector<int> numberEvent_WplusH;
std::vector<int> numberEvent_WminusH;
std::vector<int> numberEvent_ZH;
std::vector<int> numberEvent_ttH;
std::vector<std::vector<int> > numberEvent;	
std::vector<int> numberEvent_2018;

std::vector<float> crossSection_ggH;
std::vector<float> crossSection_VBF;
std::vector<float> crossSection_WplusH;
std::vector<float> crossSection_WminusH;
std::vector<float> crossSection_ZH;
std::vector<float> crossSection_ttH;
std::vector<float> crossSection_bbH;
std::vector<float> crossSection_tHq;

std::vector<std::vector<float> > crossSection;	

Color_t color;
TString nome_file;




float x_min = 60;
float x_max = 120;

float massZ_min = 60;
float massZ_max = 120;
float massZErr_min = 0.2;
float massZErr_max = 7.2;

float BW_mean_PDG = 91.1876;
float BW_sigma_PDG = 2.4952;

float BW_mean_min = 86;
float BW_mean_max = 96;
float BW_mean_DATA_min = 86;
float BW_mean_DATA_max = 96;
float BW_mean_MC_min = 86;
float BW_mean_MC_max = 96;

float CB_mean_min = -5;
float CB_mean_max = 5;
float CB_sigma_min = 0;
float CB_sigma_max = 30;
float CB_alpha_min = 0;
float CB_alpha_max = 30;
float CB_exp_min = 0;
float CB_exp_max = 30;

float DSCB_mean_min = -5;
float DSCB_mean_max = 5;
float DSCB_sigma_min = 0;
float DSCB_sigma_max = 5;
float DSCB_alphaL_min = 0;
float DSCB_alphaL_max = 20;
float DSCB_expL_min = 0;
float DSCB_expL_max = 20;
float DSCB_alphaR_min = 0;
float DSCB_alphaR_max = 10;
float DSCB_expR_min = 0;
float DSCB_expR_max = 10;

float tau_min = -1;
float tau_max = 1;	
				
float fsig_min = 0.01;
float fsig_max = 1;


void DrawResolution(TH1F* h1, TH1F* h2, TH1F* h3, std::vector<TString> histos_name, TString nome_canvas, TString save, TString x_name);
void DrawResolution(TH1F* h1, TH1F* h2, TH1F* h3, std::vector<TString> histos_name, TString nome_canvas, TString save, TString x_name){

	h1->SetStats(0);
	h2->SetStats(0);
	h3->SetStats(0);

// 	h1->Scale(1/h1->Integral());
// 	h2->Scale(1/h2->Integral());
// 	h3->Scale(1/h3->Integral());
	
	float max = 1;
// 	std::cout<<h1->GetTitle()<<"\t"<<h2->GetTitle()<<"\t"<<h3->GetTitle();
// 	std::cout<<max<<"\t"<<h1->GetMaximum()<<"\t"<<h2->GetMaximum()<<"\t"<<h3->GetMaximum();

	max = h1->GetMaximum();
	if(h2->GetMaximum() > max)
		max = h2->GetMaximum();
	if(h3->GetMaximum() > max)
		max = h3->GetMaximum();
	
	std::cout<<"\t max final = "<<max<<std::endl;

	max = 1.1*max;
	
	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11;
   	pad11 = new TPad("pad1", "pad1", 0, 0.15, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	
	h1->SetTitle(nome_canvas);
	h1->SetLineColor(kBlack);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(0.75);
	h1->SetMarkerColor(kBlack);
	h1->GetYaxis()->SetRangeUser(0, max);
// 	h1->GetYaxis()->SetRangeUser(0, 0.035);
// 	h1->GetYaxis()->SetTitle("#sigma");
	h1->GetXaxis()->SetTitle(x_name);

	h2->SetTitle(nome_canvas);
	h2->SetLineColor(kRed);
	h2->SetMarkerStyle(20);
	h2->SetMarkerSize(0.75);
	h2->SetMarkerColor(kRed);
	h2->GetYaxis()->SetRangeUser(0, max);
// 	h2->GetYaxis()->SetRangeUser(0, 0.035);
// 	h2->GetYaxis()->SetTitle("#sigma");
	h2->GetXaxis()->SetTitle(x_name);

	h3->SetTitle(nome_canvas);
	h3->SetLineColor(kGreen);
	h3->SetMarkerStyle(20);
	h3->SetMarkerSize(0.75);
	h3->SetMarkerColor(kGreen);
	h3->GetYaxis()->SetRangeUser(0, max);
// 	h3->GetYaxis()->SetRangeUser(0, 0.035);
// 	h3->GetYaxis()->SetTitle("#sigma");
	h3->GetXaxis()->SetTitle(x_name);

	TH1F *ratio_2 = (TH1F*) h2->Clone();
	TH1F *ratio_3 = (TH1F*) h3->Clone();
	TLegend *legend = new TLegend(0.65,0.15,0.9,0.35);
	
	canvas->Update();
	
	h1->Draw("E");
	h2->Draw("same E");
	h3->Draw("same E");

	legend->AddEntry(h1, histos_name.at(0));
	legend->AddEntry(h2, histos_name.at(1));
	legend->AddEntry(h3, histos_name.at(2));
	legend->Draw();
	
	
// 	float mean_1 = h1->GetMean();
// 	TLatex *tex1;
// 	TString fit_param_latex = Form("mean_1 =  %f", mean_1);
// 	tex1 = new TLatex(0.5,0.8, fit_param_latex);
// 	tex1->SetNDC();
// 	tex1->Draw();
// 
// 	mean_1 = h2->GetMean();
// 	fit_param_latex = Form("mean_2 =  %f", mean_1);
// 	tex1 = new TLatex(0.5,0.65, fit_param_latex);
// 	tex1->SetNDC();
// 	tex1->Draw();
// 
// 	mean_1 = h3->GetMean();
// 	fit_param_latex = Form("mean_3 =  %f", mean_1);
// 	tex1 = new TLatex(0.5,0.5, fit_param_latex);
// 	tex1->SetNDC();
// 	tex1->Draw();
// 	pad11->Update();
	
	
	canvas->Update();
   	canvas->cd();
   	
	TPad* pad22;
	pad22 = new TPad("pad2", "pad2", 0, 0.001, 1, 0.15);
	pad22->Draw();
	pad22->cd();
	
	ratio_2->Divide(h1);
	ratio_2->GetYaxis()->SetRangeUser(0.75, 1.25);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_2->GetYaxis()->SetTitleSize(0.05);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	TString label_Yaxis = Form("X / %s", histos_name[0].Data());
	ratio_2->GetYaxis()->SetTitle(label_Yaxis);
	ratio_2->GetXaxis()->SetTitle(x_name);
	ratio_2->Draw("E");
	ratio_3->Divide(h1);
	ratio_3->GetYaxis()->SetRangeUser(0.5, 1.05);
	ratio_3->SetTitle("");
	ratio_3->SetStats(0);
	ratio_3->GetYaxis()->SetTitleSize(0.05);
	ratio_3->GetYaxis()->SetLabelSize(0.14);
	ratio_3->GetYaxis()->SetTitle(label_Yaxis);
	ratio_3->GetXaxis()->SetTitle(x_name);
	ratio_3->Draw("same E");
	canvas->Update();

	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(1);
	line->Draw();

	canvas->Update();
   	
	canvas->Print(save);
	
}
