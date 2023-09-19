#include "Higgs_Mass_setup.h"
float dataMCErr(float pt, float eta, TH2F* hMuScaleFacUnc);

//TString USER_DIR = "/afs/cern.ch/work/f/ferrico/private/Rochester/CMSSW_12_4_3/src/";
//TString USER_DIR = "/afs/cern.ch/work/x/xzuo/UF_NanoAODTool_10619p2/src/Roch_test/";
TString USER_DIR = "/eos/user/j/jovonden/rochester/Roch_test/";
//std::cout<<USER_DIR<<std::endl;


void make_hists(std::vector<TString> samples, float lumi, float crossSection, bool Anno_2018, float reweight){
	TFile* _file0;
	TTree *tree;                                             
	Long64_t nentries;
	TRandom3 rand;  
                                                                                                                                                                                                       
	double u1;	

	std::vector<Double_t> pt_bins;
	pt_bins.push_back(5);
	pt_bins.push_back(10);
	pt_bins.push_back(15);
	pt_bins.push_back(20);
	pt_bins.push_back(25);
	pt_bins.push_back(30);
	pt_bins.push_back(35);
	pt_bins.push_back(40);
	pt_bins.push_back(45);
	pt_bins.push_back(50);
	pt_bins.push_back(55);
	pt_bins.push_back(60);
	pt_bins.push_back(65);
	pt_bins.push_back(70);
	pt_bins.push_back(75);
	pt_bins.push_back(80);
	pt_bins.push_back(100);
	pt_bins.push_back(140);
	pt_bins.push_back(200);
	
    std::vector<Double_t> eta_bins;
// 	for(int i = 0; i < 17; i++)
// 		eta_bins.push_back(-2.4 + 0.3*i);
	eta_bins.push_back(-2.4);
// eta_bins.push_back(0);
	eta_bins.push_back(2.4);


	std::vector<Double_t> phi_bins;
// 	for(int i = 0; i < 17; i++)
// 		phi_bins.push_back(-3.2 + 0.4*i);
	phi_bins.push_back(-3.2);
// phi_bins.push_back(0);
	phi_bins.push_back(3.2);
	
	std::vector<Double_t> qOverPt_bins;
	for(int i = 0; i < 21; i++)
		qOverPt_bins.push_back(0 + (float)i/200);

	std::vector<Double_t> mass_bins;
	for(int i = 0; i < 61; i++)
		mass_bins.push_back(60 + (float)i);

	TH1F* h_Zboson_mass[3];
	TH1F* h_Zboson_pT[3];
	TH1F* h_lepton_qOverPt[3];
	TH1F* h_lepton_qOverPt_scan[3][40][40];
	TH1F* h_lepton_pt[3];
	TH1F* h_lepton_eta[3];
	TH1F* h_lepton_phi[3];

	TH2F* h_Zboson_mass_vs_phi[3][2];
	TH2F* h_Zboson_mass_vs_eta[3][2];


	TH3F* h_Zboson_mass_EtaPhi[3];
	TH3F* h_lepton_qOverPt_scan_EtaPhi[3][2];
	
	TH2F* h_pTRes = new TH2F("h_pTRes","h_pTRes", pt_bins.size()-1, &pt_bins[0], 100, -0.1, 0.1);

	for(int i = 0; i < samples.size(); i++){
		TString histo_name = "h_Zboson_mass_" + samples.at(i);
		h_Zboson_mass[i] = new TH1F(histo_name, histo_name, 80, 81, 101);
		histo_name = "h_Zboson_pT_" + samples.at(i);
		h_Zboson_pT[i] = new TH1F(histo_name, histo_name, 80, 0., 160);
		histo_name = "h_lepton_qOverPt_" + samples.at(i);
		h_lepton_qOverPt[i] = new TH1F(histo_name, histo_name, 50, -0.1, 0.1);
		histo_name = "h_lepton_pt_" + samples.at(i);
		h_lepton_pt[i] = new TH1F(histo_name, histo_name, 40, 0., 200);
		histo_name = "h_lepton_eta_" + samples.at(i);
		h_lepton_eta[i] = new TH1F(histo_name, histo_name, 24, -2.4, 2.4);
		histo_name = "h_lepton_phi_" + samples.at(i);
		h_lepton_phi[i] = new TH1F(histo_name, histo_name, 1, -3.2, 3.2);
		
		histo_name = "h_Zboson_mass_EtaPhi_" + samples.at(i);
		h_Zboson_mass_EtaPhi[i] = new TH3F(histo_name, histo_name, eta_bins.size()-1, &eta_bins[0], phi_bins.size()-1, &phi_bins[0], mass_bins.size()-1, &mass_bins[0]);
		for(int q = 0; q < 2; q++){
			if(q == 0)
				histo_name = Form("h_lepton_qOverPt_scan_EtaPhi_%s_positive", samples[i].Data());
			else
				histo_name = Form("h_lepton_qOverPt_scan_EtaPhi_%s_negative", samples[i].Data());
			h_lepton_qOverPt_scan_EtaPhi[i][q] = new TH3F(histo_name, histo_name, eta_bins.size()-1, &eta_bins[0], phi_bins.size()-1, &phi_bins[0], qOverPt_bins.size()-1, &qOverPt_bins[0]);


			if(q == 0)
				histo_name = Form("h_Zboson_mass_vs_phi_%s_positive", samples[i].Data());
			else
				histo_name = Form("h_Zboson_mass_vs_phi_%s_negative", samples[i].Data());
			h_Zboson_mass_vs_phi[i][q] = new TH2F(histo_name, histo_name, phi_bins.size()-1, &phi_bins[0], 80, 81, 101);


			if(q == 0)
				histo_name = Form("h_Zboson_mass_vs_eta_%s_positive", samples[i].Data());
			else
				histo_name = Form("h_Zboson_mass_vs_eta_%s_negative", samples[i].Data());
			h_Zboson_mass_vs_eta[i][q] = new TH2F(histo_name, histo_name, eta_bins.size()-1, &eta_bins[0], 80, 81, 101);



		}

		for(int j = 0; j < 40; j++){
			for(int h = 0; h < 40; h++){
				histo_name = Form("h_lepton_qOverPt_%s_%d_%d", samples[i].Data(), j, h);
				h_lepton_qOverPt_scan[i][j][h] = new TH1F(histo_name, histo_name, 50, -0.1, 0.1);
				
			}
		}
	}
		

	std::vector<double> *muon_pt = 0; 
	std::vector<double> *muon_eta = 0; 
	std::vector<double> *muon_phi = 0; 
	std::vector<double> *muon_mass = 0; 
	std::vector<double> *muon_charge = 0; 
	std::vector<double> *muon_isLoose = 0; 
	std::vector<double> *muon_isMedium = 0; 
	std::vector<double> *muon_isTight = 0; 
	std::vector<double> *muon_trackIso = 0; 
	
	std::vector<double> *GEN_pt = 0; 
	std::vector<double> *GEN_eta = 0; 
	std::vector<double> *GEN_phi = 0; 
	std::vector<double> *GEN_mass = 0; 
	std::vector<double> *GEN_id = 0; 
	std::vector<double> *GEN_status = 0; 
	std::vector<double> *GEN_motherId = 0; 
	

	std::vector<double> *hlt_pt = 0; 
	std::vector<double> *hlt_eta = 0; 
	std::vector<double> *hlt_phi = 0; 

	bool passedTrig;

	cout << muon_pt << endl;

// 	for(int s = 0; s < samples.size(); s++){
	for(int s = 0; s < 2; s++){
        TString file1, userdir;
		userdir = "/eos/user/j/jovonden/rochester/Roch_test/";

        file1 = userdir + samples.at(s);
        file1 += ".root";
        std::cout<< "step1" <<file1<<std::endl;
		if(Anno_2018) _file0 = new TFile(file1);

		if(_file0){
			tree = (TTree*)gDirectory->Get("Ana/passedEvents");	
		}
		else{
			std::cout<<"File not found"<<std::endl;
			return;
		}

		if(!tree){
			std::cout<<"tree not found"<<std::endl;
			return;
		}

		TH1F* h_num_eventi = (TH1F*)_file0->Get("Ana/sumWeights");
		float num_event = h_num_eventi->Integral();

		tree->SetBranchAddress("muon_pt", &muon_pt);	
		tree->SetBranchAddress("muon_eta", &muon_eta);	
		tree->SetBranchAddress("muon_phi", &muon_phi);	
		tree->SetBranchAddress("muon_mass", &muon_mass);	
		tree->SetBranchAddress("muon_charge", &muon_charge);	
		tree->SetBranchAddress("muon_isLoose", &muon_isLoose);	
		tree->SetBranchAddress("muon_isMedium", &muon_isMedium);	
		tree->SetBranchAddress("muon_isTight", &muon_isTight);	
		tree->SetBranchAddress("muon_trackIso", &muon_trackIso);	

		tree->SetBranchAddress("GEN_pt", &GEN_pt);	
		tree->SetBranchAddress("GEN_eta", &GEN_eta);	
		tree->SetBranchAddress("GEN_phi", &GEN_phi);	
		tree->SetBranchAddress("GEN_mass", &GEN_mass);	
		tree->SetBranchAddress("GEN_id", &GEN_id);	
		tree->SetBranchAddress("GEN_status", &GEN_status);	
		tree->SetBranchAddress("GEN_motherId", &GEN_motherId);	

		tree->SetBranchAddress("hlt_pt", &hlt_pt);	
		tree->SetBranchAddress("hlt_eta", &hlt_eta);	
		tree->SetBranchAddress("hlt_phi", &hlt_phi);	

		tree->SetBranchAddress("passedTrig", &passedTrig);	
	
		nentries = tree->GetEntries();
		
		std::cout<<nentries<<std::endl;		
	


			
		for(int entry = 0; entry < nentries; entry++){
// 		for(int entry = 0; entry < 5000000; entry++){
// 		for(int entry = 0; entry < 1500000; entry++){
// 		for(int entry = 0; entry < 250000; entry++){
// 		for(int entry = 0; entry < 5000; entry++){
			if (entry % 1000 != 0) continue;

			tree->GetEntry(entry);
			if(entry == 0) std::cout<<weight<<std::endl;
			
			if(entry % 250000 == 0)
				std::cout<<entry<<"\t / \t"<<nentries<<std::endl;	

// 			if(!passedTrig) continue;

			if(muon_pt->size() < 2) continue;

			std::vector<int> muon_hlt_match;
			std::vector<int> muon_GEN_match;
		

			for(int j = 0; j < muon_pt->size(); j++){
				if(muon_hlt_match.size() > 1) continue;

				float deltaR = 1000;
				if(muon_pt->at(j) < 10) continue;
				if(fabs(muon_eta->at(j) > 2.4)) continue;
				if(muon_trackIso->at(j) / muon_pt->at(j) > 0.1) continue;
				if(!muon_isMedium) continue;
 				//if(!muon_isMedium) continue;
// 				if(!muon_isLoose) continue;
							
				for(int i = 0; i < hlt_pt->size(); i++){

					float deltaEta = hlt_eta->at(i) - muon_eta->at(j);
					float deltaPhi = hlt_phi->at(i) - muon_phi->at(j);
					float tmp_deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
					if(tmp_deltaR < deltaR) deltaR = tmp_deltaR;
				}
				if(deltaR < 0.35)
					muon_hlt_match.push_back(j);
			
				deltaR = 1000;
				int k = -1;
				for(int i = 0; i < GEN_pt->size(); i++){
					if(GEN_status->at(i) != 1) continue;
					if(fabs(GEN_id->at(i)) != 13) continue;
					if(GEN_motherId->at(i) != 23) continue;
										
					float deltaEta = GEN_eta->at(i) - muon_eta->at(j);
					float deltaPhi = GEN_phi->at(i) - muon_phi->at(j);
					float tmp_deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
					if(tmp_deltaR < deltaR){
						deltaR = tmp_deltaR;
						if(deltaR < 0.1)
							k = i;
					}
				}
				if(k != -1)
					muon_GEN_match.push_back(k);

			}
				
			if(muon_hlt_match.size() < 2) continue;
			if(muon_pt->at(muon_hlt_match[0]) < 25) continue;
			if(muon_charge->at(muon_hlt_match[0]) + muon_charge->at(muon_hlt_match[1]) != 0) continue;
					
			TLorentzVector lep0, lep1;
			lep0.SetPtEtaPhiM(muon_pt->at(muon_hlt_match[0]), muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]), muon_mass->at(muon_hlt_match[0]));
			lep1.SetPtEtaPhiM(muon_pt->at(muon_hlt_match[1]), muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]), muon_mass->at(muon_hlt_match[1]));
			TLorentzVector Z = lep0+lep1;
		
			if(Z.M() > 81 && Z.M() < 101){		
				if(s != 1){
					weight = lumi*crossSection/num_event;
				}
				else
					weight = 1;
				
				h_Zboson_pT[s]->Fill(Z.Pt(), weight);

				if(s == 0){
					if(muon_GEN_match.size() < 2) continue;
					if(GEN_id->at(muon_GEN_match[0]) + GEN_id->at(muon_GEN_match[1]) != 0) continue;

					TLorentzVector lep0, lep1;
					lep0.SetPtEtaPhiM(GEN_pt->at(muon_GEN_match[0]), GEN_eta->at(muon_GEN_match[0]), GEN_phi->at(muon_GEN_match[0]), GEN_mass->at(muon_GEN_match[0]));
					lep1.SetPtEtaPhiM(GEN_pt->at(muon_GEN_match[1]), GEN_eta->at(muon_GEN_match[1]), GEN_phi->at(muon_GEN_match[1]), GEN_mass->at(muon_GEN_match[1]));
					TLorentzVector Z = lep0+lep1;
	
					if(muon_charge->at(muon_hlt_match[0]) * GEN_id->at(muon_GEN_match[0]) < 0){
						h_pTRes->Fill(muon_pt->at(muon_hlt_match[0]), (GEN_pt->at(muon_GEN_match[0]) - muon_pt->at(muon_hlt_match[0])) / GEN_pt->at(muon_GEN_match[0]));
						h_pTRes->Fill(muon_pt->at(muon_hlt_match[1]), (GEN_pt->at(muon_GEN_match[1]) - muon_pt->at(muon_hlt_match[1])) / GEN_pt->at(muon_GEN_match[1]));
					}
					else{
						h_pTRes->Fill(muon_pt->at(muon_hlt_match[0]), (GEN_pt->at(muon_GEN_match[1]) - muon_pt->at(muon_hlt_match[0])) / GEN_pt->at(muon_GEN_match[1]));
						h_pTRes->Fill(muon_pt->at(muon_hlt_match[1]), (GEN_pt->at(muon_GEN_match[0]) - muon_pt->at(muon_hlt_match[1])) / GEN_pt->at(muon_GEN_match[0]));
					}				
				}
			
			}
		}
	}

	TH1F* h_resGEN = new TH1F("h_resGEN", "h_resGEN", pt_bins.size()-1, &pt_bins[0]);

	for(int xxx = 1; xxx < h_pTRes->GetNbinsX()+1; xxx++){
		TString save = Form("rochester_v1%.0d.pdf", xxx);

		TString name = Form("Res_pt_%.0f_%.0f", pt_bins.at(xxx-1), pt_bins.at(xxx));

		TH1F* h_tmp_MC	= (TH1F*)h_pTRes->ProjectionY(name, xxx, xxx);
		TCanvas* c1 = new TCanvas(name,name, 600, 600);
// 		h_tmp_MC->Draw();
		RooRealVar var("var", "var", -0.20, 0.20);
		RooRealVar MeanA("MeanA", "MeanA", 0, -0.01, 0.01);
		RooRealVar Sigma_GaussA("Sigma_GaussA", "Sigma_GaussA", 0.01, 0.0001, 0.5);
		RooGaussian GaussA("GaussA", "GaussA", var, MeanA, Sigma_GaussA);
		RooDataHist histo("prova_MC", "prova_MC", var, h_tmp_MC);
		RooPlot* xframe = var.frame(Title(name));
		histo.plotOn(xframe);
		GaussA.fitTo(histo, Range(-0.03, 0.03));
		GaussA.plotOn(xframe,RooFit::LineColor(kRed+2));
		GaussA.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));
		xframe->Draw();	
		
		h_resGEN->SetBinContent(xxx, Sigma_GaussA.getVal());
		
		c1->SaveAs(save);
	}
	TH1F* h_ZBoson_pt_correction = (TH1F*)h_Zboson_pT[1]->Clone();

}
/*
	DrawResolution(h_Zboson_pT[0], h_Zboson_pT[1], h_Zboson_pT[2], samples, "Zboson_pT", save, "pT_{2l} (GeV)");

	h_ZBoson_pt_correction->Divide(h_Zboson_pT[0]);
	TCanvas* c_ZBoson_pt_correction = new TCanvas("ZBoson_pt_correction", "ZBoson_pt_correction", 600, 600);
	h_ZBoson_pt_correction->Draw();
	c_ZBoson_pt_correction->Print(save);
	
	float ZbosonPt = 1;

	for(int s = 0; s < 3; s++){
		h_Zboson_pT[s]->Reset();
	}

	for(int s = 0; s < 2; s++){

        TString file1, userdir;
        userdir = "/eos/user/j/jovonden/rochester/Roch_test/";
        file1 = userdir + samples.at(s) + ".root";
        std::cout<< "step1" <<file1<<std::endl;
		if(Anno_2018) _file0 = new TFile(file1, "READ");

		if(_file0){
			tree = (TTree*)gDirectory->Get("Ana/passedEvents");	
		}
		else{
			std::cout<<"File not found"<<std::endl;
			return;
		}

		if(!tree){
			std::cout<<"tree not found"<<std::endl;
			return;
		}

		TH1F* h_num_eventi = (TH1F*)_file0->Get("Ana/sumWeights");
		float num_event = h_num_eventi->Integral();
		std::cout<<"Num event = "<<num_event<<std::endl;

		tree->SetBranchAddress("muon_pt", &muon_pt);	
		tree->SetBranchAddress("muon_eta", &muon_eta);	
		tree->SetBranchAddress("muon_phi", &muon_phi);	
		tree->SetBranchAddress("muon_mass", &muon_mass);	
		tree->SetBranchAddress("muon_charge", &muon_charge);	
		tree->SetBranchAddress("muon_isLoose", &muon_isLoose);	
		tree->SetBranchAddress("muon_isMedium", &muon_isMedium);	
		tree->SetBranchAddress("muon_isTight", &muon_isTight);	
		tree->SetBranchAddress("muon_trackIso", &muon_trackIso);	

		tree->SetBranchAddress("GEN_pt", &GEN_pt);	
		tree->SetBranchAddress("GEN_eta", &GEN_eta);	
		tree->SetBranchAddress("GEN_phi", &GEN_phi);	
		tree->SetBranchAddress("GEN_mass", &GEN_mass);	
		tree->SetBranchAddress("GEN_id", &GEN_id);	
		tree->SetBranchAddress("GEN_status", &GEN_status);	
		tree->SetBranchAddress("GEN_motherId", &GEN_motherId);	

		tree->SetBranchAddress("hlt_pt", &hlt_pt);	
		tree->SetBranchAddress("hlt_eta", &hlt_eta);	
		tree->SetBranchAddress("hlt_phi", &hlt_phi);	

		tree->SetBranchAddress("passedTrig", &passedTrig);	
	
		nentries = tree->GetEntries();
		
		std::cout<<nentries<<std::endl;			
			
		for(int entry = 0; entry < nentries; entry++){
// 		for(int entry = 0; entry < 2500000; entry++){
// 		for(int entry = 0; entry < 1500000; entry++){
// 		for(int entry = 0; entry < 500000; entry++){
// 		for(int entry = 0; entry < 100000; entry++){

			tree->GetEntry(entry);
			if(entry == 0) std::cout<<weight<<std::endl;
			
			if(entry % 250000 == 0)
				std::cout<<entry<<"\t / \t"<<nentries<<std::endl;	

// 			if(!passedTrig) continue;

			if(muon_pt->size() < 2) continue;

			std::vector<int> muon_hlt_match;
			std::vector<int> muon_GEN_match;
		
			for(int j = 0; j < muon_pt->size(); j++){
				if(muon_hlt_match.size() > 1) continue;

				float deltaR = 1000;
				if(muon_pt->at(j) < 10) continue;
				if(fabs(muon_eta->at(j) > 2.4)) continue;
				if(muon_trackIso->at(j) / muon_pt->at(j) > 0.1) continue;
				if(!muon_isMedium) continue;
 				//if(!muon_isMedium) continue;
// 				if(!muon_isLoose) continue;
				
			
				for(int i = 0; i < hlt_pt->size(); i++){

					float deltaEta = hlt_eta->at(i) - muon_eta->at(j);
					float deltaPhi = hlt_phi->at(i) - muon_phi->at(j);
					float tmp_deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
					if(tmp_deltaR < deltaR) deltaR = tmp_deltaR;
				}
				if(deltaR < 0.35)
					muon_hlt_match.push_back(j);
			
				if(s == 0){
					deltaR = 1000;
					int k = -1;
					for(int i = 0; i < GEN_pt->size(); i++){
						if(GEN_status->at(i) != 1) continue;
						if(fabs(GEN_id->at(i)) != 13) continue;
						if(GEN_motherId->at(i) != 23) continue;

						float deltaEta = GEN_eta->at(i) - muon_eta->at(j);
						float deltaPhi = GEN_phi->at(i) - muon_phi->at(j);
						float tmp_deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
						if(tmp_deltaR < deltaR){
							deltaR = tmp_deltaR;
							if(deltaR < 0.1)
								k = i;
						}
					}
					if(k != -1)
						muon_GEN_match.push_back(k);
				}

			}
				
			if(muon_hlt_match.size() < 2) continue;
			if(muon_charge->at(muon_hlt_match[0]) + muon_charge->at(muon_hlt_match[1]) != 0) continue;
			if(muon_pt->at(muon_hlt_match[0]) < 25) continue;

					
			TLorentzVector lep0, lep1;
			lep0.SetPtEtaPhiM(muon_pt->at(muon_hlt_match[0]), muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]), muon_mass->at(muon_hlt_match[0]));
			lep1.SetPtEtaPhiM(muon_pt->at(muon_hlt_match[1]), muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]), muon_mass->at(muon_hlt_match[1]));
			TLorentzVector Z = lep0+lep1;
		
			if(Z.M() > 81 && Z.M()<101){		
// 			if(Z.M() > 60 && Z.M()<120){		

				if(s == 0){
					ZbosonPt = h_ZBoson_pt_correction->GetBinContent(h_ZBoson_pt_correction->FindBin(Z.Pt()));
					weight = lumi*crossSection/num_event;

					reweight = ZbosonPt * weight;
				}
				else{
					weight = 1;
					reweight = 1;
				}

					
				h_Zboson_mass[s]->Fill(Z.M(), reweight);
				h_Zboson_pT[s]->Fill(Z.Pt(), reweight);

				h_lepton_qOverPt[s]->Fill(muon_charge->at(muon_hlt_match[0]) / muon_pt->at(muon_hlt_match[0]), reweight);
				h_lepton_qOverPt[s]->Fill(muon_charge->at(muon_hlt_match[1]) / muon_pt->at(muon_hlt_match[1]), reweight);

				h_lepton_pt[s]->Fill(muon_pt->at(muon_hlt_match[0]), reweight);
				h_lepton_pt[s]->Fill(muon_pt->at(muon_hlt_match[1]), reweight);

				h_lepton_eta[s]->Fill(muon_eta->at(muon_hlt_match[0]), reweight);
				h_lepton_eta[s]->Fill(muon_eta->at(muon_hlt_match[1]), reweight);

				h_lepton_phi[s]->Fill(muon_phi->at(muon_hlt_match[0]), reweight);
				h_lepton_phi[s]->Fill(muon_phi->at(muon_hlt_match[1]), reweight);


				if(muon_charge->at(muon_hlt_match[0]) > 0){
					h_lepton_qOverPt_scan_EtaPhi[s][0]->Fill(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]), (float)1 / muon_pt->at(muon_hlt_match[0]), reweight);
					h_lepton_qOverPt_scan_EtaPhi[s][1]->Fill(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]), (float)1 / muon_pt->at(muon_hlt_match[1]), reweight);

					h_Zboson_mass_vs_eta[s][0]->Fill(muon_eta->at(muon_hlt_match[0]), Z.M(), reweight);
					h_Zboson_mass_vs_eta[s][1]->Fill(muon_eta->at(muon_hlt_match[1]), Z.M(), reweight);
					h_Zboson_mass_vs_phi[s][0]->Fill(muon_phi->at(muon_hlt_match[0]), Z.M(), reweight);
					h_Zboson_mass_vs_phi[s][1]->Fill(muon_phi->at(muon_hlt_match[1]), Z.M(), reweight);


				}
				else{
					h_lepton_qOverPt_scan_EtaPhi[s][0]->Fill(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]), (float)1 / muon_pt->at(muon_hlt_match[1]), reweight);
					h_lepton_qOverPt_scan_EtaPhi[s][1]->Fill(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]), (float)1 / muon_pt->at(muon_hlt_match[0]), reweight);

					h_Zboson_mass_vs_eta[s][0]->Fill(muon_eta->at(muon_hlt_match[1]), Z.M(), reweight);
					h_Zboson_mass_vs_eta[s][1]->Fill(muon_eta->at(muon_hlt_match[0]), Z.M(), reweight);
					h_Zboson_mass_vs_phi[s][0]->Fill(muon_phi->at(muon_hlt_match[1]), Z.M(), reweight);
					h_Zboson_mass_vs_phi[s][1]->Fill(muon_phi->at(muon_hlt_match[0]), Z.M(), reweight);


				}


								
				if(s == 0){
					if(muon_GEN_match.size() < 2) continue;
					if(GEN_id->at(muon_GEN_match[0]) + GEN_id->at(muon_GEN_match[1]) != 0) continue;
					
					float GEN_pt0;
// 					GEN_pt0 = GEN_pt->at(muon_GEN_match[0]);
// 					rand.SetSeed(abs(static_cast<int>(sin(muon_phi->at(muon_hlt_match[0]))*100000)));
					GEN_pt0 = GEN_pt->at(muon_GEN_match[0]) * (1 + rand.Gaus(0, h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[0])))));
// 					GEN_pt0 = rand.Gaus(GEN_pt->at(muon_GEN_match[0]), GEN_pt->at(muon_GEN_match[0]) * h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[0]))));
// 					std::cout<<GEN_pt->at(muon_GEN_match[0])<<"\t"<<h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[0]))<<"\t"<<h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[0])))<<"\t"<<GEN_pt0<<std::endl;
					
					float GEN_pt1;
// 					GEN_pt1 = GEN_pt->at(muon_GEN_match[1]);
// 					rand.SetSeed(abs(static_cast<int>(sin(muon_phi->at(muon_hlt_match[1]))*100000)));
					GEN_pt1 = GEN_pt->at(muon_GEN_match[1]) * (1 + rand.Gaus(0, h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[1])))));
// 					GEN_pt1 = rand.Gaus(GEN_pt->at(muon_GEN_match[1]), GEN_pt->at(muon_GEN_match[1]) * h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[1]))));
// 					std::cout<<GEN_pt->at(muon_GEN_match[1])<<"\t"<<h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[1]))<<"\t"<<h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[1])))<<"\t"<<GEN_pt1<<std::endl;
					

					TLorentzVector lep0, lep1;
					lep0.SetPtEtaPhiM(GEN_pt0, GEN_eta->at(muon_GEN_match[0]), GEN_phi->at(muon_GEN_match[0]), GEN_mass->at(muon_GEN_match[0]));
					lep1.SetPtEtaPhiM(GEN_pt1, GEN_eta->at(muon_GEN_match[1]), GEN_phi->at(muon_GEN_match[1]), GEN_mass->at(muon_GEN_match[1]));
					TLorentzVector Z = lep0+lep1;
		
					h_Zboson_mass[2]->Fill(Z.M(), reweight);
					h_Zboson_pT[2]->Fill(Z.Pt(), reweight);

					h_lepton_qOverPt[2]->Fill(GEN_id->at(muon_GEN_match[0]) / ((float)13 * GEN_pt0), reweight);
					h_lepton_qOverPt[2]->Fill(GEN_id->at(muon_GEN_match[1]) / ((float)13 * GEN_pt1), reweight);

					h_lepton_pt[2]->Fill(GEN_pt0, reweight);
					h_lepton_pt[2]->Fill(GEN_pt1, reweight);

					h_lepton_eta[2]->Fill(GEN_eta->at(muon_GEN_match[0]), reweight);
					h_lepton_eta[2]->Fill(GEN_eta->at(muon_GEN_match[1]), reweight);

					h_lepton_phi[2]->Fill(GEN_phi->at(muon_GEN_match[0]), reweight);
					h_lepton_phi[2]->Fill(GEN_phi->at(muon_GEN_match[1]), reweight);
					
					if(GEN_id->at(muon_GEN_match[0]) < 0){
						h_lepton_qOverPt_scan_EtaPhi[2][0]->Fill(GEN_eta->at(muon_GEN_match[0]), GEN_phi->at(muon_GEN_match[0]), (float)1 / ( GEN_pt0), weight);
						h_lepton_qOverPt_scan_EtaPhi[2][1]->Fill(GEN_eta->at(muon_GEN_match[1]), GEN_phi->at(muon_GEN_match[1]), (float)1 / ( GEN_pt1), weight);

						h_Zboson_mass_vs_eta[2][0]->Fill(muon_eta->at(muon_hlt_match[0]), Z.M(), reweight);
						h_Zboson_mass_vs_eta[2][1]->Fill(muon_eta->at(muon_hlt_match[1]), Z.M(), reweight);
						h_Zboson_mass_vs_phi[2][0]->Fill(muon_phi->at(muon_hlt_match[0]), Z.M(), reweight);
						h_Zboson_mass_vs_phi[2][1]->Fill(muon_phi->at(muon_hlt_match[1]), Z.M(), reweight);
					}
					else{
						h_lepton_qOverPt_scan_EtaPhi[2][0]->Fill(GEN_eta->at(muon_GEN_match[1]), GEN_phi->at(muon_GEN_match[1]), (float)1 / ( GEN_pt1), weight);
						h_lepton_qOverPt_scan_EtaPhi[2][1]->Fill(GEN_eta->at(muon_GEN_match[0]), GEN_phi->at(muon_GEN_match[0]), (float)1 / ( GEN_pt0), weight);

						h_Zboson_mass_vs_eta[2][0]->Fill(muon_eta->at(muon_hlt_match[1]), Z.M(), reweight);
						h_Zboson_mass_vs_eta[2][1]->Fill(muon_eta->at(muon_hlt_match[0]), Z.M(), reweight);
						h_Zboson_mass_vs_phi[2][0]->Fill(muon_phi->at(muon_hlt_match[1]), Z.M(), reweight);
						h_Zboson_mass_vs_phi[2][1]->Fill(muon_phi->at(muon_hlt_match[0]), Z.M(), reweight);
					}

				}
				
			}

//	TFile* _fileH = new TFile(file1, "READ");
	
	//TH1F* h_num_eventi = (TH1F*)_file0->Get("Ana/sumWeights");

// 	DrawResolution(h_Zboson_mass[0], h_Zboson_mass[1], h_Zboson_mass[2], samples, "Zboson_mass", save + "[", "mass_{2l} (GeV)");
	DrawResolution(h_Zboson_mass[0], h_Zboson_mass[1], h_Zboson_mass[2], samples, "Zboson_mass", save, "mass_{2l} (GeV)");
	DrawResolution(h_Zboson_pT[0], h_Zboson_pT[1], h_Zboson_pT[2], samples, "Zboson_pT", save, "pT_{2l} (GeV)");
	DrawResolution(h_lepton_pt[0], h_lepton_pt[1], h_lepton_pt[2], samples, "h_lepton_pt", save, "p_{T} (GeV)");
	DrawResolution(h_lepton_eta[0], h_lepton_eta[1], h_lepton_eta[2], samples, "h_lepton_eta", save, "#eta");
	DrawResolution(h_lepton_phi[0], h_lepton_phi[1], h_lepton_phi[2], samples, "h_lepton_phi", save, "#phi");
	DrawResolution(h_lepton_qOverPt[0], h_lepton_qOverPt[1], h_lepton_qOverPt[2], samples, "h_lepton_qOverPt", save, "q/p_{T} (GeV)");
// 	DrawResolution(h_lepton_qOverPt[0], h_lepton_qOverPt[1], h_lepton_qOverPt[2], samples, "h_lepton_qOverPt", save + "]", "q/p_{T} (GeV)");
//     return;

}
	}}

*/

void rochester(){
	gROOT->Reset();
    gROOT->SetBatch();	
    gStyle->SetOptStat(0);
    
    bool Anno_2018 = 1;
    float lumi = 86.6855;
    if(Anno_2018) lumi = 373;
    
    float reweight = 1;   

	float crossSection = 13730; //pb (from DAS)
	float weight = 1;
	
	std::vector<TString> samples;
	samples.push_back("MC");	
	samples.push_back("DATA");	
	samples.push_back("GEN");	

	make_hists(samples, lumi, crossSection, Anno_2018, reweight);





}


/*


	TMultiGraph* Mass_vs_phi_pos = new TMultiGraph();
	TMultiGraph* Mass_vs_phi_neg = new TMultiGraph();
	TMultiGraph* Mass_vs_eta_pos = new TMultiGraph();
	TMultiGraph* Mass_vs_eta_neg = new TMultiGraph();

	for(int s = 0; s < 3; s++){
		std::vector<Double_t> mean_pos;
		std::vector<Double_t> mean_neg;
		std::vector<Double_t> eta_bins_centered;
		std::vector<Double_t> phi_bins_centered;
		for(int xxx = 1; xxx < h_Zboson_mass_vs_eta[s][0]->GetNbinsX()+1; xxx++){
			TString name = Form("Mass_eta_%.1f_%.1f_pos", eta_bins.at(xxx-1), eta_bins.at(xxx));
// 			std::cout<<"ETA_name = "<<name<<std::endl;
			
			TH1F* h_tmp_MC	= (TH1F*)h_Zboson_mass_vs_eta[s][0]->ProjectionY(name, xxx, xxx);		
    	    RooRealVar massZ_P("massZ","massZ",81,101);
			RooDataHist histo("mass Z", "mass Z", massZ_P, h_tmp_MC);
    
    	    RooRealVar PDGmZ("PDGmZ","PDGmZ", 91.19, 86, 96);
	        RooRealVar PDGwZ("PDGwZ","PDGwZ", 2.5, 1, 4);
	        RooBreitWigner PDGBW("PDGBW","PDGBW",massZ_P,PDGmZ,PDGwZ);
	        RooBreitWigner PDGBW_2("PDGBW_2","PDGBW_2",massZ_P,PDGmZ,PDGwZ);
        
			RooRealVar CB_mean("CB_mean", "CB_mean", 0, -5, 5);//125, 120, 130);
			RooRealVar Sigma("Sigma", "Sigma", 1, 0.1, 20);//1, 0, 30);//sigma[decay]);
			RooRealVar AlphaL("AlphaL", "AlphaL", 1, 0.1, 30);//alphaL[decay]);
			RooRealVar ExpL("ExpL", "ExpL", 1, 0.1, 30);//expL[decay]);
			RooRealVar AlphaR("AlphaR", "AlphaR", 1, 0.1, 30);//alphaR[decay]);
			RooRealVar ExpR("ExpR", "ExpR", 1, 0.1, 50);//expR[decay]);
			RooCrystalBall DSCB("DSCB", "DSCB", massZ_P, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    
	        RooFFTConvPdf model("CW","CW",massZ_P,PDGBW,DSCB);
	        RooFFTConvPdf model_2("CW_2","CW_2",massZ_P,PDGBW,DSCB);
	        
			if(s != 20){
				model.fitTo(histo, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			else{
				PDGBW.fitTo(histo, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			mean_pos.push_back(PDGmZ.getVal());
			eta_bins_centered.push_back(eta_bins.at(xxx-1) + fabs((eta_bins.at(1) - eta_bins.at(0))/2));


			TCanvas *c_MC = new TCanvas(name, name, 900, 700);
			c_MC->SetFrameFillColor(0);
			RooPlot* xframe = massZ_P.frame("massZ");
			xframe = massZ_P.frame(Title(name));
			histo.plotOn(xframe);
			if(s != 20){
				model.plotOn(xframe,RooFit::LineColor(kBlue));
				model.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));
			}
			else{
				PDGBW.plotOn(xframe,RooFit::LineColor(kRed));
				PDGBW.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));				
			}
			xframe->Draw();
			c_MC->Print(save);
			
			
			name = Form("Mass_eta_%.1f_%.1f_neg", eta_bins.at(xxx-1), eta_bins.at(xxx));
			TH1F* h_tmp_MC_2	= (TH1F*)h_Zboson_mass_vs_eta[s][1]->ProjectionY(name, xxx, xxx);
			RooDataHist histo_neg("mass Z", "mass Z", massZ_P, h_tmp_MC_2);
			if(s != 20){
				model_2.fitTo(histo_neg, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			else{
				PDGBW_2.fitTo(histo_neg, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			mean_neg.push_back(PDGmZ.getVal());
			
			TCanvas* c_MC_2 = new TCanvas(name, name, 900, 700);
			c_MC_2->SetFrameFillColor(0);
			xframe = massZ_P.frame("massZ");
			xframe = massZ_P.frame(Title(name));
			histo_neg.plotOn(xframe);
			if(s != 20){
				model_2.plotOn(xframe,RooFit::LineColor(kBlue));
				model_2.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));
			}
			else{
				PDGBW_2.plotOn(xframe,RooFit::LineColor(kRed));
				PDGBW_2.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));				
			}
			xframe->Draw();
			c_MC_2->Print(save);		
		}
		
		TGraph* gr_new = new TGraph(eta_bins.size()-1, &eta_bins_centered[0], &mean_pos[0]);
		if(s == 1)
			gr_new->SetMarkerColor(kRed);
		if(s == 2)
			gr_new->SetMarkerColor(kGreen);
		Mass_vs_eta_pos->Add(gr_new);
		
		gr_new = new TGraph(eta_bins.size()-1, &eta_bins_centered[0], &mean_neg[0]);
		if(s == 1)
			gr_new->SetMarkerColor(kRed);
		if(s == 2)
			gr_new->SetMarkerColor(kGreen);
		Mass_vs_eta_neg->Add(gr_new);
		
// 		for(int i = 0; i < mean_pos.size(); i++)
// 			std::cout<<mean_pos.at(i)<<"\t"<<mean_neg.at(i)<<std::endl;
// 
// 		return;
		mean_pos.clear();
		mean_neg.clear();	
					
		for(int xxx = 1; xxx < h_Zboson_mass_vs_phi[s][0]->GetNbinsX()+1; xxx++){
			TString name = Form("Mass_phi_%.1f_%.1f_pos", phi_bins.at(xxx-1), phi_bins.at(xxx));
// 			std::cout<<"PHI_name = "<<name<<std::endl;
			TH1F* h_tmp_MC	= (TH1F*)h_Zboson_mass_vs_phi[s][0]->ProjectionY(name, xxx, xxx);		
    	    RooRealVar massZ_P("massZ","massZ",81,101);
			RooDataHist histo("mass Z", "mass Z", massZ_P, h_tmp_MC);
    
    	    RooRealVar PDGmZ("PDGmZ","PDGmZ",91.19, 86, 96);
	        RooRealVar PDGwZ("PDGwZ","PDGwZ", 2.5, 1, 4);
	        RooBreitWigner PDGBW("PDGBW","PDGBW",massZ_P,PDGmZ,PDGwZ);
 	        RooBreitWigner PDGBW_2("PDGBW_2","PDGBW_2",massZ_P,PDGmZ,PDGwZ);
       
			RooRealVar CB_mean("CB_mean", "CB_mean", 0, -5, 5);//125, 120, 130);
			RooRealVar Sigma("Sigma", "Sigma", 1, 0.1, 20);//1, 0, 30);//sigma[decay]);
			RooRealVar AlphaL("AlphaL", "AlphaL", 1, 0.1, 30);//alphaL[decay]);
			RooRealVar ExpL("ExpL", "ExpL", 1, 0.1, 30);//expL[decay]);
			RooRealVar AlphaR("AlphaR", "AlphaR", 1, 0.1, 30);//alphaR[decay]);
			RooRealVar ExpR("ExpR", "ExpR", 1, 0.1, 50);//expR[decay]);
			RooCrystalBall DSCB("DSCB", "DSCB", massZ_P, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    
	        RooFFTConvPdf model("CW","CW",massZ_P,PDGBW,DSCB);
	        RooFFTConvPdf model_2("CW_2","CW_2",massZ_P,PDGBW,DSCB);
	        
			if(s != 20){
				model.fitTo(histo, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			else{
				PDGBW.fitTo(histo, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			mean_pos.push_back(PDGmZ.getVal());
			phi_bins_centered.push_back(phi_bins.at(xxx-1) + fabs((phi_bins.at(1) - phi_bins.at(0))/2));

			TCanvas *c_MC = new TCanvas(name, name, 900, 700);
			c_MC->SetFrameFillColor(0);
			RooPlot* xframe = massZ_P.frame("massZ");
			xframe = massZ_P.frame(Title(name));
			histo.plotOn(xframe);
			if(s != 20){
				model.plotOn(xframe,RooFit::LineColor(kBlue));
				model.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));
			}
			else{
				PDGBW.plotOn(xframe,RooFit::LineColor(kRed));
				PDGBW.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));				
			}
			xframe->Draw();
			c_MC->Print(save);
					
			name = Form("Mass_phi_%.1f_%.1f_neg", phi_bins.at(xxx-1), phi_bins.at(xxx));
			TH1F* h_tmp_MC_2	= (TH1F*)h_Zboson_mass_vs_phi[s][1]->ProjectionY(name, xxx, xxx);
			RooDataHist histo_neg("mass Z", "mass Z", massZ_P, h_tmp_MC_2);
			if(s != 20){
				model_2.fitTo(histo_neg, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			else{
				PDGBW_2.fitTo(histo_neg, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			mean_neg.push_back(PDGmZ.getVal());

			c_MC = new TCanvas(name, name, 900, 700);
			c_MC->SetFrameFillColor(0);
			xframe = massZ_P.frame("massZ");
			xframe = massZ_P.frame(Title(name));
			histo_neg.plotOn(xframe);
			if(s != 20){
				model_2.plotOn(xframe,RooFit::LineColor(kBlue));
				model_2.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));
			}
			else{
				PDGBW_2.plotOn(xframe,RooFit::LineColor(kRed));
				PDGBW_2.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));				
			}
			xframe->Draw();
			c_MC->Print(save);		
		}
		
		gr_new = new TGraph(phi_bins.size()-1, &phi_bins_centered[0], &mean_pos[0]);
		if(s == 1)
			gr_new->SetMarkerColor(kRed);
		if(s == 2)
			gr_new->SetMarkerColor(kGreen);
		Mass_vs_phi_pos->Add(gr_new);
		
		gr_new = new TGraph(phi_bins.size()-1, &phi_bins_centered[0], &mean_neg[0]);
		if(s == 1)
			gr_new->SetMarkerColor(kRed);
		if(s == 2)
			gr_new->SetMarkerColor(kGreen);
		Mass_vs_phi_neg->Add(gr_new);
			
	}

	for(int xxx = 1; xxx < h_Zboson_mass_vs_eta[0][0]->GetNbinsX()+1; xxx++){
		TString name = Form("Mass_eta_%.1f_%.1f_pos_MC", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_MC	= (TH1F*)h_Zboson_mass_vs_eta[0][0]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_eta_%.1f_%.1f_pos_DATA", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_DATA	= (TH1F*)h_Zboson_mass_vs_eta[1][0]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_eta_%.1f_%.1f_pos", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_GEN	= (TH1F*)h_Zboson_mass_vs_eta[2][0]->ProjectionY(name, xxx, xxx);		

		DrawResolution(h_tmp_MC, h_tmp_DATA, h_tmp_GEN, samples, name, save, "mass_{2l} (GeV)");

		name = Form("Mass_eta_%.1f_%.1f_neg_MC", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_MC_2	= (TH1F*)h_Zboson_mass_vs_eta[0][1]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_eta_%.1f_%.1f_neg_DATA", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_DATA_2	= (TH1F*)h_Zboson_mass_vs_eta[1][1]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_eta_%.1f_%.1f_neg", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_GEN_2	= (TH1F*)h_Zboson_mass_vs_eta[2][1]->ProjectionY(name, xxx, xxx);	
		name = Form("Mass_eta_%.1f_%.1f_neg", eta_bins.at(xxx-1), eta_bins.at(xxx));

		DrawResolution(h_tmp_MC_2, h_tmp_DATA_2, h_tmp_GEN_2, samples, name, save, "mass_{2l} (GeV)");
	}					
	for(int xxx = 1; xxx < h_Zboson_mass_vs_phi[0][0]->GetNbinsX()+1; xxx++){
		TString name = Form("Mass_phi_%.1f_%.1f_pos_MC", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_MC	= (TH1F*)h_Zboson_mass_vs_phi[0][0]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_phi_%.1f_%.1f_pos_DATA", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_DATA	= (TH1F*)h_Zboson_mass_vs_phi[1][0]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_phi_%.1f_%.1f_pos", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_GEN	= (TH1F*)h_Zboson_mass_vs_phi[2][0]->ProjectionY(name, xxx, xxx);		

		DrawResolution(h_tmp_MC, h_tmp_DATA, h_tmp_GEN, samples, name, save, "mass_{2l} (GeV)");

		name = Form("Mass_phi_%.1f_%.1f_neg_MC", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_MC_2	= (TH1F*)h_Zboson_mass_vs_phi[0][1]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_phi_%.1f_%.1f_neg_DATA", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_DATA_2	= (TH1F*)h_Zboson_mass_vs_phi[1][1]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_phi_%.1f_%.1f_neg", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_GEN_2	= (TH1F*)h_Zboson_mass_vs_phi[2][1]->ProjectionY(name, xxx, xxx);	
		name = Form("Mass_phi_%.1f_%.1f_neg", phi_bins.at(xxx-1), phi_bins.at(xxx));

		DrawResolution(h_tmp_MC_2, h_tmp_DATA_2, h_tmp_GEN_2, samples, name, save, "mass_{2l} (GeV)");
		
	}

	TCanvas* c_Mass_vs_eta_pos = new TCanvas("c_Mass_vs_eta_pos", "c_Mass_vs_eta_pos", 600, 600);
	Mass_vs_eta_pos->SetTitle("c_Mass_vs_eta_pos");
	Mass_vs_eta_pos->Draw("AP*");
	c_Mass_vs_eta_pos->Print(save);
	TCanvas* c_Mass_vs_eta_neg = new TCanvas("c_Mass_vs_eta_neg", "c_Mass_vs_eta_neg", 600, 600);
	Mass_vs_eta_neg->SetTitle("c_Mass_vs_eta_neg");
	Mass_vs_eta_neg->Draw("AP*");
	c_Mass_vs_eta_neg->Print(save);
	TCanvas* c_Mass_vs_phi_pos = new TCanvas("c_Mass_vs_phi_pos", "c_Mass_vs_phi_pos", 600, 600);
	Mass_vs_phi_pos->SetTitle("c_Mass_vs_phi_pos");
	Mass_vs_phi_pos->Draw("AP*");
	c_Mass_vs_phi_pos->Print(save);
	TCanvas* c_Mass_vs_phi_neg = new TCanvas("c_Mass_vs_phi_neg", "c_Mass_vs_phi_neg", 600, 600);
	Mass_vs_phi_neg->SetTitle("c_Mass_vs_phi_neg");
	Mass_vs_phi_neg->Draw("AP*");
	c_Mass_vs_phi_neg->Print(save);
// 	c_Mass_vs_phi_neg->Print(save + "]");
// 
// return;

	for(int s = 0; s < 3; s++){
		h_Zboson_mass[s]->Reset();
		h_Zboson_pT[s]->Reset();
		h_lepton_pt[s]->Reset();
		h_lepton_eta[s]->Reset();
		h_lepton_phi[s]->Reset();
		h_lepton_qOverPt[s]->Reset();	
		h_Zboson_mass_vs_phi[s][0]->Reset();	
		h_Zboson_mass_vs_phi[s][1]->Reset();	
		h_Zboson_mass_vs_eta[s][0]->Reset();	
		h_Zboson_mass_vs_eta[s][1]->Reset();	
	}

	float Corr_DATA_pos;
	float Corr_MC_pos;
	float Corr_DATA_neg;
	float Corr_MC_neg;
	float D_m_MC;
	float D_a_MC;
	float D_m_DATA;
	float D_a_DATA;
	TH2F* h_Multiplicative_MC = new TH2F("Multiplicative_MC","Multiplicative_MC", eta_bins.size()-1, &eta_bins[0], phi_bins.size()-1, &phi_bins[0]);
	TH2F* h_Additive_MC = new TH2F("Additive_MC","Additive_MC", eta_bins.size()-1, &eta_bins[0], phi_bins.size()-1, &phi_bins[0]);
	TH2F* h_Multiplicative_DATA = new TH2F("Multiplicative_DATA","Multiplicative_DATA", eta_bins.size()-1, &eta_bins[0], phi_bins.size()-1, &phi_bins[0]);
	TH2F* h_Additive_DATA = new TH2F("Additive_DATA","Additive_DATA", eta_bins.size()-1, &eta_bins[0], phi_bins.size()-1, &phi_bins[0]);
	TH2F* h_denominator_MC_pos = new TH2F("h_denominator_MC_pos","h_denominator_MC_pos", eta_bins.size()-1, &eta_bins[0], phi_bins.size()-1, &phi_bins[0]);
	TH2F* h_denominator_MC_neg = new TH2F("h_denominator_MC_neg","h_denominator_MC_neg", eta_bins.size()-1, &eta_bins[0], phi_bins.size()-1, &phi_bins[0]);
	TH2F* h_denominator_DATA_pos = new TH2F("h_denominator_DATA_pos","h_denominator_DATA_pos", eta_bins.size()-1, &eta_bins[0], phi_bins.size()-1, &phi_bins[0]);
	TH2F* h_denominator_DATA_neg = new TH2F("h_denominator_DATA_neg","h_denominator_DATA_neg", eta_bins.size()-1, &eta_bins[0], phi_bins.size()-1, &phi_bins[0]);
	

	for(int xxx = 1; xxx < h_lepton_qOverPt_scan_EtaPhi[0][0]->GetNbinsX()+1; xxx++){
		for(int yyy = 1; yyy < h_lepton_qOverPt_scan_EtaPhi[0][0]->GetNbinsY()+1; yyy++){
			TString name = Form("1overPt_eta_%.1f_%.1f_phi_%.1f_%.1f", eta_bins.at(xxx-1), eta_bins.at(xxx), phi_bins.at(yyy-1), phi_bins.at(yyy));
			std::cout<<"========"<<std::endl;
			std::cout<<name<<std::endl;
			TH1F* h_tmp_MC	= (TH1F*)h_lepton_qOverPt_scan_EtaPhi[0][0]->ProjectionZ(name + "_MC+", xxx, xxx, yyy, yyy);
			TH1F* h_MC = (TH1F*)h_tmp_MC->Clone();
			TH1F* h_tmp_DATA = (TH1F*)h_lepton_qOverPt_scan_EtaPhi[1][0]->ProjectionZ(name + "_DATA+", xxx, xxx, yyy, yyy);
			TH1F* h_DATA = (TH1F*)h_tmp_DATA->Clone();
			TH1F* h_tmp_GEN	= (TH1F*)h_lepton_qOverPt_scan_EtaPhi[2][0]->ProjectionZ(name + "_GEN+", xxx, xxx, yyy, yyy);
			TH1F* h_GEN = (TH1F*)h_tmp_GEN->Clone();
			h_denominator_MC_pos->SetBinContent(xxx, yyy, h_MC->GetMean());
			h_denominator_DATA_pos->SetBinContent(xxx, yyy, h_DATA->GetMean());

			std::cout<<"Before\nDenominator+ ="<<h_MC->GetMean()<<"\t"<<h_DATA->GetMean()<<std::endl;

			DrawResolution(h_MC, h_DATA, h_GEN, samples, name + " pos", save, "1 / p_{T} (GeV)");

			Corr_MC_pos = h_GEN->GetMean() - h_MC->GetMean();
			Corr_DATA_pos = h_GEN->GetMean() - h_DATA->GetMean();
			std::cout<<"Mean+  = "<<h_GEN->GetMean()<<"\t"<<h_MC->GetMean()<<"\t"<<h_DATA->GetMean()<<std::endl;

			h_tmp_MC->Reset();
			h_tmp_DATA->Reset();
			h_tmp_GEN->Reset();
			h_MC->Reset();
			h_DATA->Reset();
			h_GEN->Reset();

			h_tmp_MC	= (TH1F*)h_lepton_qOverPt_scan_EtaPhi[0][1]->ProjectionZ(name + "_MC-", xxx, xxx, yyy, yyy);
			h_MC = (TH1F*)h_tmp_MC->Clone();
			h_tmp_DATA	= (TH1F*)h_lepton_qOverPt_scan_EtaPhi[1][1]->ProjectionZ(name + "_DATA-", xxx, xxx, yyy, yyy);
			h_DATA = (TH1F*)h_tmp_DATA->Clone();
			h_tmp_GEN	= (TH1F*)h_lepton_qOverPt_scan_EtaPhi[2][1]->ProjectionZ(name + "_GEN-", xxx, xxx, yyy, yyy);
			h_GEN = (TH1F*)h_tmp_GEN->Clone();
			h_denominator_MC_neg->SetBinContent(xxx, yyy, h_MC->GetMean());
			h_denominator_DATA_neg->SetBinContent(xxx, yyy, h_DATA->GetMean());
			std::cout<<"Denominator- ="<<h_MC->GetMean()<<"\t"<<h_DATA->GetMean()<<std::endl;
	
			DrawResolution(h_MC, h_DATA, h_GEN, samples, name + " neg", save, "1 / p_{T} (GeV)");
		
			Corr_MC_neg = h_GEN->GetMean() - h_MC->GetMean();
			Corr_DATA_neg = h_GEN->GetMean() - h_DATA->GetMean();
			std::cout<<"Mean-  = "<<h_GEN->GetMean()<<"\t"<<h_MC->GetMean()<<"\t"<<h_DATA->GetMean()<<std::endl;
// 			
			std::cout<<"C_MC = "<<Corr_MC_pos<<"\t"<<Corr_MC_neg<<std::endl;
			std::cout<<"C_DATA = "<<Corr_DATA_pos<<"\t"<<Corr_DATA_neg<<std::endl;

			D_m_MC = (Corr_MC_pos + Corr_MC_neg) / (float)2;
			D_a_MC = (Corr_MC_pos - Corr_MC_neg) / (float)2;
			D_m_DATA = (Corr_DATA_pos + Corr_DATA_neg) / (float)2;
			D_a_DATA = (Corr_DATA_pos - Corr_DATA_neg) / (float)2;
			std::cout<<"D_m_MC = "<<D_m_MC<<std::endl;
			std::cout<<"D_m_DATA = "<<D_m_DATA<<std::endl;
			std::cout<<"D_a_MC = "<<D_a_MC<<std::endl;
			std::cout<<"D_a_DATA = "<<D_a_DATA<<std::endl;

			float M = D_m_MC;
// 			float M = 1 + 2 * D_m_MC / (h_denominator_MC_pos->GetBinContent(xxx, yyy) + h_denominator_MC_neg->GetBinContent(xxx, yyy));
			float A = D_a_MC;
			h_Additive_MC->SetBinContent(xxx, yyy, A);
			h_Multiplicative_MC->SetBinContent(xxx, yyy, M);
			std::cout<<"MC: M = "<<M<<"\t A = "<<A<<std::endl;
// 			std::cout<<"Controllo MC = "<<A<<"\t"<<D_a_MC<<"\t"<<D_m_MC<<"\t"<<numerator_MC<<"\t"<<denominator_MC<<std::endl;
// 			std::cout<<"Controllo MC = "<<A<<"\t"<<M<<"\t"<<D_a_MC<<"\t"<<D_m_MC<<"\t"<<denominator_MC<<std::endl;

			M = D_m_DATA;
//             M = 1 + 2 * D_m_DATA / (h_denominator_DATA_pos->GetBinContent(xxx, yyy) + h_denominator_DATA_neg->GetBinContent(xxx, yyy));
			A = D_a_DATA;
			h_Additive_DATA->SetBinContent(xxx, yyy, A);
			h_Multiplicative_DATA->SetBinContent(xxx, yyy, M);
			std::cout<<"DATA: M = "<<M<<"\t A = "<<A<<std::endl;
// 			std::cout<<"Controllo DATA = "<<A<<"\t"<<D_a_DATA<<"\t"<<D_m_DATA<<"\t"<<numerator_DATA<<"\t"<<denominator_DATA<<std::endl;
// 			std::cout<<"Controllo DATA = "<<A<<"\t"<<M<<"\t"<<D_a_DATA<<"\t"<<D_m_DATA<<"\t"<<denominator_DATA<<std::endl;
		}
	}			
// 	return;		

	TCanvas* c_Additive_MC = new TCanvas("c_Additive_MC", "c_Additive_MC", 600, 600);
	c_Additive_MC->SetRightMargin(0.2);
	h_Additive_MC->Draw("COLZ");
	c_Additive_MC->Print(save);
	TCanvas* c_Additive_DATA = new TCanvas("c_Additive_DATA", "c_Additive_DATA", 600, 600);
	c_Additive_DATA->SetRightMargin(0.2);
	h_Additive_DATA->Draw("COLZ");
	c_Additive_DATA->Print(save);
	TCanvas* c_Multiplicative_MC = new TCanvas("c_Multiplicative_MC", "c_Multiplicative_MC", 600, 600);
	c_Multiplicative_MC->SetRightMargin(0.2);
	h_Multiplicative_MC->Draw("COLZ");
	c_Multiplicative_MC->Print(save);
	TCanvas* c_Multiplicative_DATA = new TCanvas("c_Multiplicative_DATA", "c_Multiplicative_DATA", 600, 600);
	c_Multiplicative_DATA->SetRightMargin(0.2);
	h_Multiplicative_DATA->Draw("COLZ");
	c_Multiplicative_DATA->Print(save);

// 	c_Multiplicative_DATA->Print(save + "]");
// 	return;		


		
	for(int s = 0; s < 3; s++){
		for(int q = 0; q < 2; q++){
			h_lepton_qOverPt_scan_EtaPhi[s][q]->Reset();
		}
	}


	TH2F* h_MC_mass = new TH2F("h_MC_mass","h_MC_mass", 40, 81, 101, 40, 81, 101);
	TH2F* h_DATA_mass = new TH2F("h_DATA_mass","h_DATA_mass", 40, 81, 101, 40, 81, 101);

	TH1F* h_pTVariance_MC_pos = new TH1F("h_pTVariance_MC_pos","h_pTVariance_MC_pos", 100, -0.025, 0.025);
	TH1F* h_pTVariance_MC_neg = new TH1F("h_pTVariance_MC_neg","h_pTVariance_MC_neg", 100, -0.025, 0.025);
	TH1F* h_pTVariance_DATA_pos = new TH1F("h_pTVariance_DATA_pos","h_pTVariance_DATA_pos", 100, -0.1, 0.1);
	TH1F* h_pTVariance_DATA_neg = new TH1F("h_pTVariance_DATA_neg","h_pTVariance_DATA_neg", 100, -0.1, 0.1);


// 	for(int s = 0; s < samples.size(); s++){
	for(int s = 0; s < 2; s++){

        TString file1;
        file1 = "/eos/user/j/jovonden/rochester/Roch_test/";
        file1 += samples.at(s);
        file1 += ".root";
        std::cout<< "step1" <<file1<<std::endl;
		if(Anno_2018) _file0 = new TFile(file1);


		if(_file0){
			tree = (TTree*)gDirectory->Get("Ana/passedEvents");	
		}
		else{
			std::cout<<"File not found"<<std::endl;
			return;
		}

		if(!tree){
			std::cout<<"tree not found"<<std::endl;
			return;
		}

		TH1F* h_num_eventi = (TH1F*)_file0->Get("Ana/sumWeights");
		float num_event = h_num_eventi->Integral();

		tree->SetBranchAddress("muon_pt", &muon_pt);	
		tree->SetBranchAddress("muon_eta", &muon_eta);	
		tree->SetBranchAddress("muon_phi", &muon_phi);	
		tree->SetBranchAddress("muon_mass", &muon_mass);	
		tree->SetBranchAddress("muon_charge", &muon_charge);	
		tree->SetBranchAddress("muon_isLoose", &muon_isLoose);	
		tree->SetBranchAddress("muon_isMedium", &muon_isMedium);	
		tree->SetBranchAddress("muon_isTight", &muon_isTight);	
		tree->SetBranchAddress("muon_trackIso", &muon_trackIso);	

		tree->SetBranchAddress("GEN_pt", &GEN_pt);	
		tree->SetBranchAddress("GEN_eta", &GEN_eta);	
		tree->SetBranchAddress("GEN_phi", &GEN_phi);	
		tree->SetBranchAddress("GEN_mass", &GEN_mass);	
		tree->SetBranchAddress("GEN_id", &GEN_id);	
		tree->SetBranchAddress("GEN_status", &GEN_status);	
		tree->SetBranchAddress("GEN_motherId", &GEN_motherId);	

		tree->SetBranchAddress("hlt_pt", &hlt_pt);	
		tree->SetBranchAddress("hlt_eta", &hlt_eta);	
		tree->SetBranchAddress("hlt_phi", &hlt_phi);	

		tree->SetBranchAddress("passedTrig", &passedTrig);	
	
		nentries = tree->GetEntries();
		
		std::cout<<nentries<<std::endl;		
	for(int entry = 0; entry < nentries; entry++){
// 		for(int entry = 0; entry < 2500000; entry++){
// 		for(int entry = 0; entry < 500000; entry++){
// 		for(int entry = 0; entry < 250000; entry++){
// 		for(int entry = 0; entry < 50000; entry++){

			tree->GetEntry(entry);
			if(entry == 0) std::cout<<weight<<std::endl;
	
			if(entry % 250000 == 0)
				std::cout<<entry<<"\t / \t"<<nentries<<std::endl;	

	// 			if(!passedTrig) continue;

			if(muon_pt->size() < 2) continue;

			std::vector<int> muon_hlt_match;
			std::vector<int> muon_GEN_match;

			for(int j = 0; j < muon_pt->size(); j++){
				if(muon_hlt_match.size() > 1) continue;

				float deltaR = 1000;
				if(muon_pt->at(j) < 10) continue;
				if(fabs(muon_eta->at(j) > 2.4)) continue;
				if(muon_trackIso->at(j) / muon_pt->at(j) > 0.1) continue;
				if(!muon_isMedium) continue;
	 				//if(!muon_isMedium) continue;
	// 				if(!muon_isLoose) continue;
		
	
				for(int i = 0; i < hlt_pt->size(); i++){
					if(hlt_pt->at(i) < 10) continue;

					float deltaEta = hlt_eta->at(i) - muon_eta->at(j);
					float deltaPhi = hlt_phi->at(i) - muon_phi->at(j);
					float tmp_deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
					if(tmp_deltaR < deltaR) deltaR = tmp_deltaR;
				}
				if(deltaR < 0.35)
					muon_hlt_match.push_back(j);
	
				if(s == 0){
					deltaR = 1000;
					int k = -1;
					for(int i = 0; i < GEN_pt->size(); i++){
						if(GEN_status->at(i) != 1) continue;
						if(fabs(GEN_id->at(i)) != 13) continue;
						if(GEN_motherId->at(i) != 23) continue;

						float deltaEta = GEN_eta->at(i) - muon_eta->at(j);
						float deltaPhi = GEN_phi->at(i) - muon_phi->at(j);
						float tmp_deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
						if(tmp_deltaR < deltaR){
							deltaR = tmp_deltaR;
							if(deltaR < 0.1)
								k = i;
						}
					}
					if(k != -1)
						muon_GEN_match.push_back(k);
				}

			}
		
			if(muon_hlt_match.size() < 2) continue;
			if(muon_charge->at(muon_hlt_match[0]) + muon_charge->at(muon_hlt_match[1]) != 0) continue;
			if(muon_pt->at(muon_hlt_match[0]) < 25) continue;

			TLorentzVector lep0, lep1, Z;
			float new_pt0;
			float new_pt1;
			if(muon_charge->at(muon_hlt_match[0]) > 0){
				if(s == 0){
					float Additive0 = h_Additive_MC->GetBinContent(h_Additive_MC->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));	
					float Additive1 = h_Additive_MC->GetBinContent(h_Additive_MC->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));	
// 					float Multiplicative0 = h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));
// 					float Multiplicative1 = h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));
					float Multiplicative0 = 1 + h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]))) / h_denominator_MC_pos->GetBinContent(h_denominator_MC_pos->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));	
					float Multiplicative1 = 1 + h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]))) / h_denominator_MC_neg->GetBinContent(h_denominator_MC_neg->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));
					new_pt0 = (float) 1 / (Multiplicative0 / muon_pt->at(muon_hlt_match[0]) + Additive0);
					new_pt1 = (float) 1 / (Multiplicative1 / muon_pt->at(muon_hlt_match[1]) - Additive1);
// 					if(entry < 200){
// 						std::cout<<"MC: "<<std::endl;
// 						std::cout<<Additive0<<"\t"<<Multiplicative0<<"\t"<<new_pt0<<"\t"<<muon_pt->at(muon_hlt_match[0])<<std::endl;
// 						std::cout<<Additive1<<"\t"<<Multiplicative1<<"\t"<<new_pt1<<"\t"<<muon_pt->at(muon_hlt_match[1])<<std::endl;
// 						std::cout<<"+"<<h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])))<<"\t"<<h_denominator_MC_pos->GetBinContent(h_denominator_MC_pos->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])))<<std::endl;
// 						std::cout<<"-"<<h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])))<<"\t"<<h_denominator_MC_neg->GetBinContent(h_denominator_MC_neg->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])))<<std::endl;
// 					}
				}
				else{
					float Additive0 = h_Additive_DATA->GetBinContent(h_Additive_DATA->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));	
					float Additive1 = h_Additive_DATA->GetBinContent(h_Additive_DATA->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));	
// 					float Multiplicative0 = h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));
// 					float Multiplicative1 = h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));
					float Multiplicative0 = 1 + h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]))) / h_denominator_DATA_pos->GetBinContent(h_denominator_DATA_pos->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));	
					float Multiplicative1 = 1 + h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]))) / h_denominator_DATA_neg->GetBinContent(h_denominator_DATA_neg->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));
					new_pt0 = (float) 1 / (Multiplicative0 / muon_pt->at(muon_hlt_match[0]) + Additive0);
					new_pt1 = (float) 1 / (Multiplicative1 / muon_pt->at(muon_hlt_match[1]) - Additive1);
// 					if(entry < 200){
// 						std::cout<<"DATA: "<<std::endl;
// 						std::cout<<Additive0<<"\t"<<Multiplicative0<<"\t"<<new_pt0<<"\t"<<muon_pt->at(muon_hlt_match[0])<<std::endl;
// 						std::cout<<Additive1<<"\t"<<Multiplicative1<<"\t"<<new_pt1<<"\t"<<muon_pt->at(muon_hlt_match[1])<<std::endl;
// 						std::cout<<"+"<<h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])))<<"\t"<<h_denominator_DATA_pos->GetBinContent(h_denominator_DATA_pos->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])))<<std::endl;
// 						std::cout<<"-"<<h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])))<<"\t"<<h_denominator_DATA_neg->GetBinContent(h_denominator_DATA_neg->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])))<<std::endl;
// 					}






				}
			}
			else{
				if(s == 0){
					float Additive0 = h_Additive_MC->GetBinContent(h_Additive_MC->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));	
					float Additive1 = h_Additive_MC->GetBinContent(h_Additive_MC->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));	
// 					float Multiplicative0 = h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));
// 					float Multiplicative1 = h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));
					float Multiplicative0 = 1 + h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]))) / h_denominator_MC_neg->GetBinContent(h_denominator_MC_neg->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));	
					float Multiplicative1 = 1 + h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]))) / h_denominator_MC_pos->GetBinContent(h_denominator_MC_pos->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));
					new_pt0 = (float) 1 / (Multiplicative0 / muon_pt->at(muon_hlt_match[0]) - Additive0);
					new_pt1 = (float) 1 / (Multiplicative1 / muon_pt->at(muon_hlt_match[1]) + Additive1);
// 					if(entry < 200){
// 						std::cout<<"MC: "<<std::endl;
// 						std::cout<<Additive0<<"\t"<<Multiplicative0<<"\t"<<new_pt0<<"\t"<<muon_pt->at(muon_hlt_match[0])<<std::endl;
// 						std::cout<<Additive1<<"\t"<<Multiplicative1<<"\t"<<new_pt1<<"\t"<<muon_pt->at(muon_hlt_match[1])<<std::endl;
// 						std::cout<<"+"<<h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])))<<"\t"<<h_denominator_MC_pos->GetBinContent(h_denominator_MC_pos->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])))<<std::endl;
// 						std::cout<<"-"<<h_Multiplicative_MC->GetBinContent(h_Multiplicative_MC->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])))<<"\t"<<h_denominator_MC_neg->GetBinContent(h_denominator_MC_neg->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])))<<std::endl;
// 					}
				}
				else{
					float Additive0 = h_Additive_DATA->GetBinContent(h_Additive_DATA->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));	
					float Additive1 = h_Additive_DATA->GetBinContent(h_Additive_DATA->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));	
// 					float Multiplicative0 = h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));
// 					float Multiplicative1 = h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));
					float Multiplicative0 = 1 + h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]))) / h_denominator_DATA_neg->GetBinContent(h_denominator_DATA_neg->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])));	
					float Multiplicative1 = 1 + h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]))) / h_denominator_DATA_pos->GetBinContent(h_denominator_DATA_pos->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])));
					new_pt0 = (float) 1 / (Multiplicative0 / muon_pt->at(muon_hlt_match[0]) - Additive0);
					new_pt1 = (float) 1 / (Multiplicative1 / muon_pt->at(muon_hlt_match[1]) + Additive1);
// 					if(entry < 200){
// 						std::cout<<"DATA: "<<std::endl;
// 						std::cout<<Additive0<<"\t"<<Multiplicative0<<"\t"<<new_pt0<<"\t"<<muon_pt->at(muon_hlt_match[0])<<std::endl;
// 						std::cout<<Additive1<<"\t"<<Multiplicative1<<"\t"<<new_pt1<<"\t"<<muon_pt->at(muon_hlt_match[1])<<std::endl;
// 						std::cout<<"+"<<h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])))<<"\t"<<h_denominator_DATA_pos->GetBinContent(h_denominator_DATA_pos->FindBin(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0])))<<std::endl;
// 						std::cout<<"-"<<h_Multiplicative_DATA->GetBinContent(h_Multiplicative_DATA->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])))<<"\t"<<h_denominator_DATA_neg->GetBinContent(h_denominator_DATA_neg->FindBin(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1])))<<std::endl;
// 					}




				}
			}
// 			if(entry < 200){
// 				std::cout<<"Pt0 = "<<muon_pt->at(muon_hlt_match[0])<<"\t"<<new_pt0<<std::endl;
// 				std::cout<<"Pt1 = "<<muon_pt->at(muon_hlt_match[1])<<"\t"<<new_pt1<<std::endl;				
// 			}

// 			if(fabs(new_pt0 - muon_pt->at(muon_hlt_match[0])) / muon_pt->at(muon_hlt_match[0]) > 0.015) continue;
// 			if(fabs(new_pt1 - muon_pt->at(muon_hlt_match[1])) / muon_pt->at(muon_hlt_match[1]) > 0.015) continue;

			lep0.SetPtEtaPhiM(muon_pt->at(muon_hlt_match[0]), muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]), muon_mass->at(muon_hlt_match[0]));
			lep1.SetPtEtaPhiM(muon_pt->at(muon_hlt_match[1]), muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]), muon_mass->at(muon_hlt_match[1]));
			float unCorrectedMass = (lep0+lep1).M();
// 			Z = lep0+lep1;

	
			lep0.SetPtEtaPhiM(new_pt0, muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]), muon_mass->at(muon_hlt_match[0]));
			lep1.SetPtEtaPhiM(new_pt1, muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]), muon_mass->at(muon_hlt_match[1]));
			Z = lep0+lep1;
			
			if(Z.M() > 81 && Z.M()<101){	
// 			if(Z.M() > 60 && Z.M()<120){	

				if(s == 0){
					ZbosonPt = h_ZBoson_pt_correction->GetBinContent(h_ZBoson_pt_correction->FindBin(Z.Pt()));
					weight = lumi*crossSection/num_event;

					reweight = ZbosonPt * weight;
				}
				else{
					weight = 1;
					reweight = 1;
				}


				h_Zboson_mass[s]->Fill(Z.M(), reweight);
				h_Zboson_pT[s]->Fill(Z.Pt(), reweight);

				h_lepton_qOverPt[s]->Fill(muon_charge->at(muon_hlt_match[0]) / new_pt0, reweight);
				h_lepton_qOverPt[s]->Fill(muon_charge->at(muon_hlt_match[1]) / new_pt1, reweight);

				h_lepton_pt[s]->Fill(new_pt0, reweight);
				h_lepton_pt[s]->Fill(new_pt1, reweight);

				h_lepton_eta[s]->Fill(muon_eta->at(muon_hlt_match[0]), reweight);
				h_lepton_eta[s]->Fill(muon_eta->at(muon_hlt_match[1]), reweight);

				h_lepton_phi[s]->Fill(muon_phi->at(muon_hlt_match[0]), reweight);
				h_lepton_phi[s]->Fill(muon_phi->at(muon_hlt_match[1]), reweight);

				if(s == 0)
					h_MC_mass->Fill(unCorrectedMass, Z.M(), reweight);
				else
					h_DATA_mass->Fill(unCorrectedMass, Z.M(), reweight);


				if(muon_charge->at(muon_hlt_match[0]) > 0){
					h_lepton_qOverPt_scan_EtaPhi[s][0]->Fill(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]), (float)1 / new_pt0, reweight);
					h_lepton_qOverPt_scan_EtaPhi[s][1]->Fill(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]), (float)1 / new_pt1, reweight);

					h_Zboson_mass_vs_eta[s][0]->Fill(muon_eta->at(muon_hlt_match[0]), Z.M(), reweight);
					h_Zboson_mass_vs_eta[s][1]->Fill(muon_eta->at(muon_hlt_match[1]), Z.M(), reweight);
					h_Zboson_mass_vs_phi[s][0]->Fill(muon_phi->at(muon_hlt_match[0]), Z.M(), reweight);
					h_Zboson_mass_vs_phi[s][1]->Fill(muon_phi->at(muon_hlt_match[1]), Z.M(), reweight);
					if(s == 0){
						h_pTVariance_MC_pos->Fill((new_pt0 - muon_pt->at(muon_hlt_match[0]))/muon_pt->at(muon_hlt_match[0]), reweight);
						h_pTVariance_MC_neg->Fill((new_pt1 - muon_pt->at(muon_hlt_match[1]))/muon_pt->at(muon_hlt_match[1]), reweight);
					}
					if(s == 1){
						h_pTVariance_DATA_pos->Fill((new_pt0 - muon_pt->at(muon_hlt_match[0]))/muon_pt->at(muon_hlt_match[0]), reweight);
						h_pTVariance_DATA_neg->Fill((new_pt1 - muon_pt->at(muon_hlt_match[1]))/muon_pt->at(muon_hlt_match[1]), reweight);
					}

				}
				else{
					h_lepton_qOverPt_scan_EtaPhi[s][0]->Fill(muon_eta->at(muon_hlt_match[1]), muon_phi->at(muon_hlt_match[1]), (float)1 / new_pt1, reweight);
					h_lepton_qOverPt_scan_EtaPhi[s][1]->Fill(muon_eta->at(muon_hlt_match[0]), muon_phi->at(muon_hlt_match[0]), (float)1 / new_pt0, reweight);

					h_Zboson_mass_vs_eta[s][0]->Fill(muon_eta->at(muon_hlt_match[1]), Z.M(), reweight);
					h_Zboson_mass_vs_eta[s][1]->Fill(muon_eta->at(muon_hlt_match[0]), Z.M(), reweight);
					h_Zboson_mass_vs_phi[s][0]->Fill(muon_phi->at(muon_hlt_match[1]), Z.M(), reweight);
					h_Zboson_mass_vs_phi[s][1]->Fill(muon_phi->at(muon_hlt_match[0]), Z.M(), reweight);
					if(s == 0){
						h_pTVariance_MC_neg->Fill((new_pt0 - muon_pt->at(muon_hlt_match[0]))/muon_pt->at(muon_hlt_match[0]), reweight);
						h_pTVariance_MC_pos->Fill((new_pt1 - muon_pt->at(muon_hlt_match[1]))/muon_pt->at(muon_hlt_match[1]), reweight);
					}
					if(s == 1){
						h_pTVariance_DATA_neg->Fill((new_pt0 - muon_pt->at(muon_hlt_match[0]))/muon_pt->at(muon_hlt_match[0]), reweight);
						h_pTVariance_DATA_pos->Fill((new_pt1 - muon_pt->at(muon_hlt_match[1]))/muon_pt->at(muon_hlt_match[1]), reweight);
					}
				}	



				if(s == 0){
					if(muon_GEN_match.size() < 2) continue;
					if(GEN_id->at(muon_GEN_match[0]) + GEN_id->at(muon_GEN_match[1]) != 0) continue;

					float GEN_pt0;
					rand.SetSeed(abs(static_cast<int>(sin(muon_phi->at(muon_hlt_match[0]))*100000)));
					GEN_pt0 = GEN_pt->at(muon_GEN_match[0]) * (1 + rand.Gaus(0, h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[0])))));
// 					GEN_pt0 = rand.Gaus(GEN_pt->at(muon_GEN_match[0]), GEN_pt->at(muon_GEN_match[0]) * h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[0]))));
// 					std::cout<<GEN_pt->at(muon_GEN_match[0])<<"\t"<<h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[0]))<<"\t"<<h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[0])))<<"\t"<<GEN_pt0<<std::endl;
					
					float GEN_pt1;
					rand.SetSeed(abs(static_cast<int>(sin(muon_phi->at(muon_hlt_match[1]))*100000)));
					GEN_pt1 = GEN_pt->at(muon_GEN_match[1]) * (1 + rand.Gaus(0, h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[1])))));
// 					GEN_pt1 = rand.Gaus(GEN_pt->at(muon_GEN_match[1]), GEN_pt->at(muon_GEN_match[1]) * h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[1]))));
// 					std::cout<<GEN_pt->at(muon_GEN_match[1])<<"\t"<<h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[1]))<<"\t"<<h_resGEN->GetBinContent(h_resGEN->FindBin(GEN_pt->at(muon_GEN_match[1])))<<"\t"<<GEN_pt1<<std::endl;
			

					TLorentzVector lep0, lep1;
					lep0.SetPtEtaPhiM(GEN_pt0, GEN_eta->at(muon_GEN_match[0]), GEN_phi->at(muon_GEN_match[0]), GEN_mass->at(muon_GEN_match[0]));
					lep1.SetPtEtaPhiM(GEN_pt1, GEN_eta->at(muon_GEN_match[1]), GEN_phi->at(muon_GEN_match[1]), GEN_mass->at(muon_GEN_match[1]));
					TLorentzVector Z = lep0+lep1;

					h_Zboson_mass[2]->Fill(Z.M(), reweight);
					h_Zboson_pT[2]->Fill(Z.Pt(), reweight);

					h_lepton_qOverPt[2]->Fill(GEN_id->at(muon_GEN_match[0]) / ((float)13 * GEN_pt0), reweight);
					h_lepton_qOverPt[2]->Fill(GEN_id->at(muon_GEN_match[1]) / ((float)13 * GEN_pt1), reweight);

					h_lepton_pt[2]->Fill(GEN_pt0, reweight);
					h_lepton_pt[2]->Fill(GEN_pt1, reweight);

					h_lepton_eta[2]->Fill(GEN_eta->at(muon_GEN_match[0]), reweight);
					h_lepton_eta[2]->Fill(GEN_eta->at(muon_GEN_match[1]), reweight);

					h_lepton_phi[2]->Fill(GEN_phi->at(muon_GEN_match[0]), reweight);
					h_lepton_phi[2]->Fill(GEN_phi->at(muon_GEN_match[1]), reweight);

					if(GEN_id->at(muon_GEN_match[0]) < 0){
						h_lepton_qOverPt_scan_EtaPhi[2][0]->Fill(GEN_eta->at(muon_GEN_match[0]), GEN_phi->at(muon_GEN_match[0]), (float)1 / ( GEN_pt0), reweight);
						h_lepton_qOverPt_scan_EtaPhi[2][1]->Fill(GEN_eta->at(muon_GEN_match[1]), GEN_phi->at(muon_GEN_match[1]), (float)1 / ( GEN_pt1), reweight);

						h_Zboson_mass_vs_eta[2][0]->Fill(muon_eta->at(muon_hlt_match[0]), Z.M(), reweight);
						h_Zboson_mass_vs_eta[2][1]->Fill(muon_eta->at(muon_hlt_match[1]), Z.M(), reweight);
						h_Zboson_mass_vs_phi[2][0]->Fill(muon_phi->at(muon_hlt_match[0]), Z.M(), reweight);
						h_Zboson_mass_vs_phi[2][1]->Fill(muon_phi->at(muon_hlt_match[1]), Z.M(), reweight);
					}
					else{
						h_lepton_qOverPt_scan_EtaPhi[2][0]->Fill(GEN_eta->at(muon_GEN_match[1]), GEN_phi->at(muon_GEN_match[1]), (float)1 / ( GEN_pt1), reweight);
						h_lepton_qOverPt_scan_EtaPhi[2][1]->Fill(GEN_eta->at(muon_GEN_match[0]), GEN_phi->at(muon_GEN_match[0]), (float)1 / ( GEN_pt0), reweight);

						h_Zboson_mass_vs_eta[2][0]->Fill(muon_eta->at(muon_hlt_match[1]), Z.M(), reweight);
						h_Zboson_mass_vs_eta[2][1]->Fill(muon_eta->at(muon_hlt_match[0]), Z.M(), reweight);
						h_Zboson_mass_vs_phi[2][0]->Fill(muon_phi->at(muon_hlt_match[1]), Z.M(), reweight);
						h_Zboson_mass_vs_phi[2][1]->Fill(muon_phi->at(muon_hlt_match[0]), Z.M(), reweight);
					}		
				}
		
			}


		}	
	}
		
	
	
// 	DrawResolution(h_Zboson_mass[0], h_Zboson_mass[1], h_Zboson_mass[2], samples, "Zboson_mass", save + "[", "mass_{2l} (GeV)");
	DrawResolution(h_Zboson_mass[0], h_Zboson_mass[1], h_Zboson_mass[2], samples, "Zboson_mass", save, "mass_{2l} (GeV)");
	DrawResolution(h_Zboson_pT[0], h_Zboson_pT[1], h_Zboson_pT[2], samples, "Zboson_pT", save, "pT_{2l} (GeV)");
	DrawResolution(h_lepton_pt[0], h_lepton_pt[1], h_lepton_pt[2], samples, "h_lepton_pt", save, "p_{T} (GeV)");
	DrawResolution(h_lepton_eta[0], h_lepton_eta[1], h_lepton_eta[2], samples, "h_lepton_eta", save, "#eta");
	DrawResolution(h_lepton_phi[0], h_lepton_phi[1], h_lepton_phi[2], samples, "h_lepton_phi", save, "#phi");
	DrawResolution(h_lepton_qOverPt[0], h_lepton_qOverPt[1], h_lepton_qOverPt[2], samples, "h_lepton_qOverPt", save, "q/p_{T} (GeV)");
// 	DrawResolution(h_lepton_qOverPt[0], h_lepton_qOverPt[1], h_lepton_qOverPt[2], samples, "h_lepton_qOverPt", save + "]", "q/p_{T} (GeV)");

	TCanvas* c_MC_mass = new TCanvas("c_MC_mass", "c_MC_mass", 600, 600);
	h_MC_mass->Draw("COLZ");
	c_MC_mass->Print(save);
	TCanvas* c_DATA_mass = new TCanvas("c_DATA_mass", "c_DATA_mass", 600, 600);
	h_DATA_mass->Draw("COLZ");
	c_DATA_mass->Print(save);
// 	c_DATA_mass->Print(save + "]");


	Mass_vs_phi_pos = new TMultiGraph();
	Mass_vs_phi_neg = new TMultiGraph();
	Mass_vs_eta_pos = new TMultiGraph();
	Mass_vs_eta_neg = new TMultiGraph();

	for(int s = 0; s < 3; s++){
		std::vector<Double_t> mean_pos;
		std::vector<Double_t> mean_neg;
		std::vector<Double_t> eta_bins_centered;
		std::vector<Double_t> phi_bins_centered;
		for(int xxx = 1; xxx < h_Zboson_mass_vs_eta[s][0]->GetNbinsX()+1; xxx++){
			TString name = Form("Mass_eta_%.1f_%.1f_pos", eta_bins.at(xxx-1), eta_bins.at(xxx));
// 			std::cout<<"ETA_name = "<<name<<std::endl;

			TH1F* h_tmp_MC	= (TH1F*)h_Zboson_mass_vs_eta[s][0]->ProjectionY(name, xxx, xxx);		
    	    RooRealVar massZ_P("massZ","massZ",81,101);
			RooDataHist histo("mass Z", "mass Z", massZ_P, h_tmp_MC);
    
    	    RooRealVar PDGmZ("PDGmZ","PDGmZ", 91.19, 86, 96);
	        RooRealVar PDGwZ("PDGwZ","PDGwZ", 2.5, 1, 4);
	        RooBreitWigner PDGBW("PDGBW","PDGBW",massZ_P,PDGmZ,PDGwZ);
	        RooBreitWigner PDGBW_2("PDGBW_2","PDGBW_2",massZ_P,PDGmZ,PDGwZ);
        
			RooRealVar CB_mean("CB_mean", "CB_mean", 0, -5, 5);//125, 120, 130);
			RooRealVar Sigma("Sigma", "Sigma", 1, 0.1, 20);//1, 0, 30);//sigma[decay]);
			RooRealVar AlphaL("AlphaL", "AlphaL", 1, 0.1, 30);//alphaL[decay]);
			RooRealVar ExpL("ExpL", "ExpL", 1, 0.1, 30);//expL[decay]);
			RooRealVar AlphaR("AlphaR", "AlphaR", 1, 0.1, 30);//alphaR[decay]);
			RooRealVar ExpR("ExpR", "ExpR", 1, 0.1, 50);//expR[decay]);
			RooCrystalBall DSCB("DSCB", "DSCB", massZ_P, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    
	        RooFFTConvPdf model("CW","CW",massZ_P,PDGBW,DSCB);
	        RooFFTConvPdf model_2("CW_2","CW_2",massZ_P,PDGBW,DSCB);
	        
			if(s != 20){
				model.fitTo(histo, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			else{
				PDGBW.fitTo(histo, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			mean_pos.push_back(PDGmZ.getVal());
			eta_bins_centered.push_back(eta_bins.at(xxx-1) + fabs((eta_bins.at(1) - eta_bins.at(0))/2));


			TCanvas *c_MC = new TCanvas(name, name, 900, 700);
			c_MC->SetFrameFillColor(0);
			RooPlot* xframe = massZ_P.frame("massZ");
			xframe = massZ_P.frame(Title(name));
			histo.plotOn(xframe);
			if(s != 20){
				model.plotOn(xframe,RooFit::LineColor(kBlue));
				model.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));
			}
			else{
				PDGBW.plotOn(xframe,RooFit::LineColor(kRed));
				PDGBW.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));				
			}
			xframe->Draw();
// 			c_MC->Print(save);
			
			
			name = Form("Mass_eta_%.1f_%.1f_neg", eta_bins.at(xxx-1), eta_bins.at(xxx));
			TH1F* h_tmp_MC_2	= (TH1F*)h_Zboson_mass_vs_eta[s][1]->ProjectionY(name, xxx, xxx);
			RooDataHist histo_neg("mass Z", "mass Z", massZ_P, h_tmp_MC_2);
			if(s != 20){
				model_2.fitTo(histo_neg, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			else{
				PDGBW_2.fitTo(histo_neg, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			mean_neg.push_back(PDGmZ.getVal());
			
			TCanvas* c_MC_2 = new TCanvas(name, name, 900, 700);
			c_MC_2->SetFrameFillColor(0);
			xframe = massZ_P.frame("massZ");
			xframe = massZ_P.frame(Title(name));
			histo_neg.plotOn(xframe);
			if(s != 20){
				model_2.plotOn(xframe,RooFit::LineColor(kBlue));
				model_2.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));
			}
			else{
				PDGBW_2.plotOn(xframe,RooFit::LineColor(kRed));
				PDGBW_2.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));				
			}
			xframe->Draw();
// 			c_MC_2->Print(save);		
		}
		
		TGraph* gr_new = new TGraph(eta_bins.size()-1, &eta_bins_centered[0], &mean_pos[0]);
		if(s == 1)
			gr_new->SetMarkerColor(kRed);
		if(s == 2)
			gr_new->SetMarkerColor(kGreen);
		Mass_vs_eta_pos->Add(gr_new);
		
		gr_new = new TGraph(eta_bins.size()-1, &eta_bins_centered[0], &mean_neg[0]);
		if(s == 1)
			gr_new->SetMarkerColor(kRed);
		if(s == 2)
			gr_new->SetMarkerColor(kGreen);
		Mass_vs_eta_neg->Add(gr_new);
		
		
		mean_pos.clear();
		mean_neg.clear();	
		
		for(int xxx = 1; xxx < h_Zboson_mass_vs_phi[s][0]->GetNbinsX()+1; xxx++){
			TString name = Form("Mass_phi_%.1f_%.1f_pos", phi_bins.at(xxx-1), phi_bins.at(xxx));
// 			std::cout<<"PHI_name = "<<name<<std::endl;
			TH1F* h_tmp_MC	= (TH1F*)h_Zboson_mass_vs_phi[s][0]->ProjectionY(name, xxx, xxx);		
    	    RooRealVar massZ_P("massZ","massZ",81,101);
			RooDataHist histo("mass Z", "mass Z", massZ_P, h_tmp_MC);
    
    	    RooRealVar PDGmZ("PDGmZ","PDGmZ",91.19, 86, 96);
	        RooRealVar PDGwZ("PDGwZ","PDGwZ", 2.5, 1, 4);
	        RooBreitWigner PDGBW("PDGBW","PDGBW",massZ_P,PDGmZ,PDGwZ);
 	        RooBreitWigner PDGBW_2("PDGBW_2","PDGBW_2",massZ_P,PDGmZ,PDGwZ);
       
			RooRealVar CB_mean("CB_mean", "CB_mean", 0, -5, 5);//125, 120, 130);
			RooRealVar Sigma("Sigma", "Sigma", 1, 0.1, 20);//1, 0, 30);//sigma[decay]);
			RooRealVar AlphaL("AlphaL", "AlphaL", 1, 0.1, 30);//alphaL[decay]);
			RooRealVar ExpL("ExpL", "ExpL", 1, 0.1, 30);//expL[decay]);
			RooRealVar AlphaR("AlphaR", "AlphaR", 1, 0.1, 30);//alphaR[decay]);
			RooRealVar ExpR("ExpR", "ExpR", 1, 0.1, 50);//expR[decay]);
			RooCrystalBall DSCB("DSCB", "DSCB", massZ_P, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    
	        RooFFTConvPdf model("CW","CW",massZ_P,PDGBW,DSCB);
	        RooFFTConvPdf model_2("CW_2","CW_2",massZ_P,PDGBW,DSCB);
	        
			if(s != 20){
				model.fitTo(histo, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			else{
				PDGBW.fitTo(histo, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			mean_pos.push_back(PDGmZ.getVal());
			phi_bins_centered.push_back(phi_bins.at(xxx-1) + fabs((phi_bins.at(1) - phi_bins.at(0))/2));

			TCanvas *c_MC = new TCanvas(name, name, 900, 700);
			c_MC->SetFrameFillColor(0);
			RooPlot* xframe = massZ_P.frame("massZ");
			xframe = massZ_P.frame(Title(name));
			histo.plotOn(xframe);
			if(s != 20){
				model.plotOn(xframe,RooFit::LineColor(kBlue));
				model.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));
			}
			else{
				PDGBW.plotOn(xframe,RooFit::LineColor(kRed));
				PDGBW.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));				
			}
			xframe->Draw();
// 			c_MC->Print(save);
					
			name = Form("Mass_phi_%.1f_%.1f_neg", phi_bins.at(xxx-1), phi_bins.at(xxx));
			TH1F* h_tmp_MC_2	= (TH1F*)h_Zboson_mass_vs_phi[s][1]->ProjectionY(name, xxx, xxx);
			RooDataHist histo_neg("mass Z", "mass Z", massZ_P, h_tmp_MC_2);
			if(s != 20){
				model_2.fitTo(histo_neg, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			else{
				PDGBW_2.fitTo(histo_neg, Range(86, 96), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
			}
			mean_neg.push_back(PDGmZ.getVal());

			c_MC = new TCanvas(name, name, 900, 700);
			c_MC->SetFrameFillColor(0);
			xframe = massZ_P.frame("massZ");
			xframe = massZ_P.frame(Title(name));
			histo_neg.plotOn(xframe);
			if(s != 20){
				model_2.plotOn(xframe,RooFit::LineColor(kBlue));
				model_2.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));
			}
			else{
				PDGBW_2.plotOn(xframe,RooFit::LineColor(kRed));
				PDGBW_2.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));				
			}
			xframe->Draw();
// 			c_MC->Print(save);		
		}
		
		gr_new = new TGraph(phi_bins.size()-1, &phi_bins_centered[0], &mean_pos[0]);
		if(s == 1)
			gr_new->SetMarkerColor(kRed);
		if(s == 2)
			gr_new->SetMarkerColor(kGreen);
		Mass_vs_phi_pos->Add(gr_new);
		
		gr_new = new TGraph(phi_bins.size()-1, &phi_bins_centered[0], &mean_neg[0]);
		if(s == 1)
			gr_new->SetMarkerColor(kRed);
		if(s == 2)
			gr_new->SetMarkerColor(kGreen);
		Mass_vs_phi_neg->Add(gr_new);
			
	}

	for(int xxx = 1; xxx < h_Zboson_mass_vs_eta[0][0]->GetNbinsX()+1; xxx++){
		TString name = Form("Mass_eta_%.1f_%.1f_pos_MC", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_MC	= (TH1F*)h_Zboson_mass_vs_eta[0][0]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_eta_%.1f_%.1f_pos_DATA", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_DATA	= (TH1F*)h_Zboson_mass_vs_eta[1][0]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_eta_%.1f_%.1f_pos", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_GEN	= (TH1F*)h_Zboson_mass_vs_eta[2][0]->ProjectionY(name, xxx, xxx);		

		DrawResolution(h_tmp_MC, h_tmp_DATA, h_tmp_GEN, samples, name, save, "mass_{2l} (GeV)");

		name = Form("Mass_eta_%.1f_%.1f_neg_MC", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_MC_2	= (TH1F*)h_Zboson_mass_vs_eta[0][1]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_eta_%.1f_%.1f_neg_DATA", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_DATA_2	= (TH1F*)h_Zboson_mass_vs_eta[1][1]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_eta_%.1f_%.1f_neg", eta_bins.at(xxx-1), eta_bins.at(xxx));
		TH1F* h_tmp_GEN_2	= (TH1F*)h_Zboson_mass_vs_eta[2][1]->ProjectionY(name, xxx, xxx);	
		name = Form("Mass_eta_%.1f_%.1f_neg", eta_bins.at(xxx-1), eta_bins.at(xxx));

		DrawResolution(h_tmp_MC_2, h_tmp_DATA_2, h_tmp_GEN_2, samples, name, save, "mass_{2l} (GeV)");
	}					
	for(int xxx = 1; xxx < h_Zboson_mass_vs_phi[0][0]->GetNbinsX()+1; xxx++){
		TString name = Form("Mass_phi_%.1f_%.1f_pos_MC", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_MC	= (TH1F*)h_Zboson_mass_vs_phi[0][0]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_phi_%.1f_%.1f_pos_DATA", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_DATA	= (TH1F*)h_Zboson_mass_vs_phi[1][0]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_phi_%.1f_%.1f_pos", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_GEN	= (TH1F*)h_Zboson_mass_vs_phi[2][0]->ProjectionY(name, xxx, xxx);		

		DrawResolution(h_tmp_MC, h_tmp_DATA, h_tmp_GEN, samples, name, save, "mass_{2l} (GeV)");

		name = Form("Mass_phi_%.1f_%.1f_neg_MC", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_MC_2	= (TH1F*)h_Zboson_mass_vs_phi[0][1]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_phi_%.1f_%.1f_neg_DATA", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_DATA_2	= (TH1F*)h_Zboson_mass_vs_phi[1][1]->ProjectionY(name, xxx, xxx);		
		name = Form("Mass_phi_%.1f_%.1f_neg", phi_bins.at(xxx-1), phi_bins.at(xxx));
		TH1F* h_tmp_GEN_2	= (TH1F*)h_Zboson_mass_vs_phi[2][1]->ProjectionY(name, xxx, xxx);	
		name = Form("Mass_phi_%.1f_%.1f_neg", phi_bins.at(xxx-1), phi_bins.at(xxx));

		DrawResolution(h_tmp_MC_2, h_tmp_DATA_2, h_tmp_GEN_2, samples, name, save, "mass_{2l} (GeV)");
		
	}	

	c_Mass_vs_eta_pos = new TCanvas("c_Mass_vs_eta_pos", "c_Mass_vs_eta_pos", 600, 600);
	Mass_vs_eta_pos->SetTitle("c_Mass_vs_eta_pos");
	Mass_vs_eta_pos->Draw("AP*");
	c_Mass_vs_eta_pos->Print(save);
	c_Mass_vs_eta_neg = new TCanvas("c_Mass_vs_eta_neg", "c_Mass_vs_eta_neg", 600, 600);
	Mass_vs_eta_neg->SetTitle("c_Mass_vs_eta_neg");
	Mass_vs_eta_neg->Draw("AP*");
	c_Mass_vs_eta_neg->Print(save);
	c_Mass_vs_phi_pos = new TCanvas("c_Mass_vs_phi_pos", "c_Mass_vs_phi_pos", 600, 600);
	Mass_vs_phi_pos->SetTitle("c_Mass_vs_phi_pos");
	Mass_vs_phi_pos->Draw("AP*");
	c_Mass_vs_phi_pos->Print(save);
	c_Mass_vs_phi_neg = new TCanvas("c_Mass_vs_phi_neg", "c_Mass_vs_phi_neg", 600, 600);
	Mass_vs_phi_neg->SetTitle("c_Mass_vs_phi_neg");
	Mass_vs_phi_neg->Draw("AP*");
	c_Mass_vs_phi_neg->Print(save);
// 	c_Mass_vs_phi_neg->Print(save + "]");




	h_Additive_MC->Reset();
	h_Multiplicative_MC->Reset();
	h_Additive_DATA->Reset();
	h_Multiplicative_DATA->Reset();

	for(int xxx = 1; xxx < h_lepton_qOverPt_scan_EtaPhi[0][0]->GetNbinsX()+1; xxx++){
		for(int yyy = 1; yyy < h_lepton_qOverPt_scan_EtaPhi[0][0]->GetNbinsY()+1; yyy++){
			TString name = Form("1overPt_eta_%.1f_%.1f_phi_%.1f_%.1f", eta_bins.at(xxx-1), eta_bins.at(xxx), phi_bins.at(yyy-1), phi_bins.at(yyy));
			std::cout<<"========"<<std::endl;
			std::cout<<name<<std::endl;
			TH1F* h_tmp_MC	= (TH1F*)h_lepton_qOverPt_scan_EtaPhi[0][0]->ProjectionZ(name + "_MC+", xxx, xxx, yyy, yyy);
			TH1F* h_MC = (TH1F*)h_tmp_MC->Clone();
			TH1F* h_tmp_DATA	= (TH1F*)h_lepton_qOverPt_scan_EtaPhi[1][0]->ProjectionZ(name + "_DATA+", xxx, xxx, yyy, yyy);
			TH1F* h_DATA = (TH1F*)h_tmp_DATA->Clone();
			TH1F* h_tmp_GEN	= (TH1F*)h_lepton_qOverPt_scan_EtaPhi[2][0]->ProjectionZ(name + "_GEN+", xxx, xxx, yyy, yyy);
			TH1F* h_GEN = (TH1F*)h_tmp_GEN->Clone();
			h_denominator_MC_pos->SetBinContent(xxx, yyy, h_MC->GetMean());
			h_denominator_DATA_pos->SetBinContent(xxx, yyy, h_DATA->GetMean());

			std::cout<<"After\nDenominator+ ="<<h_MC->GetMean()<<"\t"<<h_DATA->GetMean()<<std::endl;

// 			DrawResolution(h_MC, h_DATA, h_GEN, samples, name + " pos", save, "1 / p_{T} (GeV)");

			Corr_MC_pos = h_GEN->GetMean() - h_MC->GetMean();
			Corr_DATA_pos = h_GEN->GetMean() - h_DATA->GetMean();
			std::cout<<"Mean+  = "<<h_GEN->GetMean()<<"\t"<<h_MC->GetMean()<<"\t"<<h_DATA->GetMean()<<std::endl;

			h_tmp_MC->Reset();
			h_tmp_DATA->Reset();
			h_tmp_GEN->Reset();
			h_MC->Reset();
			h_DATA->Reset();
			h_GEN->Reset();

			h_tmp_MC	= (TH1F*)h_lepton_qOverPt_scan_EtaPhi[0][1]->ProjectionZ(name + "_MC-", xxx, xxx, yyy, yyy);
			h_MC = (TH1F*)h_tmp_MC->Clone();
			h_tmp_DATA	= (TH1F*)h_lepton_qOverPt_scan_EtaPhi[1][1]->ProjectionZ(name + "_DATA-", xxx, xxx, yyy, yyy);
			h_DATA = (TH1F*)h_tmp_DATA->Clone();
			h_tmp_GEN	= (TH1F*)h_lepton_qOverPt_scan_EtaPhi[2][1]->ProjectionZ(name + "_GEN-", xxx, xxx, yyy, yyy);
			h_GEN = (TH1F*)h_tmp_GEN->Clone();
			h_denominator_MC_neg->SetBinContent(xxx, yyy, h_MC->GetMean());
			h_denominator_DATA_neg->SetBinContent(xxx, yyy, h_DATA->GetMean());
			std::cout<<"Denominator- ="<<h_MC->GetMean()<<"\t"<<h_DATA->GetMean()<<std::endl;
	
// 			DrawResolution(h_MC, h_DATA, h_GEN, samples, name + " neg", save, "1 / p_{T} (GeV)");
		
			Corr_MC_neg = h_GEN->GetMean() - h_MC->GetMean();
			Corr_DATA_neg = h_GEN->GetMean() - h_DATA->GetMean();
			std::cout<<"Mean-  = "<<h_GEN->GetMean()<<"\t"<<h_MC->GetMean()<<"\t"<<h_DATA->GetMean()<<std::endl;
// 			
			std::cout<<"C_MC = "<<Corr_MC_pos<<"\t"<<Corr_MC_neg<<std::endl;
			std::cout<<"C_DATA = "<<Corr_DATA_pos<<"\t"<<Corr_DATA_neg<<std::endl;

			D_m_MC = (Corr_MC_pos + Corr_MC_neg) / (float)2;
			D_a_MC = (Corr_MC_pos - Corr_MC_neg) / (float)2;
			D_m_DATA = (Corr_DATA_pos + Corr_DATA_neg) / (float)2;
			D_a_DATA = (Corr_DATA_pos - Corr_DATA_neg) / (float)2;
			std::cout<<"D_m_MC = "<<D_m_MC<<std::endl;
			std::cout<<"D_m_DATA = "<<D_m_DATA<<std::endl;
			std::cout<<"D_a_MC = "<<D_a_MC<<std::endl;
			std::cout<<"D_a_DATA = "<<D_a_DATA<<std::endl;

			float M = D_m_MC;
			float A = D_a_MC;
			h_Additive_MC->SetBinContent(xxx, yyy, A);
			h_Multiplicative_MC->SetBinContent(xxx, yyy, M);
			std::cout<<"MC: M = "<<M<<"\t A = "<<A<<std::endl;
// 			std::cout<<"Controllo MC = "<<A<<"\t"<<D_a_MC<<"\t"<<D_m_MC<<"\t"<<numerator_MC<<"\t"<<denominator_MC<<std::endl;
// 			std::cout<<"Controllo MC = "<<A<<"\t"<<M<<"\t"<<D_a_MC<<"\t"<<D_m_MC<<"\t"<<denominator_MC<<std::endl;

			M = D_m_DATA;
			A = D_a_DATA;
			h_Additive_DATA->SetBinContent(xxx, yyy, A);
			h_Multiplicative_DATA->SetBinContent(xxx, yyy, M);
			std::cout<<"DATA: M = "<<M<<"\t A = "<<A<<std::endl;
// 			std::cout<<"Controllo DATA = "<<A<<"\t"<<D_a_DATA<<"\t"<<D_m_DATA<<"\t"<<numerator_DATA<<"\t"<<denominator_DATA<<std::endl;
// 			std::cout<<"Controllo DATA = "<<A<<"\t"<<M<<"\t"<<D_a_DATA<<"\t"<<D_m_DATA<<"\t"<<denominator_DATA<<std::endl;
		}
	}			
// 	return;		

	c_Additive_MC = new TCanvas("c_Additive_MC", "c_Additive_MC", 600, 600);
	c_Additive_MC->SetRightMargin(0.2);
	h_Additive_MC->Draw("COLZ");
	c_Additive_MC->Print(save);
	c_Additive_DATA = new TCanvas("c_Additive_DATA", "c_Additive_DATA", 600, 600);
	c_Additive_DATA->SetRightMargin(0.2);
	h_Additive_DATA->Draw("COLZ");
	c_Additive_DATA->Print(save);
	c_Multiplicative_MC = new TCanvas("c_Multiplicative_MC", "c_Multiplicative_MC", 600, 600);
	c_Multiplicative_MC->SetRightMargin(0.2);
	h_Multiplicative_MC->Draw("COLZ");
	c_Multiplicative_MC->Print(save);
	c_Multiplicative_DATA = new TCanvas("c_Multiplicative_DATA", "c_Multiplicative_DATA", 600, 600);
	c_Multiplicative_DATA->SetRightMargin(0.2);
	h_Multiplicative_DATA->Draw("COLZ");
	c_Multiplicative_DATA->Print(save);
// 	c_Multiplicative_DATA->Print(save + "]");


	TCanvas* c_pTVariance_MC_pos = new TCanvas("c_pTVariance_MC_pos", "c_pTVariance_MC_pos", 600, 600);
	h_pTVariance_MC_pos->SetStats(1111);
	h_pTVariance_MC_pos->Draw();
	c_pTVariance_MC_pos->Print(save);
	TCanvas* c_pTVariance_MC_neg = new TCanvas("c_pTVariance_MC_neg", "c_pTVariance_MC_neg", 600, 600);
	h_pTVariance_MC_neg->SetStats(1111);
	h_pTVariance_MC_neg->Draw();
	c_pTVariance_MC_neg->Print(save);
	TCanvas* c_pTVariance_DATA_pos = new TCanvas("c_pTVariance_DATA_pos", "c_pTVariance_DATA_pos", 600, 600);
	h_pTVariance_DATA_pos->SetStats(1111);
	h_pTVariance_DATA_pos->Draw();
	c_pTVariance_DATA_pos->Print(save);
	TCanvas* c_pTVariance_DATA_neg = new TCanvas("c_pTVariance_DATA_neg", "c_pTVariance_DATA_neg", 600, 600);
	h_pTVariance_DATA_neg->SetStats(1111);
	h_pTVariance_DATA_neg->Draw();
	c_pTVariance_DATA_neg->Print(save);
	c_pTVariance_DATA_neg->Print(save + "]");
	return;


}

*/

float dataMCErr(float pt, float eta, TH2F* hMuScaleFacUnc)
{
    return hMuScaleFacUnc->GetBinContent(hMuScaleFacUnc->FindBin(eta,pt)); 
}
























