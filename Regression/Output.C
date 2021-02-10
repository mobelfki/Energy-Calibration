	#define Output_cxx
#include "Output.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
int _npar;
double FitFunction(double *x,double *par) 
{
    
      double fitval = par[0];
      for( int i = 1; i<=_npar; i++)
      {		
	fitval += par[i]*TMath::Power(TMath::Log10(x[0]),i);
      }
      return fitval;
}
void Output::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Output.C
//      root> Output t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
	std::vector<double> Rbins;
	double Rmax = 2.5;
	double Rmin = 0.0;
	int RN = 250;
	double dR=(Rmax-Rmin)/RN;
  	for (int i=0;i<=RN;++i)
	{ 
		Rbins.push_back(Rmin+i*dR);
	}
	double E [] = {20000,25000,30000,35000,40000,45000,50000,55000,60000,70000,80000,90000,100000,120000,140000,160000,180000,200000,240000,280000,320000,360000,400000,450000,500000};
	std::vector<double> E2bins;
	double Emin = 3.;
	double Emax = 7.;
	int NE2 = 400;
	double dE = (Emax-Emin)/NE2;
	for (int i = 0; i<=NE2; ++i)
	{
		E2bins.push_back(Emin+i*dE);
	} 
	float jetbinrange [] = {3.,3.05,3.1,3.15,3.2,3.25,3.3,3.35,3.4,3.45,3.5,3.55,3.6,3.65,3.7,3.75,3.8,3.85,3.9,3.95,4.,4.05,4.1,4.15,4.2,4.25,4.3,4.35,4.4,4.45,4.5,4.55,4.6,4.65,4.7,4.75,4.8,4.85,4.9,4.95,5.,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.,6.2,6.4,6.8};

	float jetbinrange2 [] = {3.25,3.3,3.35,3.4,3.45,3.5,3.55,3.6,3.65,3.7,3.75,3.8,3.85,3.9,3.95,4.,4.05,4.1,4.15,4.2,4.25,4.3,4.35,4.4,4.45,4.5,4.55,4.6,4.65,4.7,4.75,4.8,4.85,4.9,4.95,5.,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.};
	

	double Eta [] = {-2.5,-2.4,-2.3,-2.2,-2.1,-2.,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5};	
	int etabin = 50;

		//std::cout<<jetbinrange.size()<<std::endl;

		int nbins = 54;

		TH1F *hitsreco           = new TH1F("Correction_Inc_mean","",nbins-1,jetbinrange);
		TH1F *hitsreco_semi      = new TH1F("Correction_Semi_mean","",nbins-1,jetbinrange);
		TH1F *hitsreco_had       = new TH1F("Correction_Had_mean","",nbins-1,jetbinrange);

		TH1F *hitsreco_median           = new TH1F("Correction_Inc_median","",nbins-1,jetbinrange);
		TH1F *hitsreco_semi_median      = new TH1F("Correction_Semi_median","",nbins-1,jetbinrange);
		TH1F *hitsreco_had_median      = new TH1F("Correction_Had_median","",nbins-1,jetbinrange);



		TH1F* ptreco_histos_lowEta[nbins];
		TH1F* ptreco_histos_had_lowEta[nbins];
		TH1F* ptreco_histos_semi_lowEta[nbins];

		TProfile *Reco_vs_Truth_Inc     = new TProfile("Reco_vs_Truth_Inc","",E2bins.size()-1,&E2bins[0]);
		TProfile *Reco_vs_Truth_Semi  = new TProfile("Reco_vs_Truth_Semi","",E2bins.size()-1,&E2bins[0]);
		TProfile *Reco_vs_Truth_Had   = new TProfile("Reco_vs_Truth_Had","",E2bins.size()-1,&E2bins[0]);

		for (int n = 0; n < nbins; n++) {

			ptreco_histos_semi_lowEta[n]  = new TH1F(Form("hist_Semi_LogPt_bin_%2.2f_%2.2f",jetbinrange[n],jetbinrange[n+1]), Form("hist_Semi_LogPt_bin_%2.2f_%2.2f ; 1/RecoPt/TruthPt ; AU",jetbinrange[n],jetbinrange[n+1]),  40,0.,2.);
			ptreco_histos_had_lowEta[n]    = new TH1F(Form("hist_Had_LogPt_bin_%2.2f_%2.2f",jetbinrange[n],jetbinrange[n+1]),   Form("hist_Had_LogPt_bin_%2.2f_%2.2f ; 1/RecoPt/TruthPt ; AU",jetbinrange[n],jetbinrange[n+1]),    40,0.,2.);
			ptreco_histos_lowEta[n]            = new TH1F(Form("hist_LogInc_Pt_bin_%2.2f_%2.2f",jetbinrange[n],jetbinrange[n+1]),     Form("hist_Inc_LogPt_bin_%2.2f_%2.2f ; 1/RecoPt/TruthPt ; AU",jetbinrange[n],jetbinrange[n+1]),     40,0.,2.);

		}


	int NE = 25;
	std::vector<double> Ebins;
	for(int i = 0; i < NE; i++)
	{
		Ebins.push_back(E[i]);
	}
	for(int i = 0; i < NEta; i++)
	{
		R_vs_Truth[i] = new TH2F(Form("R_vs_Truth_Eta_%2.1f_%2.1f",Etabins[i],Etabins[i+1]),Form("R_vs_Truth_Eta_%2.1f_%2.1f;TruthE;RecoE/TruthE",Etabins[i],Etabins[i+1]), Ebins.size()-1, &Ebins[0], Rbins.size()-1, &Rbins[0]);
		Truth_vs_Truth[i] = new TProfile(Form("Truth_vs_Truth_Eta_%2.1f_%2.1f",Etabins[i],Etabins[i+1]),Form("Truth_vs_Truth_Eta_%2.1f_%2.1f;TruthE;TruthE",Etabins[i],Etabins[i+1]), Ebins.size()-1, &Ebins[0]);
		Truth_vs_Reco[i] = new TProfile(Form("Truth_vs_Reco_Eta_%2.1f_%2.1f",Etabins[i],Etabins[i+1]),Form("Truth_vs_Reco_Eta_%2.1f_%2.1f;TruthE;RecoE",Etabins[i],Etabins[i+1]), Ebins.size()-1, &Ebins[0]);
		Reco_vs_Reco[i] = new TProfile(Form("Reco_vs_Reco_Eta_%2.1f_%2.1f",Etabins[i],Etabins[i+1]),Form("Reco_vs_Reco_Eta_%2.1f_%2.1f;TruthE;RecoE",Etabins[i],Etabins[i+1]), Ebins.size()-1, &Ebins[0]);
		R_vs_Reco[i] = new TH2F(Form("R_vs_Reco_Eta_%2.1f_%2.1f",Etabins[i],Etabins[i+1]),Form("R_vs_Reco_Eta_%2.1f_%2.1f;RecoE;RecoE/TruthE",Etabins[i],Etabins[i+1]), Ebins.size()-1, &Ebins[0], Rbins.size()-1, &Rbins[0]);

		R_vs_Truth[i]->Sumw2();
		R_vs_Reco[i]->Sumw2();
		Truth_vs_Truth[i]->Sumw2();
		Truth_vs_Reco[i]->Sumw2();
		Reco_vs_Reco[i]->Sumw2();
	}

   Long64_t nentries = fChain->GetEntriesFast();
  // Long64_t nentries = 1e3;
			TProfile* JES = new TProfile("JES","",50,0,500,-2,2);
			TProfile* Mu = new TProfile("Mu","",50,0,500,-2,2);
			TProfile* El = new TProfile("El","",50,0,500,-2,2);	

			std::vector<TGraphErrors*> GraphVect;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
	if(jentry%10000 == 0) std::cout<< jentry <<"/"<<nentries<<std::endl;
      // if (Cut(ientry) < 0) continue;
		if(Reco_Jet_Size == 0 )continue;
		float evtweight = pu_weight*lumi_weight;
		for( int i = 0; i<Reco_Jet_Size; i++)
		{
			if(Reco_Jet_Pt->at(i)*1e-3 < 20 ) continue;
			if(Reco_Jet_DL1_bin->at(i) != 5) continue; // 5 : 60%, 4+5 : 70%, 3+4+5 : 77%, 2+3+4+5 : 85%
			TLorentzVector jet;
			jet.SetPtEtaPhiE(Reco_Jet_Pt->at(i), Reco_Jet_Eta->at(i), Reco_Jet_Phi->at(i), Reco_Jet_E->at(i));

			std::vector<TLorentzVector> truthmatch = matchTruth(jet);
			if(truthmatch.size() == 0)continue;
			TLorentzVector truth = truthmatch[0];

			TLorentzVector mu_in_jet = addMuon(jet);
			
			//TLorentzVector el_in_jet = addElectron(mu_in_jet);

			//JES->Fill(truth.Pt()*1e-3, truth.Pt()/jet.Pt());
			//Mu->Fill(truth.Pt()*1e-3, truth.Pt()/mu_in_jet.Pt());
			//El->Fill(truth.Pt()*1e-3, truth.Pt()/el_in_jet.Pt());


			TLorentzVector reco = mu_in_jet;
			
			int bin = hitsreco->FindBin(log(reco.Pt()*1e-3));
			if(bin == 0) bin = 1;
			
			ptreco_histos_lowEta[bin-1]->Fill(truth.Pt()/reco.Pt());
			Reco_vs_Truth_Inc->Fill(log(reco.Pt()*1e-3),log(truth.Pt()*1e-3));
			if(mu_in_jet.Pt() == jet.Pt())
			{
				
				int bin = hitsreco->FindBin(log(reco.Pt()*1e-3));
				if(bin == 0) bin = 1;

				ptreco_histos_had_lowEta[bin-1]->Fill(truth.Pt()/reco.Pt());
				Reco_vs_Truth_Had->Fill(log(reco.Pt()*1e-3),log(truth.Pt()*1e-3));
			}

			if(mu_in_jet.Pt() != jet.Pt())
			{

				int bin = hitsreco->FindBin(log(reco.Pt()*1e-3));
				if(bin == 0) bin = 1;
				ptreco_histos_semi_lowEta[bin-1]->Fill(truth.Pt()/reco.Pt());
				Reco_vs_Truth_Semi->Fill(log(reco.Pt()*1e-3),log(truth.Pt()*1e-3));


			}

		}
   }

for(int i = 0; i<nbins; i++)
{
	
	hitsreco->SetBinContent(i+1,ptreco_histos_lowEta[i]->GetMean());
	hitsreco->SetBinError(i+1,ptreco_histos_lowEta[i]->GetMeanError());

	hitsreco_semi->SetBinContent(i+1,ptreco_histos_semi_lowEta[i]->GetMean());
	hitsreco_semi->SetBinError(i+1,ptreco_histos_semi_lowEta[i]->GetMeanError());

	hitsreco_had->SetBinContent(i+1,ptreco_histos_had_lowEta[i]->GetMean());
	hitsreco_had->SetBinError(i+1,ptreco_histos_had_lowEta[i]->GetMeanError());

/*
	Double_t q = 0.5;
	Double_t median_inc, median_semi, median_had;

	ptreco_histos_lowEta[i]->GetQuantiles(1, &median_inc, &q);
	ptreco_histos_semi_lowEta[i]->GetQuantiles(1, &median_semi, &q);
	ptreco_histos_had_lowEta[i]->GetQuantiles(1, &median_had, &q);


	hitsreco_median->SetBinContent(i+1,median_inc);
	hitsreco_median->SetBinError(i+1,ptreco_histos_lowEta[i]->GetMeanError());

	hitsreco_semi_median->SetBinContent(i+1,median_semi);
	hitsreco_semi_median->SetBinError(i+1,ptreco_histos_semi_lowEta[i]->GetMeanError());

	hitsreco_had_median->SetBinContent(i+1,median_had);
	hitsreco_had_median->SetBinError(i+1,ptreco_histos_had_lowEta[i]->GetMeanError());
		
*/	
}


hitsreco->SetDirectory(file);
hitsreco_semi->SetDirectory(file);
hitsreco_had->SetDirectory(file);

/*
hitsreco_median->SetDirectory(file);
hitsreco_semi_median->SetDirectory(file);
hitsreco_had_median->SetDirectory(file);
*/
/*
Reco_vs_Truth_Semi->SetDirectory(file);
Reco_vs_Truth_Had->SetDirectory(file);
Reco_vs_Truth_Inc->SetDirectory(file);
*/

file->Write();
}

void Output :: Loop2(std::vector<TF1*> corrfunc)
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 4e6;
		std::vector<TGraphErrors*> RecoGraphVect;
   Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
		if(Reco_Jet_Size == 0 )continue;
		float evtweight = pu_weight*lumi_weight;
		for( int i = 0; i<Reco_Jet_Size; i++)
		{
			TLorentzVector jet;
			jet.SetPtEtaPhiE(Reco_Jet_Pt->at(i), Reco_Jet_Eta->at(i), Reco_Jet_Phi->at(i), Reco_Jet_E->at(i));

			std::vector<TLorentzVector> truthmatch = matchTruth(jet);
			if(truthmatch.size() == 0)continue;
			TLorentzVector truth = truthmatch[0];

			TLorentzVector mu_in_jet = addMuon(jet);
			
			TLorentzVector el_in_jet = addElectron(mu_in_jet);
		
			if(mu_in_jet.Pt() == jet.Pt())continue;

			TLorentzVector reco = mu_in_jet;

			TLorentzVector recoest = truth*getTruthFactor(truth.E(),reco.Eta(),corrfunc);
		
			FillRecoHistos(truth,reco,recoest);

		}
   }

	for(int i = 0; i<NEta; i++)
	{
		RecoGraphVect.push_back(constructGraph(R_vs_Reco[i],Reco_vs_Reco[i],true));
	}

	std::vector<TF1*> Reco_fitfunctions = FitGraphs(RecoGraphVect);

	for(int i = 0; i<Reco_fitfunctions.size(); i++)
	{
		std::cout<<Reco_fitfunctions[i]->GetName()<<std::endl;
		for(int j = 0; j<Reco_fitfunctions[i]->GetNpar(); j++)
		{
				
			std::cout<<Form("a%d",j)<< " = " << Reco_fitfunctions[i]->GetParameter(j) <<";"<< std::endl;
		}
	}

	std::cout<<"Loop3"<<std::endl;
	Loop3(Reco_fitfunctions);
	
}

void Output :: Loop3(std::vector<TF1*> corrfunc)
{

	  if (fChain == 0) return;


			TProfile* JES = new TProfile("JES_Loop3","",15,0,300,-2,2);
			TProfile* Mu = new TProfile("Mu_Loop3","",15,0,300,-2,2);
			TProfile* PtReco = new TProfile("PtReco_Loop3","",15,0,300,-2,2);


   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 4e3;

   Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
		if(Reco_Jet_Size == 0 )continue;
		float evtweight = pu_weight*lumi_weight;
		for( int i = 0; i<Reco_Jet_Size; i++)
		{
			TLorentzVector jet;
			jet.SetPtEtaPhiE(Reco_Jet_Pt->at(i), Reco_Jet_Eta->at(i), Reco_Jet_Phi->at(i), Reco_Jet_E->at(i));

			std::vector<TLorentzVector> truthmatch = matchTruth(jet);
			if(truthmatch.size() == 0)continue;
			TLorentzVector truth = truthmatch[0];

			TLorentzVector mu_in_jet = addMuon(jet);
			
			TLorentzVector el_in_jet = addElectron(mu_in_jet);
		
			if(mu_in_jet.Pt() == jet.Pt())continue;

			TLorentzVector reco = mu_in_jet;

			TLorentzVector recoest = reco*getRecoFactor(reco.E(),reco.Eta(),corrfunc);
		
			JES->Fill(truth.Pt()*1e-3,jet.Pt()/truth.Pt());
			Mu->Fill(truth.Pt()*1e-3,reco.Pt()/truth.Pt());
			PtReco->Fill(truth.Pt()*1e-3,recoest.Pt()/truth.Pt());

		}
   }
 /*
JES->SetDirectory(file);
Mu->SetDirectory(file);
PtReco->SetDirectory(file);
*/
file->Write();

}
double Output :: getRecoFactor(double E, double eta, std::vector<TF1*> corrfunc)
{
	double factor = 1.;
	TString FuncName;
	int bin;
	for(int i = 0; i<NEta; i++)
	{
		double etamin = Etabins[i];
		double etamax = Etabins[i+1];
		
		if( !(etamin <= eta && eta < etamax) ) continue;
		FuncName = Form("Fit_Graph_R_vs_Reco_Eta_%2.1f_%2.1f_NPar",etamin,etamax);
		bin = i;
	}	
	for(int i = 0; i<corrfunc.size(); i++)
	{
		TString Name;
		Name = corrfunc[i]->GetName();
		//std::cout<<Name<< " ** "<<FuncName<<std::endl; 
		if(!(Name.Contains(FuncName))) continue;
			if(E >= 400000){ factor = corrfunc[i]->Eval(400000);
			}else{
			factor = corrfunc[i]->Eval(E);
			}
			//std::cout<<E<<" "<<eta<<" "<<1./(factor)<<std::endl;
			break;
	}
	
	return 1./(factor);
}

void Output :: constructHisto(std::vector<TGraphErrors*> graphs)
{

	for(int i = 0; i<graphs.size(); i++)
	{
		TH1F* hist = graphs[i]->GetHistogram();
		hist->SetNameTitle(Form("Hist_%s_",graphs[i]->GetName()),Form("Hist_%s_",graphs[i]->GetName()));
		int n = graphs[i]->GetN();
		for( int ii = 0; ii<n; ii++)
		{
			double x,y;
			graphs[i]->GetPoint(ii,x,y);
			hist->SetBinContent(hist->FindBin(x),y);
		}
		hist->SetDirectory(file);
	}
}
double Output :: getTruthFactor(double E, double eta, std::vector<TF1*> corrfunc)
{
	double factor = 1.;
	TString FuncName;
	int bin;
	for(int i = 0; i<NEta; i++)
	{
		double etamin = Etabins[i];
		double etamax = Etabins[i+1];
		
		if( !(etamin <= eta && eta < etamax) ) continue;
		FuncName = Form("Fit_Graph_R_vs_Truth_Eta_%2.1f_%2.1f_NPar",etamin,etamax);
		bin = i;
	}	
	for(int i = 0; i<corrfunc.size(); i++)
	{
		TString Name;
		Name = corrfunc[i]->GetName();
		//std::cout<<Name<< " ** "<<FuncName<<std::endl; 
		if(!(Name.Contains(FuncName))) continue;
			if(E >= 400000){ factor = corrfunc[i]->Eval(400000);
			}else{
			factor = corrfunc[i]->Eval(E);
			}
			//std::cout<<E<<" "<<eta<<" "<<factor<<std::endl;
			break;
	}
	
	return factor;
}

TF1* Output :: FitGraph(TGraphErrors* g, int npar)
{
	_npar = npar;
	double xmax,ymax,ymin,xmin;
	int npoint = g->GetN();
	g->GetPoint(1,xmin,ymin);
	g->GetPoint(npoint-1,xmax,ymax);
	TF1 *func = new TF1(Form("Function_%s_",g->GetName()),FitFunction,xmin,xmax,_npar);
	g->Fit(Form("Function_%s_",g->GetName()),"RQMEB");
	return func;
}

std::vector<TF1*> Output :: FitGraphs(std::vector<TGraphErrors*> graphs)
{
	
	int NGraph = graphs.size();
	std::vector<TF1*> result;
	for( int i = 0; i<NGraph; i++)
	{
		double xmax,ymax,ymin,xmin;
		int npoint = graphs[i]->GetN();
		graphs[i]->GetPoint(0,xmin,ymin);
		graphs[i]->GetPoint(npoint-1,xmax,ymax);
		std::vector<TString> functionlist = genFunctionList();
		double bestChi2 = 10000;
		int func_j = 100;
		for(int j = 0; j<functionlist.size(); j++)
		{
			//std::cout<<"Fitting with ... "<<functionlist[j]<<std::endl;
			TF1 *func = new TF1(Form("Fun_%d_Graph_%d",j,i),functionlist[j],xmin,xmax);
			graphs[i]->Fit(func,"RQMB");
			double chi2 = (func->GetChisquare()/(Double_t)func->GetNDF());
			//std::cout<<chi2<<std::endl;
			if(chi2 < bestChi2)
			{
				bestChi2 = chi2;
				func_j = j;
			}
		}
		
		TF1* BestFitFunction = new TF1(Form("Fit_%s_NPar_%d",graphs[i]->GetName(),func_j+2),functionlist[func_j],xmin,xmax);
		graphs[i]->Fit(BestFitFunction,"RQMEB");
		result.push_back(BestFitFunction);

		/*
		int npar = 6;
		TF1* FirstFitFunction = FitGraph(graphs[i],2);
		double bestChi2Ndf = (FirstFitFunction->GetChisquare()/FirstFitFunction->GetNDF());
		std::cout<<"Graph : "<<i+1<<std::endl;
		std::cout<<"BestChi2NDf : " << bestChi2Ndf <<std::endl;
		int BestNpar = 2;
		delete FirstFitFunction;
		for( int par = 3; par<=npar; par++)
		{
			TF1* FitFunction = FitGraph(graphs[i],par);
			double Chi2Ndf = (FitFunction->GetChisquare()/FitFunction->GetNDF());
			delete FitFunction;
			if(Chi2Ndf > bestChi2Ndf) continue;
				bestChi2Ndf = Chi2Ndf;
				BestNpar = par;
		}
		std::cout<<"BestChi2NDf : " << bestChi2Ndf << " BestNPar "<< BestNpar <<std::endl;
		TF1* BestFitFunction = FitGraph(graphs[i],BestNpar);
		BestFitFunction->SetNameTitle(Form("Fit_%s_NPar_%1d",graphs[i]->GetName(),BestNpar),Form("Fit_%s_NPar_%1d",graphs[i]->GetName(),BestNpar));
		result.push_back(BestFitFunction);
		*/

		file->cd();
		graphs[i]->Write();
		//delete BestFitFunction;
	}
	return result;
}

double Output :: getEta(TH2F *hist,bool doReco)
{
	double eta = 0;
	for(int i = 0; i<NEta; i++)
	{
		double etamin = Etabins[i];
		double etamax = Etabins[i+1];
		TString string = Form("R_vs_Truth_Eta_%2.1f_%2.1f",Etabins[i],Etabins[i+1]);
		if(doReco) string = Form("R_vs_Reco_Eta_%2.1f_%2.1f",Etabins[i],Etabins[i+1]);
		if(string != hist->GetName()) continue;
		eta = (etamin+etamax)/2;	
	}
	return eta;
}
std::vector<int> Output :: getBinList(TH2F *hist,bool doReco)
{
	int NBins = hist->GetNbinsX();	
	std::vector<int> list;
	int FirstBin = 1;
  
	for(int i = FirstBin; i<=NBins; i++)
	{
		float val = -3.1*exp(1.05*TMath::Abs(getEta(hist,doReco)))+400;
		val = 0;
		if(getNEff(project1D(hist, FirstBin, i))>val)
		{
		list.push_back(FirstBin);
		list.push_back(i);
		FirstBin = i+1;
		}
	}
	return list;

}

double Output :: getNEff(TH1 *h) 
{
  if (h->GetSumw2N()==0) return 0;

  int N=h->GetNbinsX(); double w=0, w2=0;
  for (int ibin=1;ibin<=N;++ibin) {
    w  += h->GetBinContent(ibin);
    w2 += h->GetSumw2()->At(ibin);
  }

  return w2>0?w*w/w2:0;
}

void Output :: getAvg(TProfile *inputProfile, int low, int high, double &avg, double &avg_err)
{
  //Calc average
  double sum_wy=0, sum_w=0;
  for (int i=low; i<=high;++i) {
    double y=inputProfile->GetBinContent(i), err=inputProfile->GetBinError(i);
    if (err==0) continue;
    sum_w+=1.0/err/err; sum_wy+=y/err/err;
 }
  avg = sum_w>0?sum_wy/sum_w:0;

  //Calc error
  double sum_w_err=0;
  for (int i=low; i<=high;++i) {
    double err=inputProfile->GetBinError(i);
    if (err>0) sum_w_err+=1.0/err/err;
  }
  avg_err = sum_w_err>0?1.0/sqrt(sum_w_err):0;
}

TH1D* Output :: project1D(TH2F *h, int low_bin, int high_bin)
{
  return h->ProjectionY(Form("%sBin%d-%d",h->GetName(),low_bin,high_bin),
			low_bin,high_bin);
}

void Output :: FitGauss(TH1D *R, TF1 *_gauss, Double_t fitMin, Double_t fitMax, bool adjustRange)
{

  _gauss->SetRange(fitMin,fitMax);
  _gauss->SetParameter(1, R->GetBinCenter(R->GetMaximumBin()));
  R->Fit(_gauss,"REMQ");// MAX changed LREQ here

  if (adjustRange){
    double mean=_gauss->GetParameter(1), w=(fitMax-fitMin)/2;
    if (w>0.3) w=0.3;
    if (mean<fitMin) mean=fitMin;
    _gauss->SetRange(mean-w,mean+w);
    R->Fit(_gauss,"REMQ");
  }

}
double Output :: getFitMax(TH1D *inputHist, double gausFitNrms)
{
  double tmpRMS = inputHist->GetRMS(), tmpPeak = inputHist->GetBinCenter(inputHist->GetMaximumBin());
  return tmpPeak + gausFitNrms*tmpRMS;
}

double Output :: getFitMin(TH1D *inputHist, double gausFitNrms)
{
  double tmpRMS = inputHist->GetRMS(), tmpPeak = inputHist->GetBinCenter(inputHist->GetMaximumBin());
  double min = tmpPeak - gausFitNrms*tmpRMS;

  int binMin = inputHist->FindBin(min);
  double minContent = inputHist->GetBinContent(binMin);

  if(tmpPeak < 0.6)
    min = 0.25;
  return min;
}

TGraphErrors* Output :: constructGraph(TH2F* _2DHistos, TProfile* _AvrgE, bool doReco=false)
{

		TGraphErrors *Graph = new TGraphErrors();
		double EAvg, EAvgErr;
		std::vector<int> BinList = getBinList(_2DHistos,doReco);
		for (Int_t j=0; j<BinList.size(); j=j+2){
    			int lowBin=BinList[j] , highBin=BinList[j+1];
    			getAvg(_AvrgE,lowBin,highBin,EAvg,EAvgErr);
    			TH1D *tmpHist = project1D(_2DHistos,lowBin,highBin);
			TF1  *tmpHistFit = new TF1("tmpHistFit","gaus",0.,5.);
			FitGauss(tmpHist,tmpHistFit,getFitMin(tmpHist,1.8),getFitMax(tmpHist,1.8), false);

			Int_t pointN = Graph->GetN();
			
      			Graph->SetPoint(pointN,EAvg,tmpHistFit->GetParameter(1));
     			if(pointN > 1){Graph->SetPointError(pointN,EAvgErr,tmpHistFit->GetParError(1));}
			/*
			Graph->SetPoint(pointN,EAvg,tmpHist->GetMean());
     			Graph->SetPointError(pointN,EAvgErr,tmpHist->GetMeanError(1));
			*/
			delete tmpHistFit;
			tmpHist->SetDirectory(file);
		}
	Graph->SetNameTitle(Form("Graph_%s",_2DHistos->GetName()),Form("Graph_%s",_2DHistos->GetName()));	

	return Graph;
}

void Output :: FillHistos(TLorentzVector truth, TLorentzVector reco, double weight = 1)
{
	for(int i = 0; i<NEta; i++)
	{
		double etamin = Etabins[i];
		double etamax = Etabins[i+1];
		double eta = reco.Eta();
		
		if( !(etamin <= eta && eta < etamax) ) continue;
		
		R_vs_Truth[i]->Fill(truth.E(),reco.E()/truth.E(),weight);
		Truth_vs_Truth[i]->Fill(truth.E(),truth.E(),weight);
		Truth_vs_Reco[i]->Fill(truth.E(),reco.E(),weight);

	}
}
void Output :: FillRecoHistos(TLorentzVector truth, TLorentzVector reco, TLorentzVector recoest, double weight = 1)
{
	for(int i = 0; i<NEta; i++)
	{
		double etamin = Etabins[i];
		double etamax = Etabins[i+1];
		double eta = reco.Eta();
		
		if( !(etamin <= eta && eta < etamax) ) continue;
		
		R_vs_Reco[i]->Fill(truth.E(),reco.E()/truth.E(),weight);
		Reco_vs_Reco[i]->Fill(truth.E(),recoest.E(),weight);

	}
}

std::vector<TLorentzVector> Output :: matchTruth(TLorentzVector jet)
{
	std::vector<TLorentzVector> result;
	std::vector<TLorentzVector> out;
	
	for( int i = 0; i<Truth_WZJet_Size; i++)
	{
		if(Truth_WZJet_ConeTruthLabelID->at(i) != 5) continue;
		TLorentzVector truthjet;
		truthjet.SetPtEtaPhiE(Truth_WZJet_Pt->at(i), Truth_WZJet_Eta->at(i), Truth_WZJet_Phi->at(i), Truth_WZJet_E->at(i));
		
		double dr = truthjet.DeltaR(jet);
		if( dr > 0.3 ) continue;
		result.push_back(truthjet);
	}
	if( result.size() >= 1)
	{
		
		double dr0 = result[0].DeltaR(jet);
		int n = 0;
		for( int j = 0; j<result.size(); j++)
		{
			double dr = result[j].DeltaR(jet);
			if(dr > dr0)continue;
				dr0 = dr;
				n = j;
		}

		out.push_back(result[n]);
		return out;
	}

	return result;

}

TLorentzVector Output :: addMuon(TLorentzVector jet)
{
	std::vector<TLorentzVector> muon_in_jet;
	std::vector<TLorentzVector> loss_in_jet;		
	for( int i = 0; i<Reco_Muon_Size; i++)
	{
		TLorentzVector muon,loss;
		muon.SetPtEtaPhiE(Reco_Muon_Pt->at(i), Reco_Muon_Eta->at(i), Reco_Muon_Phi->at(i), Reco_Muon_E->at(i));
		loss.SetPtEtaPhiE(Reco_Muon_ELoss_Pt->at(i), Reco_Muon_ELoss_Eta->at(i), Reco_Muon_ELoss_Phi->at(i), Reco_Muon_ELoss_E->at(i));
		double dr = muon.DeltaR(jet);
		double mindr = std::min(0.4,0.04+(10./(muon.Pt()*1e-3)));
		if( dr > mindr)continue;
		muon_in_jet.push_back(muon);
		loss_in_jet.push_back(loss);
	}
	
	if(muon_in_jet.size() == 0){return jet;}
	if(muon_in_jet.size() == 1){return jet-loss_in_jet[0]+muon_in_jet[0];}

	if(muon_in_jet.size() > 1)
	{
		double pt0 = muon_in_jet[0].Pt();
		int n = 0;
		for(int j = 0; j<muon_in_jet.size(); j++)
		{
			if(pt0 > muon_in_jet[j].Pt()) continue;
				pt0 = muon_in_jet[j].Pt();
				n = j;
		}
		return jet-loss_in_jet[n]+muon_in_jet[n];
	}

	return jet;
}

TLorentzVector Output :: addElectron(TLorentzVector jet)
{
	std::vector<TLorentzVector> electron_in_jet;		
	for( int i = 0; i<Reco_Electron_Size; i++)
	{
		TLorentzVector electron;
		electron.SetPtEtaPhiE(Reco_Electron_Pt->at(i), Reco_Electron_Eta->at(i), Reco_Electron_Phi->at(i), Reco_Electron_E->at(i));
		double dr = electron.DeltaR(jet);
		double mindr = std::min(0.4,0.04+(10./(electron.Pt()*1e-3)));
		if( dr > mindr)continue;
		electron_in_jet.push_back(electron);
	}
	
	if(electron_in_jet.size() == 0){return jet;}
	if(electron_in_jet.size() == 1){return jet+electron_in_jet[0];}

	if(electron_in_jet.size() > 1)
	{
		double pt0 = electron_in_jet[0].Pt();
		int n = 0;
		for(int j = 0; j<electron_in_jet.size(); j++)
		{
			if(pt0 > electron_in_jet[j].Pt()) continue;
				pt0 = electron_in_jet[j].Pt();
				n = j;
		}
		return jet+electron_in_jet[n];
	}

	return jet;
}

std::vector<TString> Output::genFunctionList()
{
	std::vector<TString> list;
	
	TString func = "[0]+[1]*log(x)+[2]*log(x)^2";
	list.push_back(func);

	for( int i = 3; i<=10; i++)
	{
		func.Append(Form("+[%d]*log(x)^%d",i,i));
		list.push_back(func);
	
	}	

	return list;
}
