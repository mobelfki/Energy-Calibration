#include "Var.h"
#include "TSystemDirectory.h"

int NTUP()
{

	return GetNTuple(0,0,"");
}

int GetNTuple(int i, int n, TString out){
	SetOuput(out,n);
	TChain *chain = Load();
	SetBranches(chain);
	Long64_t nentries = chain->GetEntriesFast();
	_Thrust = GenerateT();
	std::cout<<"Total Number of Entries : "<<nentries<<std::endl;
	if( n == 0)
	{
		n=nentries;
	}
	std::cout<<"Event Will be processed : "<< n <<std::endl;
	for(unsigned int jentry = i; jentry<n; jentry++)
	{
		if(jentry%10000 == 0) std::cout<<" Events : "<<jentry<<"/"<<n<<std::endl;
		Int_t nb = chain->GetEntry(jentry);
		
		SetLumi(chain);	
		NewEvents();
		ClearQuantities();

		if(!Ana())
		{
			NewEvents();
			ClearQuantities();
			continue;
		}
		
		ComputeQuantities();
		_output1_tree->Fill();
		NewEvents();
		ClearQuantities();
		
		//std::cout<<"Loop finished"<<std::endl;
	}
	std::cout<<"Loop finished"<<std::endl;

	_output1_tree->AutoSave();
	_output_file->Close();
	return 1.;
}

TChain* Load()
{
	std::vector<TString> Samples;
	Samples.push_back("/eos/user/m/mobelfki/TestSample/mc16a/mc16a.aMCnloHwpp_hh_yybb_AF2.MxAODDetailed.e4419_a875_r9364_p3629.h024.root/sample.root");
	Samples.push_back("/eos/user/m/mobelfki/TestSample/mc16d/mc16d.aMCnloHwpp_hh_yybb_AF2.MxAODDetailed.e4419_a875_r10201_p3629.h024.root/sample_01.root");
	Samples.push_back("/eos/user/m/mobelfki/TestSample/mc16d/mc16d.aMCnloHwpp_hh_yybb_AF2.MxAODDetailed.e4419_a875_r10201_p3629.h024.root/sample_02.root");
	Samples.push_back("/eos/user/m/mobelfki/TestSample/mc16e/mc16e.aMCnloHwpp_hh_yybb_AF2.MxAODDetailed.e4419_a875_r10724_p3705.h024.root/sample.root");

	std::cout<<"Start Loading the Trees"<<std::endl;
	TChain *chain = GetChains(Samples, "CollectionTree");
	return chain;
}

TChain* GetChains(std::vector<TString> Samples, TString TreeName)
{
	TChain *chain = new TChain(TreeName);
	std::cout<<"Samples To Read "<<std::endl;
	for (unsigned int i = 0; i<Samples.size(); i++)
	{
		std::cout<<Samples[i]<<std::endl;
		chain->AddFile(Samples[i],-1,TreeName);
		
	}
	std::cout<<"Trees are loaded"<<std::endl;
	return chain;
}

bool Ana()
{
	
	if(!(_event_ispass)) return false;
	if( _event_N_lep > 0) return false;
	if( _event_N_j_central > 6) return false;

	if(debug)std::cout<<"Event Passed"<<std::endl;
	for(unsigned int i=0; i<_gamma_pt->size(); i++)
	{
		TLorentzVector gamma;
		gamma.SetPtEtaPhiM(_gamma_pt->at(i),_gamma_eta->at(i),_gamma_phi->at(i),0.0);
		_photons.push_back(gamma);
	}
	if(_photons.size() != 2) return false;
	
	_Hyy = _photons.at(0) + _photons.at(1);

	if(_Hyy.M()*0.001 < 105 && _Hyy.M()*0.001 > 160) return false;

	if(debug)std::cout<<"Two Photons are selected"<<std::endl;
	vector<TLorentzVector> jes_bjets, jes_nbjets;
	vector<TLorentzVector> muj_bjets, muj_nbjets;

	vector<Int_t> bjets_nmu;
	vector<Int_t> bjets_rank;
	vector<Int_t> bjets_DL1_rank;

	vector<Float_t> bjets_score;

//	std::cout<<_jes_jets_pt->size() - _muj_jets_pt->size()<<std::endl;

	for(unsigned int i=0; i<_jes_jets_pt->size(); i++)
	{

		
		Int_t rank;
		Int_t rank_DL1;
		TLorentzVector jes_bjet, jes_nbjet;
		TLorentzVector muj_bjet, muj_nbjet;
		if(_jets_btag_60->at(i))
		{
			rank = 1;
		}else if(_jets_btag_70->at(i))
		{
			rank = 2;
		}else if(_jets_btag_77->at(i))
		{
			rank = 3;
		}else if(_jets_btag_85->at(i))
		{
			rank = 4;
		}else {
			rank = 0;
		}

		if(_jets_DL1_btag_60->at(i))
		{
			rank_DL1 = 1;
		}else if(_jets_DL1_btag_70->at(i))
		{
			rank_DL1 = 2;
		}else if(_jets_DL1_btag_77->at(i))
		{
			rank_DL1 = 3;
		}else if(_jets_DL1_btag_85->at(i))
		{
			rank_DL1 = 4;
		}else {
			rank_DL1 = 0;
		}

		_jet_rank.push_back(rank);



	//	_jet_DL1_rank.push_back(rank_DL1);

		if(rank != 0 and rank_DL1 != 0)
		{
			jes_bjet.SetPtEtaPhiE(_jes_jets_pt->at(i),_jes_jets_eta->at(i),_jes_jets_phi->at(i),_jes_jets_m->at(i));
			muj_bjet.SetPtEtaPhiE(_muj_jets_pt->at(i),_muj_jets_eta->at(i),_muj_jets_phi->at(i),_muj_jets_m->at(i));
			jes_bjets.push_back(jes_bjet);
			muj_bjets.push_back(muj_bjet);
			bjets_nmu.push_back(_jets_nmu->at(i));
			bjets_DL1_rank.push_back(rank_DL1);
			bjets_rank.push_back(rank);
			bjets_score.push_back(_jets_score->at(i));
		}else{
			jes_nbjet.SetPtEtaPhiE(_jes_jets_pt->at(i),_jes_jets_eta->at(i),_jes_jets_phi->at(i),_jes_jets_m->at(i));
			muj_nbjet.SetPtEtaPhiE(_muj_jets_pt->at(i),_muj_jets_eta->at(i),_muj_jets_phi->at(i),_muj_jets_m->at(i));
			jes_nbjets.push_back(jes_nbjet);
			muj_nbjets.push_back(muj_nbjet);
		}
		
	}
	if(debug)std::cout<<"N B Jet : "<<jes_bjets.size()<< " N Jet "<< jes_bjets.size()+jes_nbjets.size() << " NJet "<< _event_njets <<std::endl;
	if(!(jes_bjets.size() >= 2)) return false;


	// sort jet by pT
	
	for(unsigned int j=1; j<=jes_bjets.size(); j++)
	{
		for(unsigned int i=0; i<jes_bjets.size()-1; i++)
		{
			if(jes_bjets.at(i).Pt() < jes_bjets.at(i+1).Pt())
			{
				float c;
				int nmu, rank, rank_DL1;
				c    = bjets_score.at(i);
				nmu  = bjets_nmu.at(i);
				rank = bjets_rank.at(i);
				rank_DL1 = bjets_DL1_rank.at(i);
				TLorentzVector jet(jes_bjets.at(i));

				bjets_score.at(i) = bjets_score.at(i+1);
				jes_bjets.at(i)   = jes_bjets.at(i+1);
				bjets_nmu.at(i)   = bjets_nmu.at(i+1);
				bjets_rank.at(i)  = bjets_rank.at(i+1);
				bjets_DL1_rank.at(i)  = bjets_DL1_rank.at(i+1);

				bjets_score.at(i+1) = c;
				bjets_nmu.at(i+1)   = nmu;
				bjets_rank.at(i+1)  = rank;
				bjets_DL1_rank.at(i+1)  = rank_DL1;
				jes_bjets.at(i+1).SetPtEtaPhiE(jet.Pt(),jet.Eta(),jet.Phi(),jet.E());
			}
		} 
	}
	

	_Bjets.push_back(jes_bjets.at(0));
	_Bjets.push_back(jes_bjets.at(1));
	
	_Bjets_muj.push_back(muj_bjets.at(0));
	_Bjets_muj.push_back(muj_bjets.at(1));

	_score_b1 = bjets_score.at(0);
	_score_b2 = bjets_score.at(1);

	_nmu_b1   = bjets_nmu.at(0);
	_nmu_b2   = bjets_nmu.at(1);

	_rank_b1   = bjets_rank.at(0);
	_rank_b2   = bjets_rank.at(1);

	_rank_DL1_b1   = bjets_DL1_rank.at(0);
	_rank_DL1_b2   = bjets_DL1_rank.at(1);

	for(unsigned int i=2; i<jes_bjets.size(); i++){_jets.push_back(jes_bjets.at(i));}
	_n_add_bjet = _jets.size();
	for(unsigned int i=0; i<jes_nbjets.size(); i++){_jets.push_back(jes_nbjets.at(i));}

	if(debug){
	for(unsigned int i=0; i<bjets_score.size(); i++)
	{
		std::cout<<bjets_score.at(i)<<std::endl;
	}
	}
	if(debug)std::cout<<"Two BJets are selected"<<std::endl;

	_Hbb     = _Bjets.at(0) + _Bjets.at(1);
	_Hbb_muj = _Bjets_muj.at(0) + _Bjets_muj.at(1);

	_HH  = _Hyy + _Hbb;
	for(int i=0; i<_jets.size(); i++) {_SumJ += _jets.at(i);}

	return true;
}

void ComputeQuantities()
{
	
	_pt_b1  = _Bjets.at(0).Pt();
	_pt_b2  = _Bjets.at(1).Pt();
	_eta_b1 = _Bjets.at(0).Eta();
	_eta_b2 = _Bjets.at(1).Eta();
	_phi_b1 = _Bjets.at(0).Phi();
	_phi_b2 = _Bjets.at(1).Phi();
	_e_b1   = _Bjets.at(0).E();
	_e_b2   = _Bjets.at(1).E();

	//std::cout<< _Bjets_muj.size()<<std::endl;
	//std::cout<< _Bjets_muj.at(0).Pt()<<std::endl;

	_pt_b1_muj  = _Bjets_muj.at(0).Pt();
	_pt_b2_muj  = _Bjets_muj.at(1).Pt();
	_eta_b1_muj = _Bjets_muj.at(0).Eta();
	_eta_b2_muj = _Bjets_muj.at(1).Eta();
	_phi_b1_muj = _Bjets_muj.at(0).Phi();
	_phi_b2_muj = _Bjets_muj.at(1).Phi();
	_e_b1_muj   = _Bjets_muj.at(0).E();
	_e_b2_muj   = _Bjets_muj.at(1).E();

	_pt_y1  = _photons.at(0).Pt();
	_pt_y2  = _photons.at(1).Pt();
	_eta_y1 = _photons.at(0).Eta();
	_eta_y2 = _photons.at(1).Eta();
	_phi_y1 = _photons.at(0).Phi();
	_phi_y2 = _photons.at(1).Phi();
	_e_y1   = _photons.at(0).E(); 
	_e_y2   = _photons.at(1).E();

	vector<TLorentzVector> _yybb, _hh;
	_yybb.push_back(_Bjets.at(0));
	_yybb.push_back(_Bjets.at(1));
	_yybb.push_back(_photons.at(0));
	_yybb.push_back(_photons.at(1));

	_hh.push_back(_Hyy);
	_hh.push_back(_Hbb);
	

	TVector3 thrust = GetThrust(_Bjets);
	FillLMomenta(_Bjets,20,thrust);

	_m_yy     = _Hyy.M();
	_m_bb     = _Hbb.M();
	_m_bb_muj = _Hbb_muj.M();
	_m_hh     = _HH.M();

	_pt_hh   = _HH.Pt();
	_pt_yy   = _Hyy.Pt();
	_pt_bb   = _Hbb.Pt();

	_eta_hh  = _HH.Eta();
	_eta_yy  = _Hyy.Eta();
	_eta_bb  = _Hbb.Eta();

	_phi_hh  = _HH.Phi();
	_phi_bb  = _Hbb.Phi();
	_phi_yy  = _Hyy.Phi();
	
	_e_hh    = _HH.E();
	_e_yy    = _Hyy.E();
	_e_bb    = _Hbb.E();

	_dr_hh    = _Hbb.DeltaR(_Hyy);
	_dr_bb    = _Bjets.at(0).DeltaR(_Bjets.at(1));
	_dr_yy    = _photons.at(0).DeltaR(_photons.at(1));
	_dr_b1yy  = _Hyy.DeltaR(_Bjets.at(0));
	_dr_b2yy  = _Hyy.DeltaR(_Bjets.at(1));
	_dr_y1bb  = _Hbb.DeltaR(_photons.at(0));
	_dr_y2bb  = _Hbb.DeltaR(_photons.at(1));

	_pt_jets  = _SumJ.Pt();
	_eta_jets = _SumJ.Eta();
	_phi_jets = _SumJ.Phi();
	_m_jets   = _SumJ.M();

	_n_add_j   = _jets.size()-_n_add_bjet;
	
	_n_jets    = _event_njets;

	_n_btag    = _event_N_j_btag;

	FillFoxMomenta(_Bjets,20);

	_B_HT_l1 = FoxMomenta_l(_Bjets,'T',1);
	_B_HT_l2 = FoxMomenta_l(_Bjets,'T',2);
	_B_HT_l3 = FoxMomenta_l(_Bjets,'T',3);
	_B_HT_l4 = FoxMomenta_l(_Bjets,'T',4);
	_B_HT_l5 = FoxMomenta_l(_Bjets,'T',5);

	_B_HP_l1 = FoxMomenta_l(_Bjets,'P',1);
	_B_HP_l2 = FoxMomenta_l(_Bjets,'P',2);
	_B_HP_l3 = FoxMomenta_l(_Bjets,'P',3);
	_B_HP_l4 = FoxMomenta_l(_Bjets,'P',4);
	_B_HP_l5 = FoxMomenta_l(_Bjets,'P',5);

	_B_HS_l1 = FoxMomenta_l(_Bjets,'S',1);
	_B_HS_l2 = FoxMomenta_l(_Bjets,'S',2);
	_B_HS_l3 = FoxMomenta_l(_Bjets,'S',3);
	_B_HS_l4 = FoxMomenta_l(_Bjets,'S',4);
	_B_HS_l5 = FoxMomenta_l(_Bjets,'S',5);

	_B_HZ_l1 = FoxMomenta_l(_Bjets,'Z',1);
	_B_HZ_l2 = FoxMomenta_l(_Bjets,'Z',2);
	_B_HZ_l3 = FoxMomenta_l(_Bjets,'Z',3);
	_B_HZ_l4 = FoxMomenta_l(_Bjets,'Z',4);
	_B_HZ_l5 = FoxMomenta_l(_Bjets,'Z',5);

	_B_HY_l1 = FoxMomenta_l(_Bjets,'Y',1);
	_B_HY_l2 = FoxMomenta_l(_Bjets,'Y',2);
	_B_HY_l3 = FoxMomenta_l(_Bjets,'Y',3);
	_B_HY_l4 = FoxMomenta_l(_Bjets,'Y',4);
	_B_HY_l5 = FoxMomenta_l(_Bjets,'Y',5);

	_B_H1_l1 = FoxMomenta_l(_Bjets,'1',1);
	_B_H1_l2 = FoxMomenta_l(_Bjets,'1',2);
	_B_H1_l3 = FoxMomenta_l(_Bjets,'1',3);
	_B_H1_l4 = FoxMomenta_l(_Bjets,'1',4);
	_B_H1_l5 = FoxMomenta_l(_Bjets,'1',5);

//PHOTONS

	_Y_HT_l1 = FoxMomenta_l(_photons,'T',1);
	_Y_HT_l2 = FoxMomenta_l(_photons,'T',2);
	_Y_HT_l3 = FoxMomenta_l(_photons,'T',3);
	_Y_HT_l4 = FoxMomenta_l(_photons,'T',4);
	_Y_HT_l5 = FoxMomenta_l(_photons,'T',5);

	_Y_HP_l1 = FoxMomenta_l(_photons,'P',1);
	_Y_HP_l2 = FoxMomenta_l(_photons,'P',2);
	_Y_HP_l3 = FoxMomenta_l(_photons,'P',3);
	_Y_HP_l4 = FoxMomenta_l(_photons,'P',4);
	_Y_HP_l5 = FoxMomenta_l(_photons,'P',5);

	_Y_HS_l1 = FoxMomenta_l(_photons,'S',1);
	_Y_HS_l2 = FoxMomenta_l(_photons,'S',2);
	_Y_HS_l3 = FoxMomenta_l(_photons,'S',3);
	_Y_HS_l4 = FoxMomenta_l(_photons,'S',4);
	_Y_HS_l5 = FoxMomenta_l(_photons,'S',5);

	_Y_HZ_l1 = FoxMomenta_l(_photons,'Z',1);
	_Y_HZ_l2 = FoxMomenta_l(_photons,'Z',2);
	_Y_HZ_l3 = FoxMomenta_l(_photons,'Z',3);
	_Y_HZ_l4 = FoxMomenta_l(_photons,'Z',4);
	_Y_HZ_l5 = FoxMomenta_l(_photons,'Z',5);

	_Y_HY_l1 = FoxMomenta_l(_photons,'Y',1);
	_Y_HY_l2 = FoxMomenta_l(_photons,'Y',2);
	_Y_HY_l3 = FoxMomenta_l(_photons,'Y',3);
	_Y_HY_l4 = FoxMomenta_l(_photons,'Y',4);
	_Y_HY_l5 = FoxMomenta_l(_photons,'Y',5);

	_Y_H1_l1 = FoxMomenta_l(_photons,'1',1);
	_Y_H1_l2 = FoxMomenta_l(_photons,'1',2);
	_Y_H1_l3 = FoxMomenta_l(_photons,'1',3);
	_Y_H1_l4 = FoxMomenta_l(_photons,'1',4);
	_Y_H1_l5 = FoxMomenta_l(_photons,'1',5);

//YY+BB


	_YB_HT_l1 = FoxMomenta_l(_yybb,'T',1);
	_YB_HT_l2 = FoxMomenta_l(_yybb,'T',2);
	_YB_HT_l3 = FoxMomenta_l(_yybb,'T',3);
	_YB_HT_l4 = FoxMomenta_l(_yybb,'T',4);
	_YB_HT_l5 = FoxMomenta_l(_yybb,'T',5);

	_YB_HP_l1 = FoxMomenta_l(_yybb,'P',1);
	_YB_HP_l2 = FoxMomenta_l(_yybb,'P',2);
	_YB_HP_l3 = FoxMomenta_l(_yybb,'P',3);
	_YB_HP_l4 = FoxMomenta_l(_yybb,'P',4);
	_YB_HP_l5 = FoxMomenta_l(_yybb,'P',5);

	_YB_HS_l1 = FoxMomenta_l(_yybb,'S',1);
	_YB_HS_l2 = FoxMomenta_l(_yybb,'S',2);
	_YB_HS_l3 = FoxMomenta_l(_yybb,'S',3);
	_YB_HS_l4 = FoxMomenta_l(_yybb,'S',4);
	_YB_HS_l5 = FoxMomenta_l(_yybb,'S',5);

	_YB_HZ_l1 = FoxMomenta_l(_yybb,'Z',1);
	_YB_HZ_l2 = FoxMomenta_l(_yybb,'Z',2);
	_YB_HZ_l3 = FoxMomenta_l(_yybb,'Z',3);
	_YB_HZ_l4 = FoxMomenta_l(_yybb,'Z',4);
	_YB_HZ_l5 = FoxMomenta_l(_yybb,'Z',5);

	_YB_HY_l1 = FoxMomenta_l(_yybb,'Y',1);
	_YB_HY_l2 = FoxMomenta_l(_yybb,'Y',2);
	_YB_HY_l3 = FoxMomenta_l(_yybb,'Y',3);
	_YB_HY_l4 = FoxMomenta_l(_yybb,'Y',4);
	_YB_HY_l5 = FoxMomenta_l(_yybb,'Y',5);

	_YB_H1_l1 = FoxMomenta_l(_yybb,'1',1);
	_YB_H1_l2 = FoxMomenta_l(_yybb,'1',2);
	_YB_H1_l3 = FoxMomenta_l(_yybb,'1',3);
	_YB_H1_l4 = FoxMomenta_l(_yybb,'1',4);
	_YB_H1_l5 = FoxMomenta_l(_yybb,'1',5);

// HH 

	_HH_HT_l1 = FoxMomenta_l(_hh,'T',1);
	_HH_HT_l2 = FoxMomenta_l(_hh,'T',2);
	_HH_HT_l3 = FoxMomenta_l(_hh,'T',3);
	_HH_HT_l4 = FoxMomenta_l(_hh,'T',4);
	_HH_HT_l5 = FoxMomenta_l(_hh,'T',5);

	_HH_HP_l1 = FoxMomenta_l(_hh,'P',1);
	_HH_HP_l2 = FoxMomenta_l(_hh,'P',2);
	_HH_HP_l3 = FoxMomenta_l(_hh,'P',3);
	_HH_HP_l4 = FoxMomenta_l(_hh,'P',4);
	_HH_HP_l5 = FoxMomenta_l(_hh,'P',5);

	_HH_HS_l1 = FoxMomenta_l(_hh,'S',1);
	_HH_HS_l2 = FoxMomenta_l(_hh,'S',2);
	_HH_HS_l3 = FoxMomenta_l(_hh,'S',3);
	_HH_HS_l4 = FoxMomenta_l(_hh,'S',4);
	_HH_HS_l5 = FoxMomenta_l(_hh,'S',5);

	_HH_HZ_l1 = FoxMomenta_l(_hh,'Z',1);
	_HH_HZ_l2 = FoxMomenta_l(_hh,'Z',2);
	_HH_HZ_l3 = FoxMomenta_l(_hh,'Z',3);
	_HH_HZ_l4 = FoxMomenta_l(_hh,'Z',4);
	_HH_HZ_l5 = FoxMomenta_l(_hh,'Z',5);

	_HH_HY_l1 = FoxMomenta_l(_hh,'Y',1);
	_HH_HY_l2 = FoxMomenta_l(_hh,'Y',2);
	_HH_HY_l3 = FoxMomenta_l(_hh,'Y',3);
	_HH_HY_l4 = FoxMomenta_l(_hh,'Y',4);
	_HH_HY_l5 = FoxMomenta_l(_hh,'Y',5);

	_HH_H1_l1 = FoxMomenta_l(_hh,'1',1);
	_HH_H1_l2 = FoxMomenta_l(_hh,'1',2);
	_HH_H1_l3 = FoxMomenta_l(_hh,'1',3);
	_HH_H1_l4 = FoxMomenta_l(_hh,'1',4);
	_HH_H1_l5 = FoxMomenta_l(_hh,'1',5);


}
