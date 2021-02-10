///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2019-2020 CERN for the benefit of the ATLAS collaboration
*/

// KinematicFitRun.cxx
// Source file for class KinematicFitRun
// Magic happen here!
// Author: BELFKIR Mohamed 
// Email : mohamed.belfkir@cern.ch
/////////////////////////////////////////////////////////////////// 


// KinematicFitRun includes
#include <AsgTools/MessageCheck.h>
#include "KinematicFitTool/KinematicFitRun.h"

class KinematicFitRun* gThis;

KinematicFitRun::KinematicFitRun() : minuit(NULL)
{ 

}

KinematicFitRun::~KinematicFitRun() {

}


StatusCode KinematicFitRun::initialize() {

  return StatusCode::SUCCESS;
}

StatusCode KinematicFitRun::finalize() {
 

  return StatusCode::SUCCESS;
}

StatusCode KinematicFitRun::RunKF(KinematicFitEvent*& Event, KinematicFitResolution*& Resolution)
{
	
	if(Event->GetAddJets().size() >= 2)
	{
		return StatusCode::SUCCESS;
	}
	NewEvent();
  	SetEventAndResolution(Event,Resolution);
	Run();
   	return StatusCode::SUCCESS; 
}

StatusCode KinematicFitRun::SetEventAndResolution(KinematicFitEvent*& Event, KinematicFitResolution*& Resolution)
{
	_Event = Event;
	_Resolution = Resolution; 
	return StatusCode::SUCCESS;
}

void KinematicFitRun::Run()
{
	gThis=this;
	isfitted = false;
	Int_t npar = GetNParameter();
	int conv = -99;
	if(minuit) delete minuit;
	minuit = new TMinuit(npar);
	conv = minuit->SetPrintLevel(-1);
	minuit->SetFCN(&KinematicFitRun::FCN);
	minuit->SetErrorDef(0.5);
	Double_t arglist[npar];
	arglist[0] = 2;
	minuit->mnexcm("SET STR", arglist, 1, g);
	SetParameterNameandRange(minuit);
	arglist[0] = 50000.;
	arglist[1] = 1.;
	minuit->mnexcm("MIGRAD",arglist, 2, conv);
	isfitted = true;
	ShareFit(minuit);
	delete minuit;
}

void KinematicFitRun::SetParameterNameandRange(TMinuit*& m)
{

	m->mnparm(0, "EFit_photon1", _Event->GetPhotonOne()->e(), .01*_Resolution->PhotonRes()*_Event->GetPhotonOne()->e(), .1*_Resolution->PhotonRes()*_Event->GetPhotonOne()->e(), 2.*_Event->GetPhotonOne()->e(), g);
	m->mnparm(1, "EFit_photon2", _Event->GetPhotonTwo()->e(), .01*_Resolution->PhotonRes()*_Event->GetPhotonTwo()->e(), .1*_Resolution->PhotonRes()*_Event->GetPhotonTwo()->e(), 2.*_Event->GetPhotonTwo()->e(), g);

	m->mnparm(2, "EFit_jet1", _Event->GetBJetOne()->e(), .01*_Event->GetBJetOne()->e()*GetJetRes(_Event->GetBJetOne()), .1*_Event->GetBJetOne()->e()*GetJetRes(_Event->GetBJetOne()), 2.*_Event->GetBJetOne()->e(), g);
	m->mnparm(3, "EFit_jet2", _Event->GetBJetTwo()->e(), .01*_Event->GetBJetTwo()->e()*GetJetRes(_Event->GetBJetTwo()), .1*_Event->GetBJetTwo()->e()*GetJetRes(_Event->GetBJetTwo()), 2.*_Event->GetBJetTwo()->e(), g);


	m->mnparm(4, "EtaFit_jet1", _Event->GetBJetOne()->eta(), _Resolution->AnglesRes(), -7., 7., g);
	m->mnparm(5, "EtaFit_jet2", _Event->GetBJetTwo()->eta(), _Resolution->AnglesRes(), -7., 7., g);
	m->mnparm(6, "EtaFit_photon1", _Event->GetPhotonOne()->eta(), _Resolution->AnglesRes(), -7., 7., g);
	m->mnparm(7, "EtaFit_photon2", _Event->GetPhotonTwo()->eta(), _Resolution->AnglesRes(), -7., 7., g);


	m->mnparm(8, "PhiFit_jet1", _Event->GetBJetOne()->phi(), _Resolution->AnglesRes(), -1*TMath::Pi(), TMath::Pi(), g);
	m->mnparm(9, "PhiFit_jet2", _Event->GetBJetTwo()->phi(), _Resolution->AnglesRes(), -1*TMath::Pi(), TMath::Pi(), g);
	m->mnparm(10, "PhiFit_photon1",_Event->GetPhotonOne()->phi(), _Resolution->AnglesRes(), -1*TMath::Pi(), TMath::Pi(), g);
	m->mnparm(11, "PhiFit_photon2",_Event->GetPhotonTwo()->phi(), _Resolution->AnglesRes(), -1*TMath::Pi(), TMath::Pi(), g);


	if(_Resolution->FixAngles() == true)
	{
		m->FixParameter(4);
		m->FixParameter(5);
		m->FixParameter(6);
		m->FixParameter(7);

		m->FixParameter(8);
		m->FixParameter(9);
		m->FixParameter(10);
		m->FixParameter(11);
	}

	if(_Event->GetAddJets().size() == FitType::ThreeJet)
	{
		m->mnparm(12, "EFit_jet3", _Event->GetAddJets().at(0)->e(), .01*GetJetRes(_Event->GetAddJets().at(0)), .1*GetJetRes(_Event->GetAddJets().at(0)), 2.*_Event->GetAddJets().at(0)->e(), g);
		m->mnparm(13, "EtaFit_jet3", _Event->GetAddJets().at(0)->eta(), _Resolution->AnglesRes(), -7., 7., g);
		m->mnparm(14, "PhiFit_jet3", _Event->GetAddJets().at(0)->phi(), _Resolution->AnglesRes(), -1*TMath::Pi(), TMath::Pi(), g);

		if(_Resolution->FixAngles() == true)
		{
			m->FixParameter(13);
			m->FixParameter(14);
		}
	}

}


double KinematicFitRun::GetLH()
{

	Double_t LH = 0., PxHH = 0., PyHH = 0., PtHH = 0.;

	LH += TMath::Power((_FitBJet1.E() - _Event->GetBJetOne()->e()*GetJetMean(_Event->GetBJetOne())),2)/TMath::Power(_Event->GetBJetOne()->e()*(GetJetRes(_Event->GetBJetOne())),2);
	LH += TMath::Power((_FitBJet2.E() - _Event->GetBJetTwo()->e()*GetJetMean(_Event->GetBJetTwo())),2)/TMath::Power(_Event->GetBJetTwo()->e()*(GetJetRes(_Event->GetBJetTwo())),2);

	LH += TMath::Power((_FitPhoton1.E() - _Event->GetPhotonOne()->e()),2)/TMath::Power(_Resolution->PhotonRes()*_Event->GetPhotonOne()->e(),2);
	LH += TMath::Power((_FitPhoton2.E() - _Event->GetPhotonTwo()->e()),2)/TMath::Power(_Resolution->PhotonRes()*_Event->GetPhotonTwo()->e(),2);

	LH += TMath::Power((_FitBJet1.Eta()-_Event->GetBJetOne()->eta()),2)/TMath::Power((_Resolution->AnglesRes()*fabs(_Event->GetBJetOne()->eta())),2);
	LH += TMath::Power((_FitBJet2.Eta()-_Event->GetBJetTwo()->eta()),2)/TMath::Power((_Resolution->AnglesRes()*fabs(_Event->GetBJetTwo()->eta())),2);
	
	LH += TMath::Power((_FitPhoton1.Eta()-_Event->GetPhotonOne()->eta()),2)/TMath::Power((_Resolution->AnglesRes()*fabs(_Event->GetPhotonOne()->eta())),2);
	LH += TMath::Power((_FitPhoton2.Eta()-_Event->GetPhotonTwo()->eta()),2)/TMath::Power((_Resolution->AnglesRes()*fabs(_Event->GetPhotonTwo()->eta())),2);
	
	LH += TMath::Power((_FitBJet1.Phi()-_Event->GetBJetOne()->phi()),2)/TMath::Power(_Resolution->AnglesRes(),2);
	LH += TMath::Power((_FitBJet2.Phi()-_Event->GetBJetTwo()->phi()),2)/TMath::Power(_Resolution->AnglesRes(),2);

	LH += TMath::Power((_FitPhoton1.Phi()-_Event->GetPhotonOne()->phi()),2)/TMath::Power(_Resolution->AnglesRes(),2);
	LH += TMath::Power((_FitPhoton2.Phi()-_Event->GetPhotonTwo()->phi()),2)/TMath::Power(_Resolution->AnglesRes(),2);

	LH += -2*log(GetResponse(_FitBJet1, _Event->GetBJetOne()));
	LH += -2*log(GetResponse(_FitBJet2, _Event->GetBJetTwo()));

	PxHH = (_FitBJet1 + _FitBJet2 + _FitPhoton1 + _FitPhoton2).Px();
	PyHH = (_FitBJet1 + _FitBJet2 + _FitPhoton1 + _FitPhoton2).Py();
	
	if(_Event->GetAddJets().size() == FitType::ThreeJet)
	{
		
		LH += TMath::Power((_FitJet3.E() - _Event->GetAddJets().at(0)->e()*GetJetMean(_Event->GetAddJets().at(0))),2)/TMath::Power(_Event->GetAddJets().at(0)->e()*(GetJetRes(_Event->GetAddJets().at(0))),2);
		LH += TMath::Power((_FitJet3.Eta() - _Event->GetAddJets().at(0)->eta()),2)/TMath::Power((_Resolution->AnglesRes()*fabs(_Event->GetAddJets().at(0)->eta())),2);
		LH += TMath::Power((_FitJet3.Phi() - _Event->GetAddJets().at(0)->phi()),2)/TMath::Power(_Resolution->AnglesRes(),2);

		PxHH = (_FitBJet1 + _FitBJet2 + _FitPhoton1 + _FitPhoton2 + _FitJet3).Px();
		PyHH = (_FitBJet1 + _FitBJet2 + _FitPhoton1 + _FitPhoton2 + _FitJet3).Py();
	}

	
	LH += TMath::Power((PxHH - 0.0),2)/TMath::Power(GetConstraintRes(_Event->GetAddJets())*GeV, 2);
	LH += TMath::Power((PyHH - 0.0),2)/TMath::Power(GetConstraintRes(_Event->GetAddJets())*GeV, 2);

	return LH;

}

void KinematicFitRun::SetParameters(double* par)
{

	
	_FitPhoton1.SetPtEtaPhiE(par[0]*sin(2*atan(exp(-1*par[6]))), par[6], par[10], par[0]);
	_FitPhoton2.SetPtEtaPhiE(par[1]*sin(2*atan(exp(-1*par[7]))), par[7], par[11], par[1]);
	_FitBJet1.SetPtEtaPhiE(par[2]*sin(2*atan(exp(-1*par[4]))), par[4], par[8], par[2]);
	_FitBJet2.SetPtEtaPhiE(par[3]*sin(2*atan(exp(-1*par[5]))), par[5], par[9], par[3]);
	_FitJet3.SetPtEtaPhiE(par[12]*sin(2*atan(exp(-1*par[13]))), par[13], par[14], par[12]);
	
}

void KinematicFitRun::FCN(int& /*npar*/, double* /*grad*/, double& fval, double* par, int /*flag*/)
{

	
	gThis->SetParameters(par);
	
	fval = gThis->GetLH(); 
	
}

double KinematicFitRun::GetConstraintRes(std::vector<  xAOD::Jet*> AddJets)
{

	int n_add_jet = AddJets.size();
	
	double Res = _Resolution->ConstraintValue(n_add_jet);	
	
	return Res;
}

double KinematicFitRun::GetJetRes( xAOD::Jet* Jet)
{

	return _Resolution->JetRes(Jet);
}

double KinematicFitRun::GetResponse(TLorentzVector FitJet,  xAOD::Jet* Jet)
{
	return _Resolution->JetResponse(FitJet, Jet);
	
}

double KinematicFitRun::GetJetMean( xAOD::Jet* Jet)
{
	return _Resolution->JetMean(Jet);

}

void KinematicFitRun::ShareFit(TMinuit*& m)
{

	double eta, phi, E, err, pt;

	m->GetParameter(0,E,err);
	m->GetParameter(6,eta,err);
	m->GetParameter(10,phi,err);
	pt = E / cosh(eta);

	_FitPhoton1.SetPtEtaPhiE(pt,eta,phi,E);

	m->GetParameter(1,E,err);
	m->GetParameter(7,eta,err);
	m->GetParameter(11,phi,err);
	pt = E / cosh(eta);

	_FitPhoton2.SetPtEtaPhiE(pt,eta,phi,E);

	m->GetParameter(2,E,err);
	m->GetParameter(4,eta,err);
	m->GetParameter(8,phi,err);
	pt = E / cosh(eta);

	_FitBJet1.SetPtEtaPhiE(pt,eta,phi,E);

	m->GetParameter(3,E,err);
	m->GetParameter(5,eta,err);
	m->GetParameter(9,phi,err);
	pt = E / cosh(eta);

	_FitBJet2.SetPtEtaPhiE(pt,eta,phi,E);

	if(_Event->GetAddJets().size() == FitType::ThreeJet)
	{
	
		m->GetParameter(12,E,err);
		m->GetParameter(13,eta,err);
		m->GetParameter(14,phi,err);
		pt = E / cosh(eta);

		_FitJet3.SetPtEtaPhiE(pt,eta,phi,E);
	}else{
		_FitJet3.SetPtEtaPhiE(0,0,0,0);
	}
	
	_Event->SetFitPhotonOne(_FitPhoton1);
	_Event->SetFitPhotonTwo(_FitPhoton2);

	_Event->SetFitBJetOne(_FitBJet1);
	_Event->SetFitBJetTwo(_FitBJet2);

	_Event->SetFitAddJet(_FitJet3);

}
void KinematicFitRun::NewEvent()
{
	_FitPhoton1.SetPtEtaPhiE(0., 0., 0., 0.);
	_FitPhoton2.SetPtEtaPhiE(0., 0., 0., 0.);
	_FitBJet1.SetPtEtaPhiE(0., 0., 0., 0.);
	_FitBJet2.SetPtEtaPhiE(0., 0., 0., 0.);
	_FitJet3.SetPtEtaPhiE(0., 0., 0., 0.);
}

Int_t KinematicFitRun::GetNParameter()
{

	int n = 0;
	if(_Event->GetAddJets().size() == FitType::ThreeJet)
	{	
		n = 15;
		
		if(_Resolution->FixAngles() == true)
		{
			n = 5;
		}
	}

	if(_Event->GetAddJets().size() == FitType::TwoJet)
	{	
		n = 12;
		
		if(_Resolution->FixAngles() == true)
		{
			n = 4;
		}
	}
	
	return n;
}

Bool_t KinematicFitRun::isFitted()
{
	return isfitted;
}


