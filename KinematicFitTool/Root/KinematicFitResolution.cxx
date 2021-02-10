///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2019-2020 CERN for the benefit of the ATLAS collaboration
*/

// KinematicFitResolution.cxx
// Source file for class KinematicFitResolution
// Author: BELFKIR Mohamed 
// Email : mohamed.belfkir@cern.ch
/////////////////////////////////////////////////////////////////// 


// KinematicFitResolution includes
#include <AsgTools/MessageCheck.h>
#include "KinematicFitTool/KinematicFitResolution.h"

KinematicFitResolution::KinematicFitResolution() : asg::AsgTool ("KinematicFitResolution"), m_photon_Res(0), m_angles_Res(0), m_2jet_constraint(0), m_3jet_constraint(0), m_EnergyResolutionName(""), m_EnergyResponseName("")
{
	declareProperty( "Photon Resolution", m_photon_Res);
	declareProperty( "Angles Resolution", m_angles_Res);
	declareProperty( "Fix Angles Fit", m_isFixAngles);
	declareProperty( "2 Jets Contraint", m_2jet_constraint);
	declareProperty( "3 Jets Contraint", m_3jet_constraint);
	declareProperty( "EnergyResolutionFile", m_EnergyResolutionName = ""); 
	declareProperty( "EnergyResponseFile", m_EnergyResponseName = "" ); 
}

KinematicFitResolution::~KinematicFitResolution() {

}

StatusCode KinematicFitResolution::initialize()
{

	if(m_photon_Res == 0) return StatusCode::FAILURE;
	if(m_angles_Res == 0) return StatusCode::FAILURE;
	if(m_2jet_constraint == 0 ) return StatusCode::FAILURE;
	if(m_3jet_constraint == 0 ) return StatusCode::FAILURE;
	if(m_EnergyResolutionName == "") return StatusCode::FAILURE;
	if(m_EnergyResponseName == "") return StatusCode::FAILURE;

	m_EnergyResolution = new TFile(m_EnergyResolutionName, "READ");
	m_EnergyResponse = new TFile(m_EnergyResponseName, "READ");

	ANA_MSG_INFO("KinematicFitResolution::initialize() : KinematicFitResolution is initialized!");
	return StatusCode::SUCCESS;
}
double KinematicFitResolution::PhotonRes()
{

	return m_photon_Res;
}

double KinematicFitResolution::AnglesRes()
{

	return m_angles_Res;
}

bool KinematicFitResolution::FixAngles()
{
	return m_isFixAngles;
}

double KinematicFitResolution::ConstraintValue(Int_t n)
{
	if(n == 0)
	{

		return m_2jet_constraint;
	}else{
		return m_3jet_constraint;
	}
}

double KinematicFitResolution::JetRes( xAOD::Jet* Jet)
{

	double rms = 0.1;
	double pt = log(Jet->pt()*1e-3);

        		
	float logPT_bin [] = {3.25,3.3,3.35,3.4,3.45,3.5,3.55,3.6,3.65,3.7,3.75,3.8,3.85,3.9,3.95,4.,4.05,4.1,4.15,4.2,4.25,4.3,4.35,4.4,4.45,4.5,4.55,4.6,4.65,4.7,4.75,4.8,4.85,4.9,4.95,5.,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.};
	int nbins = 46;
	int bin = -99;
	for(int i = 0; i<nbins; i++)
	{
		if(pt >= logPT_bin[i] && pt < logPT_bin[i+1])
		{
			bin = i;
			break;
		}
	}
	if(pt > 6.)
	{
		bin = nbins - 2;
	}
	if(pt < 3.25)
	{
		bin = 0;
	}
	
	TString histoname = Form("EReco_Inc_LogPt_bin_%2.2f_%2.2f",logPT_bin[bin],logPT_bin[bin+1]);
	TH1F *hist = (TH1F*) m_EnergyResolution->Get(histoname);
	rms = hist->GetRMS();

	delete hist;
	return ((rms != 0) ? rms : 0.1);
	
}

double KinematicFitResolution::JetMean( xAOD::Jet* Jet)
{

	double rms = 1.;
	double pt = log(Jet->pt()*1e-3);

        		
	float logPT_bin [] = {3.25,3.3,3.35,3.4,3.45,3.5,3.55,3.6,3.65,3.7,3.75,3.8,3.85,3.9,3.95,4.,4.05,4.1,4.15,4.2,4.25,4.3,4.35,4.4,4.45,4.5,4.55,4.6,4.65,4.7,4.75,4.8,4.85,4.9,4.95,5.,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.};
	int nbins = 46;
	int bin = -99;
	for(int i = 0; i<nbins; i++)
	{
		if(pt >= logPT_bin[i] && pt < logPT_bin[i+1])
		{
			bin = i;
			break;
		}
	}
	if(pt > 6.)
	{
		bin = nbins - 2;
	}
	if(pt < 3.25)
	{
		bin = 0;
	}
	
	TString histoname = Form("EReco_Inc_LogPt_bin_%2.2f_%2.2f",logPT_bin[bin],logPT_bin[bin+1]);
	TH1F *hist = (TH1F*) m_EnergyResolution->Get(histoname);
	rms = hist->GetMean();
	delete hist;
	return ((rms != 0) ? rms : 1.);
	
}

double KinematicFitResolution::JetResponse(TLorentzVector FitJet,  xAOD::Jet* Jet)
{

	double pt = log(JetMean(Jet)*Jet->pt()*1e-3);
        		
	float logPT_bin [] = {3.25,3.3,3.35,3.4,3.45,3.5,3.55,3.6,3.65,3.7,3.75,3.8,3.85,3.9,3.95,4.,4.05,4.1,4.15,4.2,4.25,4.3,4.35,4.4,4.45,4.5,4.55,4.6,4.65,4.7,4.75,4.8,4.85,4.9,4.95,5.,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.};
	int nbins = 46;
	int bin = -99;
	for(int i = 0; i<nbins; i++)
	{
		if(pt >= logPT_bin[i] && pt < logPT_bin[i+1])
		{
			bin = i;
			break;
		}
	}
	if(pt > 6.)
	{
		bin = nbins - 2;
	}
	if(pt < 3.25)
	{
		bin = 0;
	}
	TString histoname = Form("PTReco_hist_Inc_LogPt_bin_%2.2f_%2.2f",logPT_bin[bin],logPT_bin[bin+1]);
	TH1F *hist = (TH1F*) m_EnergyResponse->Get(histoname);
	double r = (Jet->pt()*JetMean(Jet))/FitJet.Pt();
	double prob = hist->Interpolate(r);
	prob *= 100;
	delete hist;
        return ((prob != 0) ? prob : 1e-20);

}



