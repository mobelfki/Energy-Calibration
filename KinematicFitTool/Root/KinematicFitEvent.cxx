///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2019-2020 CERN for the benefit of the ATLAS collaboration
*/

// KinematicFitEvent.cxx
// Source file for class KinematicFitEvent
// Author: BELFKIR Mohamed 
// Email : mohamed.belfkir@cern.ch
/////////////////////////////////////////////////////////////////// 


// KinematicFitEvent includes
#include <AsgTools/MessageCheck.h>
#include "KinematicFitTool/KinematicFitEvent.h"

// Constructors
////////////////

KinematicFitEvent::KinematicFitEvent() : asg::AsgTool ("KinematicFitEvent"), m_Jet_Min_Pt(0), m_BtaggingWP("")
{ 

	declareProperty( "Jet Min Pt", m_Jet_Min_Pt);
	declareProperty( "BTagging WP", m_BtaggingWP);
}

KinematicFitEvent::~KinematicFitEvent() {

}


StatusCode KinematicFitEvent::initialize() {

	if(m_Jet_Min_Pt == 0) return StatusCode::FAILURE;
	if(m_BtaggingWP == "") return StatusCode::FAILURE;
	ANA_MSG_INFO("KinematicFitEvent::initialize() : KinematicFitEvent is initialized!");
  return StatusCode::SUCCESS;
}

StatusCode KinematicFitEvent::finalize() {
 

  return StatusCode::SUCCESS;
}

StatusCode KinematicFitEvent::applySelection(const xAOD::PhotonContainer& photons, xAOD::JetContainer& jets) {

	int n_photons = photons.size();
	int n_jets    = jets.size();
	if(n_photons < 2) 
	{
		return StatusCode::FAILURE;
	}
	if(n_jets < 2)
	{
		return StatusCode::FAILURE;
	}
	
	m_photon1 = photons.get(0);
	m_photon2 = photons.get(1);

	std::vector<  xAOD::Jet*> b_jets;
	std::vector<  xAOD::Jet*> nb_jets;
	
	for( auto jet : jets)
	{
		if( jet->pt()*1e-3 < m_Jet_Min_Pt) continue;
		
		if(jet->auxdata<int>("DL1rbin") > 0 )
		{
			b_jets.push_back(jet);
			
		}else{
			nb_jets.push_back(jet);
		}		
	}
	if(b_jets.size() < 2)
	{
		return StatusCode::FAILURE;
	}

	std::sort(b_jets.begin(), b_jets.end(), []( auto & i,  auto & j) { return i->pt() > j->pt(); });
	std::sort(nb_jets.begin(), nb_jets.end(), []( auto & i,  auto & j) { return i->pt() > j->pt(); });
	
	m_bjet1 = b_jets.at(0);
	m_bjet2 = b_jets.at(1);

	for(unsigned i = 2; i<b_jets.size(); i++)
	{
		m_addJ.push_back(b_jets.at(i));

	}

	for(unsigned i = 0; i<nb_jets.size(); i++)
	{
		m_addJ.push_back(nb_jets.at(i));
	}
	 
  return StatusCode::SUCCESS; 
}

const xAOD::Photon* KinematicFitEvent::GetPhotonOne()
{
	return m_photon1;
}
const xAOD::Photon* KinematicFitEvent::GetPhotonTwo()
{

	return m_photon2;
}

 xAOD::Jet* KinematicFitEvent::GetBJetTwo()
{

	return m_bjet2;
}
 xAOD::Jet* KinematicFitEvent::GetBJetOne()
{

	return m_bjet1;
}

std::vector<  xAOD::Jet*> KinematicFitEvent::GetAddJets()
{
	return m_addJ;
}


TLorentzVector KinematicFitEvent::GetFitPhotonOne()
{
	return m_photon1_Fit;
}

TLorentzVector KinematicFitEvent::GetFitPhotonTwo()
{

	return m_photon2_Fit;
}

TLorentzVector KinematicFitEvent::GetFitBJetTwo()
{

	return m_bjet2_Fit;
}
TLorentzVector KinematicFitEvent::GetFitBJetOne()
{

	return m_bjet1_Fit;
}

TLorentzVector KinematicFitEvent::GetFitAddJet()
{
	return m_addJ_Fit;
}

void KinematicFitEvent::SetFitPhotonOne(TLorentzVector tlv)
{
	m_photon1_Fit = tlv;
}
void KinematicFitEvent::SetFitPhotonTwo(TLorentzVector tlv)
{

	m_photon2_Fit = tlv;
}
void KinematicFitEvent::SetFitBJetTwo(TLorentzVector tlv)
{

	m_bjet2_Fit = tlv;
}
void KinematicFitEvent::SetFitBJetOne(TLorentzVector tlv)
{

	m_bjet1_Fit = tlv;
}

void KinematicFitEvent::SetFitAddJet(TLorentzVector tlv)
{
	m_addJ_Fit = tlv;
}

void KinematicFitEvent::Clear()
{
	
	m_addJ.clear();

}
