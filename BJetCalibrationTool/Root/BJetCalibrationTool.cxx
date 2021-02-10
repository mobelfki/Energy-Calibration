///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2019-2020 CERN for the benefit of the ATLAS collaboration
*/

// BJetCalibrationTool.cxx
// Source file for class BJetCalibrationTool
// Author: BELFKIR Mohamed 
// Email : mohamed.belfkir@cern.ch
/////////////////////////////////////////////////////////////////// 


// BJetCalibrationTool includes
#include <AsgTools/MessageCheck.h>
#include "BJetCalibrationTool/BJetCalibrationTool.h"
#include "PathResolver/PathResolver.h"

// Constructors
////////////////

BJetCalibrationTool::BJetCalibrationTool(const std::string& name) : asg::AsgTool (name), m_JetAlgo(""), m_Jet_Min_Pt(), m_Jet_Min_Eta(), m_MuonContainer_Name(""), m_Muon_Quality(), m_Muon_Min_Pt(), m_Muon_Jet_dR() , m_doVR(true), m_doPtReco(true), m_PtReco(""), m_BJet_Tagger("")
{ 

	declareProperty( "JetCollection", m_JetAlgo = "AntiKt4EMPFlow" );
	declareProperty( "Jet_Min_Pt", m_Jet_Min_Pt = 20. );
	declareProperty( "Jet_Min_Eta", m_Jet_Min_Eta = 2.5 );
        declareProperty( "MuonContainer", m_MuonContainer_Name = "Muons" );
	declareProperty( "MuonQuality", m_Muon_Quality = xAOD::Muon::Medium);
	declareProperty( "Muon_Min_Pt", m_Muon_Min_Pt = 3.);
	declareProperty( "Jet_Muon_dR", m_Muon_Jet_dR = 0.4);
	declareProperty( "BJet_Tagger", m_BJet_Tagger = "DL1r");
	declareProperty( "doVR", m_doVR = true);
        declareProperty( "doPtReco", m_doPtReco = true);
	declareProperty( "PtRecoFile", m_PtReco = "PtReco_Correction_MV2c10_DL1r_05032020.root");

}


// Destructor
///////////////
BJetCalibrationTool::~BJetCalibrationTool() {


}


StatusCode BJetCalibrationTool::initialize() {
	
  ANA_MSG_INFO ("Initializing BJetCalibrationTool ... " << name() << " ... "); 
  ANA_CHECK( this->initializeTool( name() ) );
  return StatusCode::SUCCESS;
}

StatusCode BJetCalibrationTool::finalize() {
 

  return StatusCode::SUCCESS;
}


StatusCode BJetCalibrationTool::initializeTool(const std::string& name) {

  if(m_JetAlgo == "")
  {
	ANA_MSG_FATAL("BJetCalibrationTool::initialize() : Please set JetCollection"); 
	return StatusCode::FAILURE;
  }

  if(m_PtReco == "")
  {
	ANA_MSG_FATAL("BJetCalibrationTool::initialize() : Please set PtReco file name");
	return StatusCode::FAILURE;
  }
  m_muonSelection.setName("MuonSelectionToolForBJetCalibrationTool");
  ANA_CHECK( m_muonSelection.setProperty( "MaxEta", 2.5 ));
  ANA_CHECK( m_muonSelection.setProperty( "MuQuality",(int) m_Muon_Quality));
  ANA_CHECK( m_muonSelection.initialize());

  if(m_doPtReco)
  {
  	TString filename = PathResolverFindCalibFile("BJetCalibrationTool/"+m_JetAlgo+"_"+m_PtReco);
	m_PtRecoFile = new TFile(filename,"READ");
	if(!m_PtRecoFile)
        {
		ANA_MSG_FATAL("BJetCalibrationTool::initialize() : PtReco file is not found" << filename );
		return StatusCode::FAILURE;
        }
	m_Semi_Histo = (TH1F*) m_PtRecoFile->Get("Correction_SemiLeptonic_"+m_BJet_Tagger+"_ttbar_mean"); 
	m_Had_Histo = (TH1F*) m_PtRecoFile->Get("Correction_Hadronic_"+m_BJet_Tagger+"_ttbar_mean"); 	
	
	if(m_Semi_Histo == NULL or m_Had_Histo == NULL)
	{
		ANA_MSG_FATAL("BJetCalibrationTool::initialize() : Please check the histograms names");
		return StatusCode::FAILURE;
	}

  }
  ANA_MSG_INFO("BJetCalibrationTool::initialize() : " << name << " is initialized!");
  return StatusCode::SUCCESS;
}

StatusCode BJetCalibrationTool::initializeMuonContainer(std::vector< const xAOD::Muon* >& muon)
{
	muon.clear();
	const xAOD::MuonContainer* MuonContainer = 0;	
	if( evtStore()->retrieve( MuonContainer, m_MuonContainer_Name).isFailure() )
	{
		ANA_MSG_FATAL("BJetCalibrationTool::initializeMuonContainer() : Please check MuonContainer name");
		return StatusCode::FAILURE;
	}
	for( auto mu : *MuonContainer )
	{
		xAOD::Muon::Quality muonQualtiy = m_muonSelection.getQuality(*mu);
		if( mu->pt()*1e-3 > m_Muon_Min_Pt && muonQualtiy <= m_Muon_Quality)
		{
			muon.push_back(mu);
		}
	}
	return StatusCode::SUCCESS;
}

StatusCode BJetCalibrationTool::applyBJetCalibration(xAOD::Jet& jet) { 

  std::vector< const xAOD::Muon* > muons;
  ANA_CHECK(initializeMuonContainer(muons));
  std::vector< const xAOD::Muon* > muons_in_jet = getMuonInJet(jet, muons);

  addMuon(jet, muons_in_jet);
 
  if(m_doPtReco)
  {
	applyPtReco(jet);
  }

  muons.clear();
  muons_in_jet.clear();
  return StatusCode::SUCCESS; 
}

void BJetCalibrationTool::applyPtReco(xAOD::Jet& jet)
{
	static SG::AuxElement::Decorator<Float_t> PtRecoF("PtReco_SF");	

	Int_t nmu = jet.auxdata<Int_t>("n_muons");
	Float_t PtRecoFactor = 1.;
	TLorentzVector j = jet.p4();

	if(nmu == 0)
	{
		PtRecoFactor = m_Had_Histo->Interpolate(log(j.Pt() * 0.001));

	}else if(nmu > 0){
	
		PtRecoFactor = m_Semi_Histo->Interpolate(log(j.Pt() * 0.001));
	}	
		
	j *= PtRecoFactor;
	
	xAOD::JetFourMom_t new_jet(j.Pt(),j.Eta(),j.Phi(),j.M());
	jet.setJetP4(new_jet);
	PtRecoF(jet) = PtRecoFactor;

}

std::vector< const xAOD::Muon* > BJetCalibrationTool::getMuonInJet(xAOD::Jet& jet, std::vector< const xAOD::Muon* > muons)
{

	std::vector< const xAOD::Muon* > m_muon_in_jet;

	for( unsigned int i = 0; i<muons.size(); i++)
	{
		TLorentzVector mu = muons[i]->p4();	
		TLorentzVector j = jet.p4();	
		Double_t dR = j.DeltaR(mu);
			
		if(m_doVR)
		{
		   m_Muon_Jet_dR = std::min(0.4, 0.04 + (10/ (((mu).Pt() * 0.001))));	
		} 		
		
		if(dR > m_Muon_Jet_dR ) continue;

		m_muon_in_jet.push_back(muons[i]);
		
	}
	return m_muon_in_jet;
}
void BJetCalibrationTool::addMuon(xAOD::Jet& jet, std::vector< const xAOD::Muon* > muons)
{
	static SG::AuxElement::Decorator<Int_t> NMu("n_muons");	
	
	NMu(jet) = muons.size();	
		
	if(jet.auxdata<Int_t>("n_muons") == 0)
	{
		return;
	}

	if(jet.auxdata<Int_t>("n_muons") == 1)
	{
		
		float eLoss = 0.0;
		muons[0]->parameter(eLoss,xAOD::Muon::EnergyLoss);
		TLorentzVector mu = muons[0]->p4();
		double theta=mu.Theta();
  		double phi=mu.Phi();
  		double eLossX=eLoss*sin(theta)*cos(phi);
  		double eLossY=eLoss*sin(theta)*sin(phi);
  		double eLossZ=eLoss*cos(theta);
		TLorentzVector Loss = TLorentzVector(eLossX,eLossY,eLossZ,eLoss);
		
		TLorentzVector j = jet.p4();

		j=j-Loss+mu;		
		
		xAOD::JetFourMom_t new_jet(j.Pt(),j.Eta(),j.Phi(),j.M());

		jet.setJetP4(new_jet);
		return;		
	}
	
	if(jet.auxdata<Int_t>("n_muons") > 1)
	{
		const xAOD::Muon* muon = getMuon(muons);
	
		TLorentzVector mu = muon->p4();
		TLorentzVector j = jet.p4();
		float eLoss = 0.0;
		muon->parameter(eLoss,xAOD::Muon::EnergyLoss);
		double theta=mu.Theta();
  		double phi=mu.Phi();
  		double eLossX=eLoss*sin(theta)*cos(phi);
  		double eLossY=eLoss*sin(theta)*sin(phi);
  		double eLossZ=eLoss*cos(theta);
		TLorentzVector Loss = TLorentzVector(eLossX,eLossY,eLossZ,eLoss);

		j= j-Loss+mu;		
		
		xAOD::JetFourMom_t new_jet(j.Pt(),j.Eta(),j.Phi(),j.M());

		jet.setJetP4(new_jet);

		return;		
	}
}
const xAOD::Muon* BJetCalibrationTool::getMuon(std::vector< const xAOD::Muon* > muons)
{
	Int_t muon_size = muons.size();
	Double_t max_muon_pt = muons[0]->pt();
	Int_t index = 0;	
	for(Int_t i = 0; i<muon_size; i++)
	{
		if(max_muon_pt < muons[i]->pt())
		{
			max_muon_pt = muons[i]->pt();
			index = i;
		}
	}
		return muons[index];
 
}
