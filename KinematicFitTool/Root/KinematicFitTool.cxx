///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2019-2020 CERN for the benefit of the ATLAS collaboration
*/

// KinematicFitTool.cxx
// Source file for class KinematicFitTool
// Author: BELFKIR Mohamed 
// Email : mohamed.belfkir@cern.ch
/////////////////////////////////////////////////////////////////// 


// KinematicFitTool includes
#include <AsgTools/MessageCheck.h>

#include "KinematicFitTool/KinematicFitTool.h"
#include "PathResolver/PathResolver.h"

KinematicFitTool::KinematicFitTool(const std::string& name) : asg::AsgTool (name), m_JetAlgo(""), m_Jet_Min_Pt(20.), m_photon_Res(0.01), m_angles_Res(0.01), m_isFixAngles(true), m_2jet_constraint(14.), m_3jet_constraint(16.), m_file_name(), m_BtaggingWP()
{ 
	declareProperty( "Jet Collection", m_JetAlgo = "AntiKt4EMPFlow" );
	declareProperty( "Jet Min Pt", m_Jet_Min_Pt = 20. );
	declareProperty( "BTagging WP", m_BtaggingWP = "MV2c10_FixedCutBEff_70" );
	declareProperty( "Photon Resolution", m_photon_Res = .01 );
	declareProperty( "Angles Resolution", m_angles_Res = .01 );
	declareProperty( "Fix Angles Fit", m_isFixAngles = true );
	declareProperty( "2 Jets Contraint", m_2jet_constraint = 14. );
	declareProperty( "3 Jets Contraint", m_3jet_constraint = 16. );
	declareProperty( "Files Name", m_file_name = "KinematicFitTool.root" );
}

KinematicFitTool::~KinematicFitTool() {


}


StatusCode KinematicFitTool::initialize() {
	
  ANA_MSG_INFO ("Initializing KinematicFit ... " << name() << " ... "); 
  ANA_CHECK( this->initializeTool( name() ) );
  return StatusCode::SUCCESS;
}

StatusCode KinematicFitTool::finalize() {
 

  return StatusCode::SUCCESS;
}


StatusCode KinematicFitTool::initializeTool(const std::string& name) {

  if(m_JetAlgo == "")
  {
	ANA_MSG_FATAL("KinematicFitTool::initialize() : Please set JetCollection"); 
	return StatusCode::FAILURE;
  }

  std::string m_EnergyResolution = PathResolverFindCalibFile("KinematicFitTool/EnergyResolution_"+m_JetAlgo+"_"+m_file_name);
  std::string m_EnergyResponse   = PathResolverFindCalibFile("KinematicFitTool/EnergyResponse_"+m_JetAlgo+"_"+m_file_name);


  fResolution = new KinematicFitResolution();
 	ANA_CHECK(fResolution->setProperty("Photon Resolution", m_photon_Res));
	ANA_CHECK(fResolution->setProperty("Angles Resolution", m_angles_Res));
	ANA_CHECK(fResolution->setProperty("Fix Angles Fit", m_isFixAngles));
	ANA_CHECK(fResolution->setProperty("2 Jets Contraint", m_2jet_constraint));
	ANA_CHECK(fResolution->setProperty("3 Jets Contraint", m_3jet_constraint));
	ANA_CHECK(fResolution->setProperty("EnergyResolutionFile", m_EnergyResolution)); 
	ANA_CHECK(fResolution->setProperty("EnergyResponseFile", m_EnergyResponse));
	ANA_CHECK(fResolution->initialize());

 fEvent = new KinematicFitEvent();
	ANA_CHECK(fEvent->setProperty("Jet Min Pt", m_Jet_Min_Pt)); 
	ANA_CHECK(fEvent->setProperty("BTagging WP", m_BtaggingWP));
	ANA_CHECK(fEvent->initialize());

  ANA_MSG_INFO("KinematicFitTool::initialize() : " << name << " is initialized!");

  ANA_MSG_INFO("KinematicFitTool::initialize() : KinematicFitRun is Ready!");
  return StatusCode::SUCCESS;
}


StatusCode KinematicFitTool::applyKF(const xAOD::PhotonContainer& photons, xAOD::JetContainer& jets) 
{

	static SG::AuxElement::Decorator<Float_t> PT("KF_PT");	
	static SG::AuxElement::Decorator<Float_t> ETA("KF_ETA");	
	static SG::AuxElement::Decorator<Float_t> PHI("KF_PHI");	
	static SG::AuxElement::Decorator<Float_t> M("KF_M");	
	static SG::AuxElement::Decorator<Char_t> isB("KF_isB");	

	for( auto jet : jets)
	{
		PT(*jet)    = jet->pt();
		ETA(*jet)   = jet->eta();
		PHI(*jet)   = jet->phi();
		M(*jet)     = jet->m();
		isB(*jet)   = false;
	}
			
	for( auto ph : photons)
	{
		PT(*ph)    = ph->pt();
		ETA(*ph)   = ph->eta();
		PHI(*ph)   = ph->phi();
	}
	if(!fEvent->applySelection(photons, jets)) return StatusCode::SUCCESS;
	

	fProcessor = new KinematicFitRun();

	fProcessor->RunKF(fEvent, fResolution);

	DecorateEvent(fEvent, fProcessor);

	fEvent->Clear();
	delete fProcessor;
  return StatusCode::SUCCESS; 
}

void KinematicFitTool::DecorateEvent(KinematicFitEvent*& Event, KinematicFitRun*& Run)
{
	static SG::AuxElement::Decorator<Float_t> PT("KF_PT");	
	static SG::AuxElement::Decorator<Float_t> ETA("KF_ETA");	
	static SG::AuxElement::Decorator<Float_t> PHI("KF_PHI");	
	static SG::AuxElement::Decorator<Float_t> M("KF_M");	
	static SG::AuxElement::Decorator<Char_t> isB("KF_isB");	

	bool isFitted = Run->isFitted();
	if(isFitted)
	{

		TLorentzVector jet1;
		TLorentzVector jet2;

		jet1.SetPtEtaPhiE((Event->GetFitBJetOne()).Pt(),(Event->GetFitBJetOne()).Eta(),(Event->GetFitBJetOne()).Phi(),(Event->GetFitBJetOne()).E());
		jet2.SetPtEtaPhiE((Event->GetFitBJetTwo()).Pt(),(Event->GetFitBJetTwo()).Eta(),(Event->GetFitBJetTwo()).Phi(),(Event->GetFitBJetTwo()).E());

		xAOD::JetFourMom_t jet14vec(jet1.Pt(),jet1.Eta(),jet1.Phi(),jet1.M());
		xAOD::JetFourMom_t jet24vec(jet2.Pt(),jet2.Eta(),jet2.Phi(),jet2.M());

		(Event->GetBJetOne())->setJetP4(jet14vec);
		(Event->GetBJetTwo())->setJetP4(jet24vec);

		PT(*(Event->GetBJetOne()))    = Event->GetFitBJetOne().Pt();
		ETA(*(Event->GetBJetOne()))   = Event->GetFitBJetOne().Eta();
		PHI(*(Event->GetBJetOne()))   = Event->GetFitBJetOne().Phi();
		M(*(Event->GetBJetOne()))     = Event->GetFitBJetOne().M();
		PT(*(Event->GetBJetTwo()))    = Event->GetFitBJetTwo().Pt();
		ETA(*(Event->GetBJetTwo()))   = Event->GetFitBJetTwo().Eta();
		PHI(*(Event->GetBJetTwo()))   = Event->GetFitBJetTwo().Phi();
		M(*(Event->GetBJetTwo()))     = Event->GetFitBJetTwo().M();
		isB(*(Event->GetBJetTwo()))   = true;
		isB(*(Event->GetBJetOne()))   = true;

		PT(*(Event->GetPhotonOne()))  = Event->GetFitPhotonOne().Pt();
		ETA(*(Event->GetPhotonOne())) = Event->GetFitPhotonOne().Eta();
		PHI(*(Event->GetPhotonOne())) = Event->GetFitPhotonOne().Phi();
		PT(*(Event->GetPhotonTwo()))  = Event->GetFitPhotonTwo().Pt();
		ETA(*(Event->GetPhotonTwo())) = Event->GetFitPhotonTwo().Eta();
		PHI(*(Event->GetPhotonTwo())) = Event->GetFitPhotonTwo().Phi();
		

		if(Event->GetAddJets().size() == FitType::ThreeJet)
		{

			TLorentzVector jet;
			jet.SetPtEtaPhiE((Event->GetFitAddJet()).Pt(),(Event->GetFitAddJet()).Eta(),(Event->GetFitAddJet()).Phi(),(Event->GetFitAddJet()).E());
		
			xAOD::JetFourMom_t jet4vec(jet.Pt(),jet.Eta(),jet.Phi(),jet.M());
			
			(Event->GetAddJets().at(0))->setJetP4(jet4vec);
		
			PT(*(Event->GetAddJets().at(0)))  = Event->GetFitAddJet().Pt();
			ETA(*(Event->GetAddJets().at(0))) = Event->GetFitAddJet().Eta();
			PHI(*(Event->GetAddJets().at(0))) = Event->GetFitAddJet().Phi();
			M(*(Event->GetAddJets().at(0)))   = Event->GetFitAddJet().M();
			isB(*(Event->GetAddJets().at(0))) = false;
		}
	}else{

		TLorentzVector jet1;
		TLorentzVector jet2;

		jet1.SetPtEtaPhiE((Event->GetBJetOne())->pt(),(Event->GetBJetOne())->eta(),(Event->GetBJetOne())->phi(),(Event->GetBJetOne())->e());
		jet2.SetPtEtaPhiE((Event->GetBJetTwo())->pt(),(Event->GetBJetTwo())->eta(),(Event->GetBJetTwo())->phi(),(Event->GetBJetTwo())->e());

		xAOD::JetFourMom_t jet14vec(jet1.Pt(),jet1.Eta(),jet1.Phi(),jet1.M());
		xAOD::JetFourMom_t jet24vec(jet2.Pt(),jet2.Eta(),jet2.Phi(),jet2.M());

		(Event->GetBJetOne())->setJetP4(jet14vec);
		(Event->GetBJetTwo())->setJetP4(jet24vec);

		PT(*(Event->GetBJetOne()))  = Event->GetBJetOne()->pt();
		ETA(*(Event->GetBJetOne())) = Event->GetBJetOne()->eta();
		PHI(*(Event->GetBJetOne())) = Event->GetBJetOne()->phi();
		M(*(Event->GetBJetOne()))   = Event->GetBJetOne()->m();
		PT(*(Event->GetBJetTwo()))  = Event->GetBJetTwo()->pt();
		ETA(*(Event->GetBJetTwo())) = Event->GetBJetTwo()->eta();
		PHI(*(Event->GetBJetTwo())) = Event->GetBJetTwo()->phi();
		M(*(Event->GetBJetTwo()))   = Event->GetBJetTwo()->m();

		PT(*(Event->GetPhotonOne()))  = Event->GetPhotonOne()->pt();
		ETA(*(Event->GetPhotonOne())) = Event->GetPhotonOne()->eta();
		PHI(*(Event->GetPhotonOne())) = Event->GetPhotonOne()->phi();
		PT(*(Event->GetPhotonTwo()))  = Event->GetPhotonTwo()->pt();
		ETA(*(Event->GetPhotonTwo())) = Event->GetPhotonTwo()->eta();
		PHI(*(Event->GetPhotonTwo())) = Event->GetPhotonTwo()->phi();
		

	}

}



