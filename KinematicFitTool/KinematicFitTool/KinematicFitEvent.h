///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2019-2020 CERN for the benefit of the ATLAS collaboration
*/

// KinematicFitEvent.h
// header file for class KinematicFitEvent
// Author: BELFKIR Mohamed 
// Email : mohamed.belfkir@cern.ch
///////////////////////////////////////////////////////////////////

#ifndef KinematicFitEvent_H
#define KinematicFitEvent_H 1

#include <string.h>

#include <TString.h>
#include <TEnv.h>

// KinematicFitEvent includes
#include "AsgTools/AsgTool.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/Photon.h"
#include "xAODEgamma/PhotonContainer.h"

class KinematicFitEvent : public asg::AsgTool {

	
public:

  /// Constructor with parameter name: 
  KinematicFitEvent();

  /// Destructor: 
  virtual            ~KinematicFitEvent(); 
  StatusCode         initialize();
  StatusCode         finalize();
  StatusCode applySelection(const xAOD::PhotonContainer& photons, xAOD::JetContainer& jets);
  const xAOD::Photon*        GetPhotonOne();
  const xAOD::Photon*        GetPhotonTwo();
   xAOD::Jet*           GetBJetTwo();
   xAOD::Jet*           GetBJetOne();
  std::vector<  xAOD::Jet*>  GetAddJets();
  TLorentzVector GetFitPhotonOne();
  TLorentzVector GetFitPhotonTwo();
  TLorentzVector GetFitBJetOne();
  TLorentzVector GetFitBJetTwo();
  TLorentzVector GetFitAddJet();
  void SetFitPhotonOne(TLorentzVector tlv);
  void SetFitPhotonTwo(TLorentzVector tlv);
  void SetFitBJetOne(TLorentzVector tlv);
  void SetFitBJetTwo(TLorentzVector tlv);
  void SetFitAddJet(TLorentzVector tlv);
  void Clear();

protected:


private:
  
 Double_t    m_Jet_Min_Pt;
 std::string m_BtaggingWP;

 const xAOD::Photon* m_photon1;
 const xAOD::Photon* m_photon2;

  xAOD::Jet*    m_bjet1;
  xAOD::Jet*    m_bjet2;
 std::vector<  xAOD::Jet*> m_addJ;

 TLorentzVector m_photon1_Fit;
 TLorentzVector m_photon2_Fit;

 TLorentzVector m_bjet1_Fit;
 TLorentzVector m_bjet2_Fit;
 TLorentzVector m_addJ_Fit;

//Private members
private:

}; 

#endif //> !KinematicFitEvent_H
