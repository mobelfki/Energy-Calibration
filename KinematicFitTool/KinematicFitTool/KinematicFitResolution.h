///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2019-2020 CERN for the benefit of the ATLAS collaboration
*/

// KinematicFitResolution.h
// header file for class KinematicFitResolution
// Author: BELFKIR Mohamed 
// Email : mohamed.belfkir@cern.ch
///////////////////////////////////////////////////////////////////

#ifndef KinematicFitResolution_H
#define KinematicFitResolution_H 1


// KinematicFitResolution includes
#include <string.h>
#include <TString.h>
#include <TEnv.h>
#include "AsgTools/AsgTool.h"
#include <TH1F.h>
#include <TFile.h>
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/Photon.h"
#include "xAODEgamma/PhotonContainer.h"

class KinematicFitResolution : public asg::AsgTool {

public:

  /// Constructor with parameter name: 
  KinematicFitResolution();

  /// Destructor: 
  virtual            ~KinematicFitResolution(); 
  StatusCode initialize();

  double PhotonRes();
  double AnglesRes();
  bool   FixAngles();
  double ConstraintValue(Int_t n);
  double JetRes( xAOD::Jet* Jet);
  double JetMean( xAOD::Jet* Jet);
  double JetResponse(TLorentzVector FitJet,  xAOD::Jet* Jet);

protected:
  
private:

private:

  //Variables for configuration

  double m_photon_Res;
  double m_angles_Res;
  bool   m_isFixAngles;
  double m_2jet_constraint;
  double m_3jet_constraint;
  Double_t    m_Jet_Min_Pt;
  std::string m_BtaggingWP;
  TString m_EnergyResolutionName;
  TString m_EnergyResponseName;
  TFile* m_EnergyResolution;
  TFile* m_EnergyResponse;

 
}; 

#endif //> !KinematicFitResolution_H
