///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2019-2020 CERN for the benefit of the ATLAS collaboration
*/

// KinematicFitTool.h
// header file for class KinematicFitTool
// Author: BELFKIR Mohamed 
// Email : mohamed.belfkir@cern.ch
///////////////////////////////////////////////////////////////////

#ifndef KinematicFitTool_H
#define KinematicFitTool_H 1

// KinematicFitTool includes
#include <string.h>
#include <TString.h>
#include <TEnv.h>
#include "AsgTools/AsgTool.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/Photon.h"
#include "xAODEgamma/PhotonContainer.h"
#include "KinematicFitTool/KinematicFitEvent.h"
#include "KinematicFitTool/KinematicFitResolution.h"
#include "KinematicFitTool/KinematicFitRun.h"
#include "TLorentzVector.h"

class KinematicFitTool : public asg::AsgTool {

	ASG_TOOL_CLASS0(KinematicFitTool)

public:

  /// Constructor with parameter name: 
  KinematicFitTool(const std::string& name = "KinematicFitTool");

  /// Destructor: 
  virtual ~KinematicFitTool(); 
  virtual StatusCode initializeTool(const std::string& name);
  StatusCode initialize();
  StatusCode finalize();
  StatusCode applyKF(const xAOD::PhotonContainer& photons, xAOD::JetContainer& jets);
  
  inline std::string GetBTaggingWP()
	{return m_BtaggingWP;}
  inline double GetMinJetPt()
	{return m_Jet_Min_Pt;}
	
  
protected:

   void DecorateEvent(KinematicFitEvent*& Event, KinematicFitRun*& Run);

private:
  
  enum FitType {TwoJet = 0, ThreeJet = 1};
  const int GeV     = 1e3;
  KinematicFitResolution* fResolution;
  KinematicFitEvent* fEvent;
  KinematicFitRun* fProcessor;

  std::string m_JetAlgo;
  std::string m_file_name;

  double m_photon_Res;
  double m_angles_Res;
  bool   m_isFixAngles;
  double m_2jet_constraint;
  double m_3jet_constraint;
  Double_t    m_Jet_Min_Pt;
  std::string m_BtaggingWP;

}; 

#endif //> !KinematicFitTool_H
