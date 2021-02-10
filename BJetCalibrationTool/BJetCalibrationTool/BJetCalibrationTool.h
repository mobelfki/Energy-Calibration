///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2019-2020 CERN for the benefit of the ATLAS collaboration
*/

// BJetCalibrationTool.cxx
// Source file for class BJetCalibrationTool
// Author: BELFKIR Mohamed 
// Email : mohamed.belfkir@cern.ch
///////////////////////////////////////////////////////////////////

#ifndef BJetCalibrationTool_APPLYBJETCALIBRATION_H
#define BJetCalibrationTool_APPLYBJETCALIBRATION_H 1

#include <string.h>

#include <TString.h>
#include <TEnv.h>

#include "AsgTools/AsgTool.h"

// BJetCalibrationTool includes
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"
#include "MuonSelectorTools/MuonSelectionTool.h"

class BJetCalibrationTool : public asg::AsgTool {

	ASG_TOOL_CLASS0(BJetCalibrationTool)

public:

  /// Constructor with parameter name: 
  BJetCalibrationTool(const std::string& name = "BJetCalibrationTool");

  /// Destructor: 
  virtual ~BJetCalibrationTool(); 
  virtual StatusCode initializeTool(const std::string& name);
  StatusCode initialize();
  StatusCode finalize();
  virtual StatusCode applyBJetCalibration(xAOD::Jet& jet);
  StatusCode initializeMuonContainer(std::vector< const xAOD::Muon* >& muon);
protected:
  /// This is where the actual calibration code goes.
 

//Private methods
private:
  
  std::vector< const xAOD::Muon* > getSelectedMuon(const xAOD::MuonContainer* &MuonContainer);
  std::vector< const xAOD::Muon* > getMuonInJet(xAOD::Jet& jet, std::vector< const xAOD::Muon* > muons);
  void addMuon(xAOD::Jet& jet, std::vector< const xAOD::Muon* > muons);
  const xAOD::Muon* getMuon(std::vector< const xAOD::Muon* > muons);
  void applyPtReco(xAOD::Jet& jet);

//Private members
private:

  //Variables for configuration
  std::string m_JetAlgo;
  Double_t    m_Jet_Min_Pt;
  Double_t    m_Jet_Min_Eta;
  
  std::string m_MuonContainer_Name;
  xAOD::Muon::Quality m_Muon_Quality;
  Double_t    m_Muon_Min_Pt;
  Double_t    m_Muon_Jet_dR;
  bool        m_doVR;
  bool 	      m_doPtReco;
  
  std::string m_PtReco;
  TString m_BJet_Tagger;
 
  CP::MuonSelectionTool m_muonSelection;
  
  TFile *m_PtRecoFile;
  TH1F  *m_Semi_Histo;
  TH1F  *m_Had_Histo;
 
 
}; 

#endif //> !BJetCalibrationTool_APPLYBJETCALIBRATION_H
