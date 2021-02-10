///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2019-2020 CERN for the benefit of the ATLAS collaboration
*/

// KinematicFitRun.h
// header file for class KinematicFitRun
// Author: BELFKIR Mohamed 
// Email : mohamed.belfkir@cern.ch
///////////////////////////////////////////////////////////////////

#ifndef KinematicFitRun_H
#define KinematicFitRun_H 1

// KinematicFitRun includes
#include <string.h>
#include <TString.h>
#include <TEnv.h>
#include "AsgTools/AsgTool.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/Photon.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/Photon.h"
#include "xAODEgamma/PhotonContainer.h"
#include "TMinuit.h"
#include "KinematicFitTool/KinematicFitEvent.h"
#include "KinematicFitTool/KinematicFitResolution.h"

class KinematicFitRun  {

public:

  /// Constructor with parameter name: 
  KinematicFitRun();

  /// Destructor: 
  virtual            ~KinematicFitRun(); 
  StatusCode         initialize();
  StatusCode         finalize();
  StatusCode RunKF(KinematicFitEvent*& Event, KinematicFitResolution*& Resolution);
  Bool_t isFitted();
protected:

	StatusCode SetEventAndResolution(KinematicFitEvent*& Event, KinematicFitResolution*& Resolution);
	void Run();
	void SetParameterNameandRange(TMinuit*& m);
	static void FCN(int& /*npar*/, double* /*grad*/, double& fval, double* par, int /*flag*/);
	double GetConstraintRes(std::vector<  xAOD::Jet*> AddJets);
	double GetJetRes( xAOD::Jet* Jet);
	double GetResponse(TLorentzVector FitJet,  xAOD::Jet* Jet);
	double GetJetMean( xAOD::Jet* Jet);
	Int_t GetNParameter();
	void NewEvent();
	void ShareFit(TMinuit*& m);
	double GetLH();
	void SetParameters(double* par);
	
private:

  enum FitType {TwoJet = 0, ThreeJet = 1};
  const int GeV     = 1e3;
  Int_t g = 2;
  Bool_t isfitted = false;
  TLorentzVector _FitPhoton1;
  TLorentzVector _FitPhoton2;
  TLorentzVector _FitBJet1;
  TLorentzVector _FitBJet2;
  TLorentzVector _FitJet3;
  KinematicFitEvent* _Event;
  KinematicFitResolution* _Resolution;
  TMinuit *minuit;
  
 
}; 

#endif //> !KinematicFitRun_H
