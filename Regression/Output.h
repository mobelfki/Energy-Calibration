//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 21 15:43:53 2019 by ROOT version 6.12/06
// from TTree Output/
// found on file: mc16_13TeV.root
//////////////////////////////////////////////////////////

#ifndef Output_h
#define Output_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
double Etabins [] = {-2.5,-2.4,-2.3,-2.2,-2.1,-2.,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5};
const int NEta = 50; 
TH2F* R_vs_Truth[NEta];
TH2F* R_vs_Reco[NEta];
TProfile* Truth_vs_Truth[NEta];
TProfile* Truth_vs_Reco[NEta];
TProfile* Reco_vs_Reco[NEta];
TFile* file = new TFile("output/New_Topo_DL1r_FixedCut60.root","RECREATE");

/*
TFile* ptReco = new TFile("AntiKt4EMTopo_PtReco_Correction_10092019.root","READ");

TH1F* Corr_Semi = (TH1F*) ptReco->Get("Correction_SemiLeptonic_ttbar_mean");
TH1F* Corr_Had  = (TH1F*) ptReco->Get("Correction_Hadronic_ttbar_mean");
TH1F* Corr_Inc  = (TH1F*) ptReco->Get("Correction_Inclusive_ttbar_mean");
*/

class Output {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
	TLorentzVector addMuon(TLorentzVector jet);
	TLorentzVector addElectron(TLorentzVector jet);
	std::vector<TLorentzVector> matchTruth(TLorentzVector jet);
	void FillHistos(TLorentzVector truth, TLorentzVector reco , double weight = 1);
	TGraphErrors* constructGraph(TH2F* _2DHistos, TProfile* _AvrgE, bool doReco = false );
	double getFitMin(TH1D *inputHist, double gausFitNrms);
	double getFitMax(TH1D *inputHist, double gausFitNrms);
	void FitGauss(TH1D *R, TF1 *_gauss, Double_t fitMin, Double_t fitMax, bool adjustRange);
	TH1D* project1D(TH2F *h, int low_bin, int high_bin);
	void getAvg(TProfile *inputProfile, int low, int high, double &avg, double &avg_err);
	std::vector<int> getBinList (TH2F *hist, bool doReco);
	double getNEff(TH1 *h);
	void Loop2(std::vector<TF1*> corrfunc);
	void  FillRecoHistos(TLorentzVector truth, TLorentzVector reco, TLorentzVector recoest, double weight = 1);
	double getTruthFactor(double E, double eta, std::vector<TF1*> corrfunc);
	TF1* FitGraph(TGraphErrors* g, int npar);
	std::vector<TF1*> FitGraphs(std::vector<TGraphErrors*> graphs);
	void constructHisto(std::vector<TGraphErrors*> graphs);
	double getEta(TH2F *hist,bool doReco);
	std::vector<TString> genFunctionList();
	double getRecoFactor(double E, double eta, std::vector<TF1*> corrfunc);
	void Loop3(std::vector<TF1*> corrfunc);
   // Declaration of leaf types
   vector<float>   *Truth_WZJet_Pt;
   vector<float>   *Truth_WZJet_Eta;
   vector<float>   *Truth_WZJet_Phi;
   vector<float>   *Truth_WZJet_E;
   vector<float>   *Truth_WZJet_ConeTruthLabelID;
   Float_t         Truth_WZJet_Size;
   vector<float>   *Reco_Jet_Pt;
   vector<float>   *Reco_Jet_Eta;
   vector<float>   *Reco_Jet_Phi;
   vector<float>   *Reco_Jet_E;

   vector<int>     *Reco_Jet_MV2c10_bin;
   vector<int>     *Reco_Jet_DL1_bin;

   Float_t         Reco_Jet_Size;
   vector<float>   *Reco_Muon_Pt;
   vector<float>   *Reco_Muon_Eta;
   vector<float>   *Reco_Muon_Phi;
   vector<float>   *Reco_Muon_E;
   vector<float>   *Reco_Muon_ELoss_Pt;
   vector<float>   *Reco_Muon_ELoss_Eta;
   vector<float>   *Reco_Muon_ELoss_Phi;
   vector<float>   *Reco_Muon_ELoss_E;
   Float_t         Reco_Muon_Size;
   vector<float>   *Reco_Electron_Pt;
   vector<float>   *Reco_Electron_Eta;
   vector<float>   *Reco_Electron_Phi;
   vector<float>   *Reco_Electron_E;
   Float_t         Reco_Electron_Size;
   Float_t         PU;
   Float_t         mc_weight;
   Float_t         pu_weight;
   Float_t         lumi_weight;

   // List of branches
   TBranch        *b_Truth_WZJet_Pt;   //!
   TBranch        *b_Truth_WZJet_Eta;   //!
   TBranch        *b_Truth_WZJet_Phi;   //!
   TBranch        *b_Truth_WZJet_E;   //!
   TBranch        *b_Truth_WZJet_ConeTruthLabelID;   //!
   TBranch        *b_Truth_WZJet_Size;   //!
   TBranch        *b_Reco_Jet_Pt;   //!
   TBranch        *b_Reco_Jet_Eta;   //!
   TBranch        *b_Reco_Jet_Phi;   //!
   TBranch        *b_Reco_Jet_E;   //!
   TBranch        *b_Reco_Jet_MV2c10_bin;   //!
   TBranch        *b_Reco_Jet_DL1_bin;   //!
   TBranch        *b_Reco_Jet_Size;   //!
   TBranch        *b_Reco_Muon_Pt;   //!
   TBranch        *b_Reco_Muon_Eta;   //!
   TBranch        *b_Reco_Muon_Phi;   //!
   TBranch        *b_Reco_Muon_E;   //!
   TBranch        *b_Reco_Muon_ELoss_Pt;   //!
   TBranch        *b_Reco_Muon_ELoss_Eta;   //!
   TBranch        *b_Reco_Muon_ELoss_Phi;   //!
   TBranch        *b_Reco_Muon_ELoss_E;   //!
   TBranch        *b_Reco_Muon_Size;   //!
   TBranch        *b_Reco_Electron_Pt;   //!
   TBranch        *b_Reco_Electron_Eta;   //!
   TBranch        *b_Reco_Electron_Phi;   //!
   TBranch        *b_Reco_Electron_E;   //!
   TBranch        *b_Reco_Electron_Size;   //!
   TBranch        *b_PU;   //!
   TBranch        *b_mc_weight;   //!
   TBranch        *b_pu_weight;   //!
   TBranch        *b_lumi_weight;   //!

   Output(TTree *tree=0);
   virtual ~Output();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(vector<TFile*> files);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif

#ifdef Output_cxx
Output::Output(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
	vector<TFile*> files;
   if (tree == 0) {
      TFile *f_1 = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/user/m/mobelfki/user.mobelfki.mc16_13TeV.410470.ttbar_topo_TaggerBin.Ntuple.e6337_e5984_s3126_r10724_r10726_p3759_H4D4_ANALYSIS.root/user.mobelfki.21569473._000001.ANALYSIS.root");
     
      if (!f_1 || !f_1->IsOpen()) {
         f_1 = new TFile("/eos/user/m/mobelfki/user.mobelfki.mc16_13TeV.410470.ttbar_topo_TaggerBin.Ntuple.e6337_e5984_s3126_r10724_r10726_p3759_H4D4_ANALYSIS.root/user.mobelfki.21569473._000001.ANALYSIS.root");
      }
   //   f->GetObject("Output",tree);
   files.push_back(f_1);
   }
   Init(files);
}

Output::~Output()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Output::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Output::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Output::Init(vector<TFile*> files)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Truth_WZJet_Pt = 0;
   Truth_WZJet_Eta = 0;
   Truth_WZJet_Phi = 0;
   Truth_WZJet_E = 0;
   Truth_WZJet_ConeTruthLabelID = 0;
   Reco_Jet_Pt = 0;
   Reco_Jet_Eta = 0;
   Reco_Jet_Phi = 0;
   Reco_Jet_E = 0;
   Reco_Muon_Pt = 0;
   Reco_Muon_Eta = 0;
   Reco_Muon_Phi = 0;
   Reco_Muon_E = 0;
   Reco_Muon_ELoss_Pt = 0;
   Reco_Muon_ELoss_Eta = 0;
   Reco_Muon_ELoss_Phi = 0;
   Reco_Muon_ELoss_E = 0;
   Reco_Electron_Pt = 0;
   Reco_Electron_Eta = 0;
   Reco_Electron_Phi = 0;
   Reco_Electron_E = 0;
   // Set branch addresses and branch pointers
  // if (!tree) return;
   fChain = new TChain("Output");


   fChain->AddFile(files[0]->GetName(),-1,"Output");
 //  fChain->AddFile(files[1]->GetName(),-1,"Output");
  // fChain->AddFile(files[2]->GetName(),-1,"Output");



   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Truth_WZJet_Pt", &Truth_WZJet_Pt, &b_Truth_WZJet_Pt);
   fChain->SetBranchAddress("Truth_WZJet_Eta", &Truth_WZJet_Eta, &b_Truth_WZJet_Eta);
   fChain->SetBranchAddress("Truth_WZJet_Phi", &Truth_WZJet_Phi, &b_Truth_WZJet_Phi);
   fChain->SetBranchAddress("Truth_WZJet_E", &Truth_WZJet_E, &b_Truth_WZJet_E);
   fChain->SetBranchAddress("Truth_WZJet_ConeTruthLabelID", &Truth_WZJet_ConeTruthLabelID, &b_Truth_WZJet_ConeTruthLabelID);
   fChain->SetBranchAddress("Truth_WZJet_Size", &Truth_WZJet_Size, &b_Truth_WZJet_Size);
   fChain->SetBranchAddress("Reco_Jet_Pt", &Reco_Jet_Pt, &b_Reco_Jet_Pt);
   fChain->SetBranchAddress("Reco_Jet_Eta", &Reco_Jet_Eta, &b_Reco_Jet_Eta);
   fChain->SetBranchAddress("Reco_Jet_Phi", &Reco_Jet_Phi, &b_Reco_Jet_Phi);
   fChain->SetBranchAddress("Reco_Jet_E", &Reco_Jet_E, &b_Reco_Jet_E);
   fChain->SetBranchAddress("Reco_Jet_MV2c10_bin", &Reco_Jet_MV2c10_bin, &b_Reco_Jet_MV2c10_bin);
   fChain->SetBranchAddress("Reco_Jet_DL1_bin", &Reco_Jet_DL1_bin, &b_Reco_Jet_DL1_bin);
   fChain->SetBranchAddress("Reco_Jet_Size", &Reco_Jet_Size, &b_Reco_Jet_Size);
   fChain->SetBranchAddress("Reco_Muon_Pt", &Reco_Muon_Pt, &b_Reco_Muon_Pt);
   fChain->SetBranchAddress("Reco_Muon_Eta", &Reco_Muon_Eta, &b_Reco_Muon_Eta);
   fChain->SetBranchAddress("Reco_Muon_Phi", &Reco_Muon_Phi, &b_Reco_Muon_Phi);
   fChain->SetBranchAddress("Reco_Muon_E", &Reco_Muon_E, &b_Reco_Muon_E);
   fChain->SetBranchAddress("Reco_Muon_ELoss_Pt", &Reco_Muon_ELoss_Pt, &b_Reco_Muon_ELoss_Pt);
   fChain->SetBranchAddress("Reco_Muon_ELoss_Eta", &Reco_Muon_ELoss_Eta, &b_Reco_Muon_ELoss_Eta);
   fChain->SetBranchAddress("Reco_Muon_ELoss_Phi", &Reco_Muon_ELoss_Phi, &b_Reco_Muon_ELoss_Phi);
   fChain->SetBranchAddress("Reco_Muon_ELoss_E", &Reco_Muon_ELoss_E, &b_Reco_Muon_ELoss_E);
   fChain->SetBranchAddress("Reco_Muon_Size", &Reco_Muon_Size, &b_Reco_Muon_Size);
   fChain->SetBranchAddress("Reco_Electron_Pt", &Reco_Electron_Pt, &b_Reco_Electron_Pt);
   fChain->SetBranchAddress("Reco_Electron_Eta", &Reco_Electron_Eta, &b_Reco_Electron_Eta);
   fChain->SetBranchAddress("Reco_Electron_Phi", &Reco_Electron_Phi, &b_Reco_Electron_Phi);
   fChain->SetBranchAddress("Reco_Electron_E", &Reco_Electron_E, &b_Reco_Electron_E);
   fChain->SetBranchAddress("Reco_Electron_Size", &Reco_Electron_Size, &b_Reco_Electron_Size);
   fChain->SetBranchAddress("PU", &PU, &b_PU);
   fChain->SetBranchAddress("mc_weight", &mc_weight, &b_mc_weight);
   fChain->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
   fChain->SetBranchAddress("lumi_weight", &lumi_weight, &b_lumi_weight);
   Notify();
}

Bool_t Output::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Output::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Output::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Output_cxx
