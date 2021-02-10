****BJetCalibrationTool

Definition :
    - The BJetCalibrationTool is a tool to apply calibration to B-Jet, this tool is applied after JES and B-tagging to any B-Jet.
    
Methode    :
    - The tool use applyBJetCorrection methode to apply corrections to B-Jet. The first correction is adding muons reconstructed inside the jet cone back to jet, only one muons with highest pT is added to jet if more than one muon founded in the jet. After the tool apply a scale factor to jet depends on the jet pT and jet type (Hadronic (0 muon founded) or Semileptonic (1 or more muons founded during adding muons stage)).
    - The jet will be decorated with two new variables an integer which is the number of muons inside the jet called "n_mouns" and a float number which is the scale factor applied to the jet called "PtReco_SF".

Initializtion :

    - Check the TWiki page : https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BJetCorrectionsHowTo
    
