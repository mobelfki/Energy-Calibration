# The name of the package:
    cmake_minimum_required(VERSION 3.11)

atlas_subdir (BJetCalibrationTool)

# Add the shared library:
atlas_add_library (BJetCalibrationToolLib
  BJetCalibrationTool/*.h Root/*.cxx
  PUBLIC_HEADERS BJetCalibrationTool
  LINK_LIBRARIES AnaAlgorithmLib xAODEventInfo xAODJet xAODMuon PathResolver MuonSelectorToolsLib)

if (XAOD_STANDALONE)
 # Add the dictionary (for AnalysisBase only):
 atlas_add_dictionary (BJetCalibrationToolDict
  BJetCalibrationTool/BJetCalibrationToolDict.h
  BJetCalibrationTool/selection.xml
  LINK_LIBRARIES BJetCalibrationToolLib MuonSelectorToolsLib)
endif ()

if (NOT XAOD_STANDALONE)
  # Add a component library for AthAnalysis only:
  atlas_add_component (BJetCalibrationTool
    src/components/*.cxx
    LINK_LIBRARIES BJetCalibrationToolLib MuonSelectorToolsLib)
endif ()

# Install files from the package:
atlas_install_joboptions( share/*.py )
atlas_install_data( data/* )
