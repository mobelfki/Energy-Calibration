# The name of the package:
    cmake_minimum_required(VERSION 3.11)

atlas_subdir (KinematicFitTool)

find_package( ROOT COMPONENTS Core Tree MathCore Hist Physics RIO Minuit )

# Add the shared library:
atlas_add_library (KinematicFitToolLib
  KinematicFitTool/*.h Root/*.cxx ${_dictionarySource}
  PUBLIC_HEADERS KinematicFitTool
  LINK_LIBRARIES AnaAlgorithmLib xAODEventInfo xAODJet PathResolver xAODEgamma xAODJet ${ROOT_LIBRARIES}
  INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
)

if (XAOD_STANDALONE)
 # Add the dictionary (for AnalysisBase only):
 atlas_add_dictionary (KinematicFitToolDict
  KinematicFitTool/KinematicFitToolDict.h
  KinematicFitTool/selection.xml
  LINK_LIBRARIES KinematicFitToolLib MuonSelectorToolsLib)
endif ()

if (NOT XAOD_STANDALONE)
  # Add a component library for AthAnalysis only:
  atlas_add_component (KinematicFitTool
    src/components/*.cxx
    LINK_LIBRARIES KinematicFitToolLib )
endif ()

# Install files from the package:
atlas_install_joboptions( share/*.py )
atlas_install_data( data/* )
