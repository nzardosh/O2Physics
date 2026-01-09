// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file jetFinderTemplateInstantiations.cxx
/// \brief Explicit Template Instantiation file for all JetFinder templates
///
/// This file contains explicit template instantiations for all JetFinder template classes
/// to improve compilation performance by compiling templates once instead of in every
/// translation unit that uses them.
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Florian Jonas <florian.jonas@cern.ch>

// Prevent the main function from being defined when including headers
#define O2_NO_WORKFLOW_MAIN

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/JetFinders/jetFinder.h"
#include "PWGJE/JetFinders/jetFinderHF.h"
#include "PWGJE/JetFinders/jetFinderHFHFBar.h"
#include "PWGJE/JetFinders/jetFinderV0.h"

using namespace o2;
using namespace o2::aod;

// ============================================================================
// Explicit Template Instantiations for JetFinderTask
// ============================================================================

// Charged jets
template struct JetFinderTask<o2::aod::ChargedJets, o2::aod::ChargedJetConstituents, o2::aod::ChargedEventWiseSubtractedJets, o2::aod::ChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderTask<o2::aod::ChargedMCDetectorLevelJets, o2::aod::ChargedMCDetectorLevelJetConstituents, o2::aod::ChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderTask<o2::aod::ChargedMCParticleLevelJets, o2::aod::ChargedMCParticleLevelJetConstituents, o2::aod::ChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::ChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// Neutral jets
template struct JetFinderTask<o2::aod::NeutralJets, o2::aod::NeutralJetConstituents, o2::aod::NeutralEventWiseSubtractedJets, o2::aod::NeutralEventWiseSubtractedJetConstituents>;
template struct JetFinderTask<o2::aod::NeutralMCDetectorLevelJets, o2::aod::NeutralMCDetectorLevelJetConstituents, o2::aod::NeutralMCDetectorLevelEventWiseSubtractedJets, o2::aod::NeutralMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderTask<o2::aod::NeutralMCParticleLevelJets, o2::aod::NeutralMCParticleLevelJetConstituents, o2::aod::NeutralMCParticleLevelEventWiseSubtractedJets, o2::aod::NeutralMCParticleLevelEventWiseSubtractedJetConstituents>;

// Full jets
template struct JetFinderTask<o2::aod::FullJets, o2::aod::FullJetConstituents, o2::aod::FullEventWiseSubtractedJets, o2::aod::FullEventWiseSubtractedJetConstituents>;
template struct JetFinderTask<o2::aod::FullMCDetectorLevelJets, o2::aod::FullMCDetectorLevelJetConstituents, o2::aod::FullMCDetectorLevelEventWiseSubtractedJets, o2::aod::FullMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderTask<o2::aod::FullMCParticleLevelJets, o2::aod::FullMCParticleLevelJetConstituents, o2::aod::FullMCParticleLevelEventWiseSubtractedJets, o2::aod::FullMCParticleLevelEventWiseSubtractedJetConstituents>;

// duplicates
template struct JetFinderTask<aod::Charged1Jets, aod::Charged1JetConstituents, aod::Charged1EventWiseSubtractedJets, aod::Charged1EventWiseSubtractedJetConstituents>;
template struct JetFinderTask<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents, aod::Charged1MCDetectorLevelEventWiseSubtractedJets, aod::Charged1MCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderTask<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents, aod::Charged1MCParticleLevelEventWiseSubtractedJets, aod::Charged1MCParticleLevelEventWiseSubtractedJetConstituents>;

// D0 instantiations
template struct JetFinderHFTask<o2::aod::CandidatesD0Data, o2::aod::CandidatesD0MCD, o2::aod::CandidatesD0MCP, o2::aod::JetTracksSubD0, o2::aod::JetParticlesSubD0, o2::aod::D0ChargedJets, o2::aod::D0ChargedJetConstituents, o2::aod::D0ChargedEventWiseSubtractedJets, o2::aod::D0ChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesD0Data, o2::aod::CandidatesD0MCD, o2::aod::CandidatesD0MCP, o2::aod::JetTracksSubD0, o2::aod::JetParticlesSubD0, o2::aod::D0ChargedMCDetectorLevelJets, o2::aod::D0ChargedMCDetectorLevelJetConstituents, o2::aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesD0Data, o2::aod::CandidatesD0MCD, o2::aod::CandidatesD0MCP, o2::aod::JetTracksSubD0, o2::aod::JetParticlesSubD0, o2::aod::D0ChargedMCParticleLevelJets, o2::aod::D0ChargedMCParticleLevelJetConstituents, o2::aod::D0ChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::D0ChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// Dplus instantiations
template struct JetFinderHFTask<o2::aod::CandidatesDplusData, o2::aod::CandidatesDplusMCD, o2::aod::CandidatesDplusMCP, o2::aod::JetTracksSubDplus, o2::aod::JetParticlesSubDplus, o2::aod::DplusChargedJets, o2::aod::DplusChargedJetConstituents, o2::aod::DplusChargedEventWiseSubtractedJets, o2::aod::DplusChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesDplusData, o2::aod::CandidatesDplusMCD, o2::aod::CandidatesDplusMCP, o2::aod::JetTracksSubDplus, o2::aod::JetParticlesSubDplus, o2::aod::DplusChargedMCDetectorLevelJets, o2::aod::DplusChargedMCDetectorLevelJetConstituents, o2::aod::DplusChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DplusChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesDplusData, o2::aod::CandidatesDplusMCD, o2::aod::CandidatesDplusMCP, o2::aod::JetTracksSubDplus, o2::aod::JetParticlesSubDplus, o2::aod::DplusChargedMCParticleLevelJets, o2::aod::DplusChargedMCParticleLevelJetConstituents, o2::aod::DplusChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::DplusChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// Ds instantiations
template struct JetFinderHFTask<o2::aod::CandidatesDsData, o2::aod::CandidatesDsMCD, o2::aod::CandidatesDsMCP, o2::aod::JetTracksSubDs, o2::aod::JetParticlesSubDs, o2::aod::DsChargedJets, o2::aod::DsChargedJetConstituents, o2::aod::DsChargedEventWiseSubtractedJets, o2::aod::DsChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesDsData, o2::aod::CandidatesDsMCD, o2::aod::CandidatesDsMCP, o2::aod::JetTracksSubDs, o2::aod::JetParticlesSubDs, o2::aod::DsChargedMCDetectorLevelJets, o2::aod::DsChargedMCDetectorLevelJetConstituents, o2::aod::DsChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DsChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesDsData, o2::aod::CandidatesDsMCD, o2::aod::CandidatesDsMCP, o2::aod::JetTracksSubDs, o2::aod::JetParticlesSubDs, o2::aod::DsChargedMCParticleLevelJets, o2::aod::DsChargedMCParticleLevelJetConstituents, o2::aod::DsChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::DsChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// Dstar instantiations
template struct JetFinderHFTask<o2::aod::CandidatesDstarData, o2::aod::CandidatesDstarMCD, o2::aod::CandidatesDstarMCP, o2::aod::JetTracksSubDstar, o2::aod::JetParticlesSubDstar, o2::aod::DstarChargedJets, o2::aod::DstarChargedJetConstituents, o2::aod::DstarChargedEventWiseSubtractedJets, o2::aod::DstarChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesDstarData, o2::aod::CandidatesDstarMCD, o2::aod::CandidatesDstarMCP, o2::aod::JetTracksSubDstar, o2::aod::JetParticlesSubDstar, o2::aod::DstarChargedMCDetectorLevelJets, o2::aod::DstarChargedMCDetectorLevelJetConstituents, o2::aod::DstarChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DstarChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesDstarData, o2::aod::CandidatesDstarMCD, o2::aod::CandidatesDstarMCP, o2::aod::JetTracksSubDstar, o2::aod::JetParticlesSubDstar, o2::aod::DstarChargedMCParticleLevelJets, o2::aod::DstarChargedMCParticleLevelJetConstituents, o2::aod::DstarChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::DstarChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// Lc instantiations
template struct JetFinderHFTask<o2::aod::CandidatesLcData, o2::aod::CandidatesLcMCD, o2::aod::CandidatesLcMCP, o2::aod::JetTracksSubLc, o2::aod::JetParticlesSubLc, o2::aod::LcChargedJets, o2::aod::LcChargedJetConstituents, o2::aod::LcChargedEventWiseSubtractedJets, o2::aod::LcChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesLcData, o2::aod::CandidatesLcMCD, o2::aod::CandidatesLcMCP, o2::aod::JetTracksSubLc, o2::aod::JetParticlesSubLc, o2::aod::LcChargedMCDetectorLevelJets, o2::aod::LcChargedMCDetectorLevelJetConstituents, o2::aod::LcChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::LcChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesLcData, o2::aod::CandidatesLcMCD, o2::aod::CandidatesLcMCP, o2::aod::JetTracksSubLc, o2::aod::JetParticlesSubLc, o2::aod::LcChargedMCParticleLevelJets, o2::aod::LcChargedMCParticleLevelJetConstituents, o2::aod::LcChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::LcChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// B0 instantiations
template struct JetFinderHFTask<o2::aod::CandidatesB0Data, o2::aod::CandidatesB0MCD, o2::aod::CandidatesB0MCP, o2::aod::JetTracksSubB0, o2::aod::JetParticlesSubB0, o2::aod::B0ChargedJets, o2::aod::B0ChargedJetConstituents, o2::aod::B0ChargedEventWiseSubtractedJets, o2::aod::B0ChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesB0Data, o2::aod::CandidatesB0MCD, o2::aod::CandidatesB0MCP, o2::aod::JetTracksSubB0, o2::aod::JetParticlesSubB0, o2::aod::B0ChargedMCDetectorLevelJets, o2::aod::B0ChargedMCDetectorLevelJetConstituents, o2::aod::B0ChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::B0ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesB0Data, o2::aod::CandidatesB0MCD, o2::aod::CandidatesB0MCP, o2::aod::JetTracksSubB0, o2::aod::JetParticlesSubB0, o2::aod::B0ChargedMCParticleLevelJets, o2::aod::B0ChargedMCParticleLevelJetConstituents, o2::aod::B0ChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::B0ChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// Bplus instantiations
template struct JetFinderHFTask<o2::aod::CandidatesBplusData, o2::aod::CandidatesBplusMCD, o2::aod::CandidatesBplusMCP, o2::aod::JetTracksSubBplus, o2::aod::JetParticlesSubBplus, o2::aod::BplusChargedJets, o2::aod::BplusChargedJetConstituents, o2::aod::BplusChargedEventWiseSubtractedJets, o2::aod::BplusChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesBplusData, o2::aod::CandidatesBplusMCD, o2::aod::CandidatesBplusMCP, o2::aod::JetTracksSubBplus, o2::aod::JetParticlesSubBplus, o2::aod::BplusChargedMCDetectorLevelJets, o2::aod::BplusChargedMCDetectorLevelJetConstituents, o2::aod::BplusChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::BplusChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesBplusData, o2::aod::CandidatesBplusMCD, o2::aod::CandidatesBplusMCP, o2::aod::JetTracksSubBplus, o2::aod::JetParticlesSubBplus, o2::aod::BplusChargedMCParticleLevelJets, o2::aod::BplusChargedMCParticleLevelJetConstituents, o2::aod::BplusChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::BplusChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// XicToXiPiPi instantiations
template struct JetFinderHFTask<o2::aod::CandidatesXicToXiPiPiData, o2::aod::CandidatesXicToXiPiPiMCD, o2::aod::CandidatesXicToXiPiPiMCP, o2::aod::JetTracksSubXicToXiPiPi, o2::aod::JetParticlesSubXicToXiPiPi, o2::aod::XicToXiPiPiChargedJets, o2::aod::XicToXiPiPiChargedJetConstituents, o2::aod::XicToXiPiPiChargedEventWiseSubtractedJets, o2::aod::XicToXiPiPiChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesXicToXiPiPiData, o2::aod::CandidatesXicToXiPiPiMCD, o2::aod::CandidatesXicToXiPiPiMCP, o2::aod::JetTracksSubXicToXiPiPi, o2::aod::JetParticlesSubXicToXiPiPi, o2::aod::XicToXiPiPiChargedMCDetectorLevelJets, o2::aod::XicToXiPiPiChargedMCDetectorLevelJetConstituents, o2::aod::XicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::XicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesXicToXiPiPiData, o2::aod::CandidatesXicToXiPiPiMCD, o2::aod::CandidatesXicToXiPiPiMCP, o2::aod::JetTracksSubXicToXiPiPi, o2::aod::JetParticlesSubXicToXiPiPi, o2::aod::XicToXiPiPiChargedMCParticleLevelJets, o2::aod::XicToXiPiPiChargedMCParticleLevelJetConstituents, o2::aod::XicToXiPiPiChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::XicToXiPiPiChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// Dielectron instantiations
template struct JetFinderHFTask<o2::aod::CandidatesDielectronData, o2::aod::CandidatesDielectronMCD, o2::aod::CandidatesDielectronMCP, o2::aod::JetTracksSubDielectron, o2::aod::JetParticlesSubDielectron, o2::aod::DielectronChargedJets, o2::aod::DielectronChargedJetConstituents, o2::aod::DielectronChargedEventWiseSubtractedJets, o2::aod::DielectronChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesDielectronData, o2::aod::CandidatesDielectronMCD, o2::aod::CandidatesDielectronMCP, o2::aod::JetTracksSubDielectron, o2::aod::JetParticlesSubDielectron, o2::aod::DielectronChargedMCDetectorLevelJets, o2::aod::DielectronChargedMCDetectorLevelJetConstituents, o2::aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderHFTask<o2::aod::CandidatesDielectronData, o2::aod::CandidatesDielectronMCD, o2::aod::CandidatesDielectronMCP, o2::aod::JetTracksSubDielectron, o2::aod::JetParticlesSubDielectron, o2::aod::DielectronChargedMCParticleLevelJets, o2::aod::DielectronChargedMCParticleLevelJetConstituents, o2::aod::DielectronChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::DielectronChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// D0D0Bar instantiations
template struct JetFinderHFHFBarTask<o2::aod::CandidatesD0Data, o2::aod::CandidatesD0MCD, o2::aod::CandidatesD0MCP, o2::aod::JetTracksSubD0, o2::aod::JetParticlesSubD0, o2::aod::D0ChargedJets, o2::aod::D0ChargedJetConstituents, o2::aod::D0ChargedEventWiseSubtractedJets, o2::aod::D0ChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderHFHFBarTask<o2::aod::CandidatesD0Data, o2::aod::CandidatesD0MCD, o2::aod::CandidatesD0MCP, o2::aod::JetTracksSubD0, o2::aod::JetParticlesSubD0, o2::aod::D0ChargedMCDetectorLevelJets, o2::aod::D0ChargedMCDetectorLevelJetConstituents, o2::aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderHFHFBarTask<o2::aod::CandidatesD0Data, o2::aod::CandidatesD0MCD, o2::aod::CandidatesD0MCP, o2::aod::JetTracksSubD0, o2::aod::JetParticlesSubD0, o2::aod::D0ChargedMCParticleLevelJets, o2::aod::D0ChargedMCParticleLevelJetConstituents, o2::aod::D0ChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::D0ChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// DplusDminus instantiations
template struct JetFinderHFHFBarTask<o2::aod::CandidatesDplusData, o2::aod::CandidatesDplusMCD, o2::aod::CandidatesDplusMCP, o2::aod::JetTracksSubDplus, o2::aod::JetParticlesSubDplus, o2::aod::DplusChargedJets, o2::aod::DplusChargedJetConstituents, o2::aod::DplusChargedEventWiseSubtractedJets, o2::aod::DplusChargedEventWiseSubtractedJetConstituents>;
template struct JetFinderHFHFBarTask<o2::aod::CandidatesDplusData, o2::aod::CandidatesDplusMCD, o2::aod::CandidatesDplusMCP, o2::aod::JetTracksSubDplus, o2::aod::JetParticlesSubDplus, o2::aod::DplusChargedMCDetectorLevelJets, o2::aod::DplusChargedMCDetectorLevelJetConstituents, o2::aod::DplusChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DplusChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
template struct JetFinderHFHFBarTask<o2::aod::CandidatesDplusData, o2::aod::CandidatesDplusMCD, o2::aod::CandidatesDplusMCP, o2::aod::JetTracksSubDplus, o2::aod::JetParticlesSubDplus, o2::aod::DplusChargedMCParticleLevelJets, o2::aod::DplusChargedMCParticleLevelJetConstituents, o2::aod::DplusChargedMCParticleLevelEventWiseSubtractedJets, o2::aod::DplusChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

// V0 instantiations
template struct JetFinderV0Task<o2::aod::CandidatesV0Data, o2::aod::CandidatesV0MCD, o2::aod::CandidatesV0MCP, o2::aod::V0ChargedJets, o2::aod::V0ChargedJetConstituents>;
template struct JetFinderV0Task<o2::aod::CandidatesV0Data, o2::aod::CandidatesV0MCD, o2::aod::CandidatesV0MCP, o2::aod::V0ChargedMCDetectorLevelJets, o2::aod::V0ChargedMCDetectorLevelJetConstituents>;
template struct JetFinderV0Task<o2::aod::CandidatesV0Data, o2::aod::CandidatesV0MCD, o2::aod::CandidatesV0MCP, o2::aod::V0ChargedMCParticleLevelJets, o2::aod::V0ChargedMCParticleLevelJetConstituents>;
