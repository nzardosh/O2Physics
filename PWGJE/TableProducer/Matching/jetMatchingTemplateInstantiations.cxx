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

/// \file jetMatchingTemplateInstantiations.cxx
/// \brief Explicit Template Instantiation file for all JetMatching templates
///
/// This file contains explicit template instantiations for all JetMatching template classes
/// to improve compilation performance by compiling templates once instead of in every
/// translation unit that uses them.
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Florian Jonas <florian.jonas@cern.ch>

// Prevent the main function from being defined when including headers
#define O2_NO_WORKFLOW_MAIN

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/TableProducer/Matching/Duplicates/jetMatchingDuplicates.h"
#include "PWGJE/TableProducer/Matching/jetMatchingMC.h"
#include "PWGJE/TableProducer/Matching/jetMatchingMCSub.h"
#include "PWGJE/TableProducer/Matching/jetMatchingSub.h"

// ============================================================================
// Explicit Template Instantiations for JetMatching
// ============================================================================

template struct JetMatchingMc<o2::soa::Join<o2::aod::ChargedMCDetectorLevelJets, o2::aod::ChargedMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::ChargedMCParticleLevelJets, o2::aod::ChargedMCParticleLevelJetConstituents>,
                              o2::aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets,
                              o2::aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets,
                              o2::aod::JCollisions,
                              o2::aod::JMcCollisions,
                              o2::aod::JDummys>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::NeutralMCDetectorLevelJets, o2::aod::NeutralMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::NeutralMCParticleLevelJets, o2::aod::NeutralMCParticleLevelJetConstituents>,
                              o2::aod::NeutralMCDetectorLevelJetsMatchedToNeutralMCParticleLevelJets,
                              o2::aod::NeutralMCParticleLevelJetsMatchedToNeutralMCDetectorLevelJets,
                              o2::aod::JCollisions,
                              o2::aod::JMcCollisions,
                              o2::aod::JetClustersMCD>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::FullMCDetectorLevelJets, o2::aod::FullMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::FullMCParticleLevelJets, o2::aod::FullMCParticleLevelJetConstituents>,
                              o2::aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets,
                              o2::aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets,
                              o2::aod::JCollisions,
                              o2::aod::JMcCollisions,
                              o2::aod::JetClustersMCD>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::D0ChargedMCDetectorLevelJets, o2::aod::D0ChargedMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::D0ChargedMCParticleLevelJets, o2::aod::D0ChargedMCParticleLevelJetConstituents>,
                              o2::aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets,
                              o2::aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets,
                              o2::aod::CandidatesD0MCD,
                              o2::aod::CandidatesD0MCP,
                              o2::aod::JDummys>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::DplusChargedMCDetectorLevelJets, o2::aod::DplusChargedMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::DplusChargedMCParticleLevelJets, o2::aod::DplusChargedMCParticleLevelJetConstituents>,
                              o2::aod::DplusChargedMCDetectorLevelJetsMatchedToDplusChargedMCParticleLevelJets,
                              o2::aod::DplusChargedMCParticleLevelJetsMatchedToDplusChargedMCDetectorLevelJets,
                              o2::aod::CandidatesDplusMCD,
                              o2::aod::CandidatesDplusMCP,
                              o2::aod::JDummys>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::DsChargedMCDetectorLevelJets, o2::aod::DsChargedMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::DsChargedMCParticleLevelJets, o2::aod::DsChargedMCParticleLevelJetConstituents>,
                              o2::aod::DsChargedMCDetectorLevelJetsMatchedToDsChargedMCParticleLevelJets,
                              o2::aod::DsChargedMCParticleLevelJetsMatchedToDsChargedMCDetectorLevelJets,
                              o2::aod::CandidatesDsMCD,
                              o2::aod::CandidatesDsMCP,
                              o2::aod::JDummys>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::DstarChargedMCDetectorLevelJets, o2::aod::DstarChargedMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::DstarChargedMCParticleLevelJets, o2::aod::DstarChargedMCParticleLevelJetConstituents>,
                              o2::aod::DstarChargedMCDetectorLevelJetsMatchedToDstarChargedMCParticleLevelJets,
                              o2::aod::DstarChargedMCParticleLevelJetsMatchedToDstarChargedMCDetectorLevelJets,
                              o2::aod::CandidatesDstarMCD,
                              o2::aod::CandidatesDstarMCP,
                              o2::aod::JDummys>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::LcChargedMCDetectorLevelJets, o2::aod::LcChargedMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::LcChargedMCParticleLevelJets, o2::aod::LcChargedMCParticleLevelJetConstituents>,
                              o2::aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets,
                              o2::aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets,
                              o2::aod::CandidatesLcMCD,
                              o2::aod::CandidatesLcMCP,
                              o2::aod::JDummys>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::B0ChargedMCDetectorLevelJets, o2::aod::B0ChargedMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::B0ChargedMCParticleLevelJets, o2::aod::B0ChargedMCParticleLevelJetConstituents>,
                              o2::aod::B0ChargedMCDetectorLevelJetsMatchedToB0ChargedMCParticleLevelJets,
                              o2::aod::B0ChargedMCParticleLevelJetsMatchedToB0ChargedMCDetectorLevelJets,
                              o2::aod::CandidatesB0MCD,
                              o2::aod::CandidatesB0MCP,
                              o2::aod::JDummys>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::BplusChargedMCDetectorLevelJets, o2::aod::BplusChargedMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::BplusChargedMCParticleLevelJets, o2::aod::BplusChargedMCParticleLevelJetConstituents>,
                              o2::aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets,
                              o2::aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets,
                              o2::aod::CandidatesBplusMCD,
                              o2::aod::CandidatesBplusMCP,
                              o2::aod::JDummys>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::XicToXiPiPiChargedMCDetectorLevelJets, o2::aod::XicToXiPiPiChargedMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::XicToXiPiPiChargedMCParticleLevelJets, o2::aod::XicToXiPiPiChargedMCParticleLevelJetConstituents>,
                              o2::aod::XicToXiPiPiChargedMCDetectorLevelJetsMatchedToXicToXiPiPiChargedMCParticleLevelJets,
                              o2::aod::XicToXiPiPiChargedMCParticleLevelJetsMatchedToXicToXiPiPiChargedMCDetectorLevelJets,
                              o2::aod::CandidatesXicToXiPiPiMCD,
                              o2::aod::CandidatesXicToXiPiPiMCP,
                              o2::aod::JDummys>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::DielectronChargedMCDetectorLevelJets, o2::aod::DielectronChargedMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::DielectronChargedMCParticleLevelJets, o2::aod::DielectronChargedMCParticleLevelJetConstituents>,
                              o2::aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCParticleLevelJets,
                              o2::aod::DielectronChargedMCParticleLevelJetsMatchedToDielectronChargedMCDetectorLevelJets,
                              o2::aod::CandidatesDielectronMCD,
                              o2::aod::CandidatesDielectronMCP,
                              o2::aod::JDummys>;

template struct JetMatchingMc<o2::soa::Join<o2::aod::V0ChargedMCDetectorLevelJets, o2::aod::V0ChargedMCDetectorLevelJetConstituents>,
                              o2::soa::Join<o2::aod::V0ChargedMCParticleLevelJets, o2::aod::V0ChargedMCParticleLevelJetConstituents>,
                              o2::aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets,
                              o2::aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets,
                              o2::aod::CandidatesV0MCD,
                              o2::aod::CandidatesV0MCP,
                              o2::aod::JDummys>;

template struct JetMatchingMcSub<o2::soa::Join<o2::aod::ChargedMCDetectorLevelJets, o2::aod::ChargedMCDetectorLevelJetConstituents>,
                                 o2::soa::Join<o2::aod::ChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                 o2::aod::ChargedMCDetectorLevelJetsMatchedToChargedMCDetectorLevelEventWiseSubtractedJets,
                                 o2::aod::ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToChargedMCDetectorLevelJets,
                                 o2::aod::JDummys>;

template struct JetMatchingMcSub<o2::soa::Join<o2::aod::D0ChargedMCDetectorLevelJets, o2::aod::D0ChargedMCDetectorLevelJetConstituents>,
                                 o2::soa::Join<o2::aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                 o2::aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCDetectorLevelEventWiseSubtractedJets,
                                 o2::aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToD0ChargedMCDetectorLevelJets,
                                 o2::aod::CandidatesD0MCD>;

template struct JetMatchingMcSub<o2::soa::Join<o2::aod::DplusChargedMCDetectorLevelJets, o2::aod::DplusChargedMCDetectorLevelJetConstituents>,
                                 o2::soa::Join<o2::aod::DplusChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DplusChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                 o2::aod::DplusChargedMCDetectorLevelJetsMatchedToDplusChargedMCDetectorLevelEventWiseSubtractedJets,
                                 o2::aod::DplusChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToDplusChargedMCDetectorLevelJets,
                                 o2::aod::CandidatesDplusMCD>;

template struct JetMatchingMcSub<o2::soa::Join<o2::aod::DsChargedMCDetectorLevelJets, o2::aod::DsChargedMCDetectorLevelJetConstituents>,
                                 o2::soa::Join<o2::aod::DsChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DsChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                 o2::aod::DsChargedMCDetectorLevelJetsMatchedToDsChargedMCDetectorLevelEventWiseSubtractedJets,
                                 o2::aod::DsChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToDsChargedMCDetectorLevelJets,
                                 o2::aod::CandidatesDsMCD>;

template struct JetMatchingMcSub<o2::soa::Join<o2::aod::DstarChargedMCDetectorLevelJets, o2::aod::DstarChargedMCDetectorLevelJetConstituents>,
                                 o2::soa::Join<o2::aod::DstarChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DstarChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                 o2::aod::DstarChargedMCDetectorLevelJetsMatchedToDstarChargedMCDetectorLevelEventWiseSubtractedJets,
                                 o2::aod::DstarChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToDstarChargedMCDetectorLevelJets,
                                 o2::aod::CandidatesDstarMCD>;

template struct JetMatchingMcSub<o2::soa::Join<o2::aod::LcChargedMCDetectorLevelJets, o2::aod::LcChargedMCDetectorLevelJetConstituents>,
                                 o2::soa::Join<o2::aod::LcChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::LcChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                 o2::aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCDetectorLevelEventWiseSubtractedJets,
                                 o2::aod::LcChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToLcChargedMCDetectorLevelJets,
                                 o2::aod::CandidatesLcMCD>;

template struct JetMatchingMcSub<o2::soa::Join<o2::aod::B0ChargedMCDetectorLevelJets, o2::aod::B0ChargedMCDetectorLevelJetConstituents>,
                                 o2::soa::Join<o2::aod::B0ChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::B0ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                 o2::aod::B0ChargedMCDetectorLevelJetsMatchedToB0ChargedMCDetectorLevelEventWiseSubtractedJets,
                                 o2::aod::B0ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToB0ChargedMCDetectorLevelJets,
                                 o2::aod::CandidatesB0MCD>;

template struct JetMatchingMcSub<o2::soa::Join<o2::aod::BplusChargedMCDetectorLevelJets, o2::aod::BplusChargedMCDetectorLevelJetConstituents>,
                                 o2::soa::Join<o2::aod::BplusChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::BplusChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                 o2::aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCDetectorLevelEventWiseSubtractedJets,
                                 o2::aod::BplusChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToBplusChargedMCDetectorLevelJets,
                                 o2::aod::CandidatesBplusMCD>;

template struct JetMatchingMcSub<o2::soa::Join<o2::aod::XicToXiPiPiChargedMCDetectorLevelJets, o2::aod::XicToXiPiPiChargedMCDetectorLevelJetConstituents>,
                                 o2::soa::Join<o2::aod::XicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::XicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                 o2::aod::XicToXiPiPiChargedMCDetectorLevelJetsMatchedToXicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJets,
                                 o2::aod::XicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToXicToXiPiPiChargedMCDetectorLevelJets,
                                 o2::aod::CandidatesXicToXiPiPiMCD>;

template struct JetMatchingMcSub<o2::soa::Join<o2::aod::DielectronChargedMCDetectorLevelJets, o2::aod::DielectronChargedMCDetectorLevelJetConstituents>,
                                 o2::soa::Join<o2::aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                 o2::aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCDetectorLevelEventWiseSubtractedJets,
                                 o2::aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToDielectronChargedMCDetectorLevelJets,
                                 o2::aod::CandidatesDielectronMCD>;

template struct JetMatchingSub<o2::soa::Join<o2::aod::ChargedJets, o2::aod::ChargedJetConstituents>,
                               o2::soa::Join<o2::aod::ChargedEventWiseSubtractedJets, o2::aod::ChargedEventWiseSubtractedJetConstituents>,
                               o2::aod::ChargedJetsMatchedToChargedEventWiseSubtractedJets,
                               o2::aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets,
                               o2::aod::JTrackSubs,
                               o2::aod::JDummys>;

template struct JetMatchingSub<o2::soa::Join<o2::aod::D0ChargedJets, o2::aod::D0ChargedJetConstituents>,
                               o2::soa::Join<o2::aod::D0ChargedEventWiseSubtractedJets, o2::aod::D0ChargedEventWiseSubtractedJetConstituents>,
                               o2::aod::D0ChargedJetsMatchedToD0ChargedEventWiseSubtractedJets,
                               o2::aod::D0ChargedEventWiseSubtractedJetsMatchedToD0ChargedJets,
                               o2::aod::JTrackD0Subs,
                               o2::aod::CandidatesD0Data>;

template struct JetMatchingSub<o2::soa::Join<o2::aod::DplusChargedJets, o2::aod::DplusChargedJetConstituents>,
                               o2::soa::Join<o2::aod::DplusChargedEventWiseSubtractedJets, o2::aod::DplusChargedEventWiseSubtractedJetConstituents>,
                               o2::aod::DplusChargedJetsMatchedToDplusChargedEventWiseSubtractedJets,
                               o2::aod::DplusChargedEventWiseSubtractedJetsMatchedToDplusChargedJets,
                               o2::aod::JTrackDplusSubs,
                               o2::aod::CandidatesDplusData>;

template struct JetMatchingSub<o2::soa::Join<o2::aod::DsChargedJets, o2::aod::DsChargedJetConstituents>,
                               o2::soa::Join<o2::aod::DsChargedEventWiseSubtractedJets, o2::aod::DsChargedEventWiseSubtractedJetConstituents>,
                               o2::aod::DsChargedJetsMatchedToDsChargedEventWiseSubtractedJets,
                               o2::aod::DsChargedEventWiseSubtractedJetsMatchedToDsChargedJets,
                               o2::aod::JTrackDsSubs,
                               o2::aod::CandidatesDsData>;

template struct JetMatchingSub<o2::soa::Join<o2::aod::DstarChargedJets, o2::aod::DstarChargedJetConstituents>,
                               o2::soa::Join<o2::aod::DstarChargedEventWiseSubtractedJets, o2::aod::DstarChargedEventWiseSubtractedJetConstituents>,
                               o2::aod::DstarChargedJetsMatchedToDstarChargedEventWiseSubtractedJets,
                               o2::aod::DstarChargedEventWiseSubtractedJetsMatchedToDstarChargedJets,
                               o2::aod::JTrackDstarSubs,
                               o2::aod::CandidatesDstarData>;

template struct JetMatchingSub<o2::soa::Join<o2::aod::LcChargedJets, o2::aod::LcChargedJetConstituents>,
                               o2::soa::Join<o2::aod::LcChargedEventWiseSubtractedJets, o2::aod::LcChargedEventWiseSubtractedJetConstituents>,
                               o2::aod::LcChargedJetsMatchedToLcChargedEventWiseSubtractedJets,
                               o2::aod::LcChargedEventWiseSubtractedJetsMatchedToLcChargedJets,
                               o2::aod::JTrackLcSubs,
                               o2::aod::CandidatesLcData>;

template struct JetMatchingSub<o2::soa::Join<o2::aod::B0ChargedJets, o2::aod::B0ChargedJetConstituents>,
                               o2::soa::Join<o2::aod::B0ChargedEventWiseSubtractedJets, o2::aod::B0ChargedEventWiseSubtractedJetConstituents>,
                               o2::aod::B0ChargedJetsMatchedToB0ChargedEventWiseSubtractedJets,
                               o2::aod::B0ChargedEventWiseSubtractedJetsMatchedToB0ChargedJets,
                               o2::aod::JTrackB0Subs,
                               o2::aod::CandidatesB0Data>;

template struct JetMatchingSub<o2::soa::Join<o2::aod::BplusChargedJets, o2::aod::BplusChargedJetConstituents>,
                               o2::soa::Join<o2::aod::BplusChargedEventWiseSubtractedJets, o2::aod::BplusChargedEventWiseSubtractedJetConstituents>,
                               o2::aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets,
                               o2::aod::BplusChargedEventWiseSubtractedJetsMatchedToBplusChargedJets,
                               o2::aod::JTrackBplusSubs,
                               o2::aod::CandidatesBplusData>;

template struct JetMatchingSub<o2::soa::Join<o2::aod::XicToXiPiPiChargedJets, o2::aod::XicToXiPiPiChargedJetConstituents>,
                               o2::soa::Join<o2::aod::XicToXiPiPiChargedEventWiseSubtractedJets, o2::aod::XicToXiPiPiChargedEventWiseSubtractedJetConstituents>,
                               o2::aod::XicToXiPiPiChargedJetsMatchedToXicToXiPiPiChargedEventWiseSubtractedJets,
                               o2::aod::XicToXiPiPiChargedEventWiseSubtractedJetsMatchedToXicToXiPiPiChargedJets,
                               o2::aod::JTrackXicToXiPiPiSubs,
                               o2::aod::CandidatesXicToXiPiPiData>;

template struct JetMatchingSub<o2::soa::Join<o2::aod::DielectronChargedJets, o2::aod::DielectronChargedJetConstituents>,
                               o2::soa::Join<o2::aod::DielectronChargedEventWiseSubtractedJets, o2::aod::DielectronChargedEventWiseSubtractedJetConstituents>,
                               o2::aod::DielectronChargedJetsMatchedToDielectronChargedEventWiseSubtractedJets,
                               o2::aod::DielectronChargedEventWiseSubtractedJetsMatchedToDielectronChargedJets,
                               o2::aod::JTrackDielectronSubs,
                               o2::aod::CandidatesDielectronData>;

template struct JetMatchingDuplicates<o2::soa::Join<o2::aod::ChargedJets, o2::aod::ChargedJetConstituents>,
                                      o2::soa::Join<o2::aod::Charged1Jets, o2::aod::Charged1JetConstituents>,
                                      o2::aod::ChargedJetsMatchedToCharged1Jets,
                                      o2::aod::Charged1JetsMatchedToChargedJets,
                                      o2::aod::JTracks,
                                      o2::aod::JDummys>;

template struct JetMatchingDuplicates<o2::soa::Join<o2::aod::ChargedMCDetectorLevelJets, o2::aod::ChargedMCDetectorLevelJetConstituents>,
                                      o2::soa::Join<o2::aod::Charged1MCDetectorLevelJets, o2::aod::Charged1MCDetectorLevelJetConstituents>,
                                      o2::aod::ChargedMCDetectorLevelJetsMatchedToCharged1MCDetectorLevelJets,
                                      o2::aod::Charged1MCDetectorLevelJetsMatchedToChargedMCDetectorLevelJets,
                                      o2::aod::JTracks,
                                      o2::aod::JDummys>;

template struct JetMatchingDuplicates<o2::soa::Join<o2::aod::ChargedMCParticleLevelJets, o2::aod::ChargedMCParticleLevelJetConstituents>,
                                      o2::soa::Join<o2::aod::Charged1MCParticleLevelJets, o2::aod::Charged1MCParticleLevelJetConstituents>,
                                      o2::aod::ChargedMCParticleLevelJetsMatchedToCharged1MCParticleLevelJets,
                                      o2::aod::Charged1MCParticleLevelJetsMatchedToChargedMCParticleLevelJets,
                                      o2::aod::JMcParticles,
                                      o2::aod::JDummys>;
