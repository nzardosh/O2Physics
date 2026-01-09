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

/// \file jetSubstructureMatchingTemplateInstantiations.cxx
/// \brief Explicit Template Instantiation file for all JetSubstructureMatching templates
///
/// This file contains explicit template instantiations for all JetSubstructureMatching template classes
/// to improve compilation performance by compiling templates once instead of in every
/// translation unit that uses them.
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Florian Jonas <florian.jonas@cern.ch>

// Prevent the main function from being defined when including headers
#define O2_NO_WORKFLOW_MAIN

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/TableProducer/Matching/Substructure/jetSubstructureMatchingMC.h"
#include "PWGJE/TableProducer/Matching/Substructure/jetSubstructureMatchingSub.h"

// ============================================================================
// Explicit Template Instantiations for JetSubstructureMatching
// ============================================================================

template struct JetSubstructureMatchingMC<o2::soa::Join<o2::aod::ChargedMCDetectorLevelJets, o2::aod::ChargedMCDetectorLevelJetConstituents, o2::aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>,
                                          o2::soa::Join<o2::aod::ChargedMCParticleLevelJets, o2::aod::ChargedMCParticleLevelJetConstituents, o2::aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>,
                                          o2::aod::ChargedMCDetectorLevelSPsMatchedToChargedMCParticleLevelSPs,
                                          o2::aod::ChargedMCParticleLevelSPsMatchedToChargedMCDetectorLevelSPs,
                                          o2::aod::ChargedMCDetectorLevelPRsMatchedToChargedMCParticleLevelPRs,
                                          o2::aod::ChargedMCParticleLevelPRsMatchedToChargedMCDetectorLevelPRs,
                                          o2::aod::ChargedMCDetectorLevelSPs,
                                          o2::aod::ChargedMCParticleLevelSPs,
                                          o2::aod::ChargedMCDetectorLevelPRs,
                                          o2::aod::ChargedMCParticleLevelPRs,
                                          o2::aod::JCollisions,
                                          o2::aod::JMcCollisions,
                                          o2::aod::JetTracksMCD,
                                          o2::aod::JetParticles,
                                          o2::aod::JDummys>;

template struct JetSubstructureMatchingMC<o2::soa::Join<o2::aod::D0ChargedMCDetectorLevelJets, o2::aod::D0ChargedMCDetectorLevelJetConstituents, o2::aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>,
                                          o2::soa::Join<o2::aod::D0ChargedMCParticleLevelJets, o2::aod::D0ChargedMCParticleLevelJetConstituents, o2::aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>,
                                          o2::aod::D0ChargedMCDetectorLevelSPsMatchedToD0ChargedMCParticleLevelSPs,
                                          o2::aod::D0ChargedMCParticleLevelSPsMatchedToD0ChargedMCDetectorLevelSPs,
                                          o2::aod::D0ChargedMCDetectorLevelPRsMatchedToD0ChargedMCParticleLevelPRs,
                                          o2::aod::D0ChargedMCParticleLevelPRsMatchedToD0ChargedMCDetectorLevelPRs,
                                          o2::aod::D0ChargedMCDetectorLevelSPs,
                                          o2::aod::D0ChargedMCParticleLevelSPs,
                                          o2::aod::D0ChargedMCDetectorLevelPRs,
                                          o2::aod::D0ChargedMCParticleLevelPRs,
                                          o2::aod::CandidatesD0MCD,
                                          o2::aod::CandidatesD0MCP,
                                          o2::aod::JetTracksMCD,
                                          o2::aod::JetParticles,
                                          o2::aod::JDummys>;

template struct JetSubstructureMatchingMC<o2::soa::Join<o2::aod::DplusChargedMCDetectorLevelJets, o2::aod::DplusChargedMCDetectorLevelJetConstituents, o2::aod::DplusChargedMCDetectorLevelJetsMatchedToDplusChargedMCParticleLevelJets>,
                                          o2::soa::Join<o2::aod::DplusChargedMCParticleLevelJets, o2::aod::DplusChargedMCParticleLevelJetConstituents, o2::aod::DplusChargedMCParticleLevelJetsMatchedToDplusChargedMCDetectorLevelJets>,
                                          o2::aod::DplusChargedMCDetectorLevelSPsMatchedToDplusChargedMCParticleLevelSPs,
                                          o2::aod::DplusChargedMCParticleLevelSPsMatchedToDplusChargedMCDetectorLevelSPs,
                                          o2::aod::DplusChargedMCDetectorLevelPRsMatchedToDplusChargedMCParticleLevelPRs,
                                          o2::aod::DplusChargedMCParticleLevelPRsMatchedToDplusChargedMCDetectorLevelPRs,
                                          o2::aod::DplusChargedMCDetectorLevelSPs,
                                          o2::aod::DplusChargedMCParticleLevelSPs,
                                          o2::aod::DplusChargedMCDetectorLevelPRs,
                                          o2::aod::DplusChargedMCParticleLevelPRs,
                                          o2::aod::CandidatesDplusMCD,
                                          o2::aod::CandidatesDplusMCP,
                                          o2::aod::JetTracksMCD,
                                          o2::aod::JetParticles,
                                          o2::aod::JDummys>;

template struct JetSubstructureMatchingMC<o2::soa::Join<o2::aod::DsChargedMCDetectorLevelJets, o2::aod::DsChargedMCDetectorLevelJetConstituents, o2::aod::DsChargedMCDetectorLevelJetsMatchedToDsChargedMCParticleLevelJets>,
                                          o2::soa::Join<o2::aod::DsChargedMCParticleLevelJets, o2::aod::DsChargedMCParticleLevelJetConstituents, o2::aod::DsChargedMCParticleLevelJetsMatchedToDsChargedMCDetectorLevelJets>,
                                          o2::aod::DsChargedMCDetectorLevelSPsMatchedToDsChargedMCParticleLevelSPs,
                                          o2::aod::DsChargedMCParticleLevelSPsMatchedToDsChargedMCDetectorLevelSPs,
                                          o2::aod::DsChargedMCDetectorLevelPRsMatchedToDsChargedMCParticleLevelPRs,
                                          o2::aod::DsChargedMCParticleLevelPRsMatchedToDsChargedMCDetectorLevelPRs,
                                          o2::aod::DsChargedMCDetectorLevelSPs,
                                          o2::aod::DsChargedMCParticleLevelSPs,
                                          o2::aod::DsChargedMCDetectorLevelPRs,
                                          o2::aod::DsChargedMCParticleLevelPRs,
                                          o2::aod::CandidatesDsMCD,
                                          o2::aod::CandidatesDsMCP,
                                          o2::aod::JetTracksMCD,
                                          o2::aod::JetParticles,
                                          o2::aod::JDummys>;

template struct JetSubstructureMatchingMC<o2::soa::Join<o2::aod::DstarChargedMCDetectorLevelJets, o2::aod::DstarChargedMCDetectorLevelJetConstituents, o2::aod::DstarChargedMCDetectorLevelJetsMatchedToDstarChargedMCParticleLevelJets>,
                                          o2::soa::Join<o2::aod::DstarChargedMCParticleLevelJets, o2::aod::DstarChargedMCParticleLevelJetConstituents, o2::aod::DstarChargedMCParticleLevelJetsMatchedToDstarChargedMCDetectorLevelJets>,
                                          o2::aod::DstarChargedMCDetectorLevelSPsMatchedToDstarChargedMCParticleLevelSPs,
                                          o2::aod::DstarChargedMCParticleLevelSPsMatchedToDstarChargedMCDetectorLevelSPs,
                                          o2::aod::DstarChargedMCDetectorLevelPRsMatchedToDstarChargedMCParticleLevelPRs,
                                          o2::aod::DstarChargedMCParticleLevelPRsMatchedToDstarChargedMCDetectorLevelPRs,
                                          o2::aod::DstarChargedMCDetectorLevelSPs,
                                          o2::aod::DstarChargedMCParticleLevelSPs,
                                          o2::aod::DstarChargedMCDetectorLevelPRs,
                                          o2::aod::DstarChargedMCParticleLevelPRs,
                                          o2::aod::CandidatesDstarMCD,
                                          o2::aod::CandidatesDstarMCP,
                                          o2::aod::JetTracksMCD,
                                          o2::aod::JetParticles,
                                          o2::aod::JDummys>;

template struct JetSubstructureMatchingMC<o2::soa::Join<o2::aod::LcChargedMCDetectorLevelJets, o2::aod::LcChargedMCDetectorLevelJetConstituents, o2::aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets>,
                                          o2::soa::Join<o2::aod::LcChargedMCParticleLevelJets, o2::aod::LcChargedMCParticleLevelJetConstituents, o2::aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets>,
                                          o2::aod::LcChargedMCDetectorLevelSPsMatchedToLcChargedMCParticleLevelSPs,
                                          o2::aod::LcChargedMCParticleLevelSPsMatchedToLcChargedMCDetectorLevelSPs,
                                          o2::aod::LcChargedMCDetectorLevelPRsMatchedToLcChargedMCParticleLevelPRs,
                                          o2::aod::LcChargedMCParticleLevelPRsMatchedToLcChargedMCDetectorLevelPRs,
                                          o2::aod::LcChargedMCDetectorLevelSPs,
                                          o2::aod::LcChargedMCParticleLevelSPs,
                                          o2::aod::LcChargedMCDetectorLevelPRs,
                                          o2::aod::LcChargedMCParticleLevelPRs,
                                          o2::aod::CandidatesLcMCD,
                                          o2::aod::CandidatesLcMCP,
                                          o2::aod::JetTracksMCD,
                                          o2::aod::JetParticles,
                                          o2::aod::JDummys>;

template struct JetSubstructureMatchingMC<o2::soa::Join<o2::aod::B0ChargedMCDetectorLevelJets, o2::aod::B0ChargedMCDetectorLevelJetConstituents, o2::aod::B0ChargedMCDetectorLevelJetsMatchedToB0ChargedMCParticleLevelJets>,
                                          o2::soa::Join<o2::aod::B0ChargedMCParticleLevelJets, o2::aod::B0ChargedMCParticleLevelJetConstituents, o2::aod::B0ChargedMCParticleLevelJetsMatchedToB0ChargedMCDetectorLevelJets>,
                                          o2::aod::B0ChargedMCDetectorLevelSPsMatchedToB0ChargedMCParticleLevelSPs,
                                          o2::aod::B0ChargedMCParticleLevelSPsMatchedToB0ChargedMCDetectorLevelSPs,
                                          o2::aod::B0ChargedMCDetectorLevelPRsMatchedToB0ChargedMCParticleLevelPRs,
                                          o2::aod::B0ChargedMCParticleLevelPRsMatchedToB0ChargedMCDetectorLevelPRs,
                                          o2::aod::B0ChargedMCDetectorLevelSPs,
                                          o2::aod::B0ChargedMCParticleLevelSPs,
                                          o2::aod::B0ChargedMCDetectorLevelPRs,
                                          o2::aod::B0ChargedMCParticleLevelPRs,
                                          o2::aod::CandidatesB0MCD,
                                          o2::aod::CandidatesB0MCP,
                                          o2::aod::JetTracksMCD,
                                          o2::aod::JetParticles,
                                          o2::aod::JDummys>;

template struct JetSubstructureMatchingMC<o2::soa::Join<o2::aod::BplusChargedMCDetectorLevelJets, o2::aod::BplusChargedMCDetectorLevelJetConstituents, o2::aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets>,
                                          o2::soa::Join<o2::aod::BplusChargedMCParticleLevelJets, o2::aod::BplusChargedMCParticleLevelJetConstituents, o2::aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets>,
                                          o2::aod::BplusChargedMCDetectorLevelSPsMatchedToBplusChargedMCParticleLevelSPs,
                                          o2::aod::BplusChargedMCParticleLevelSPsMatchedToBplusChargedMCDetectorLevelSPs,
                                          o2::aod::BplusChargedMCDetectorLevelPRsMatchedToBplusChargedMCParticleLevelPRs,
                                          o2::aod::BplusChargedMCParticleLevelPRsMatchedToBplusChargedMCDetectorLevelPRs,
                                          o2::aod::BplusChargedMCDetectorLevelSPs,
                                          o2::aod::BplusChargedMCParticleLevelSPs,
                                          o2::aod::BplusChargedMCDetectorLevelPRs,
                                          o2::aod::BplusChargedMCParticleLevelPRs,
                                          o2::aod::CandidatesBplusMCD,
                                          o2::aod::CandidatesBplusMCP,
                                          o2::aod::JetTracksMCD,
                                          o2::aod::JetParticles,
                                          o2::aod::JDummys>;

template struct JetSubstructureMatchingMC<o2::soa::Join<o2::aod::XicToXiPiPiChargedMCDetectorLevelJets, o2::aod::XicToXiPiPiChargedMCDetectorLevelJetConstituents, o2::aod::XicToXiPiPiChargedMCDetectorLevelJetsMatchedToXicToXiPiPiChargedMCParticleLevelJets>,
                                          o2::soa::Join<o2::aod::XicToXiPiPiChargedMCParticleLevelJets, o2::aod::XicToXiPiPiChargedMCParticleLevelJetConstituents, o2::aod::XicToXiPiPiChargedMCParticleLevelJetsMatchedToXicToXiPiPiChargedMCDetectorLevelJets>,
                                          o2::aod::XicToXiPiPiChargedMCDetectorLevelSPsMatchedToXicToXiPiPiChargedMCParticleLevelSPs,
                                          o2::aod::XicToXiPiPiChargedMCParticleLevelSPsMatchedToXicToXiPiPiChargedMCDetectorLevelSPs,
                                          o2::aod::XicToXiPiPiChargedMCDetectorLevelPRsMatchedToXicToXiPiPiChargedMCParticleLevelPRs,
                                          o2::aod::XicToXiPiPiChargedMCParticleLevelPRsMatchedToXicToXiPiPiChargedMCDetectorLevelPRs,
                                          o2::aod::XicToXiPiPiChargedMCDetectorLevelSPs,
                                          o2::aod::XicToXiPiPiChargedMCParticleLevelSPs,
                                          o2::aod::XicToXiPiPiChargedMCDetectorLevelPRs,
                                          o2::aod::XicToXiPiPiChargedMCParticleLevelPRs,
                                          o2::aod::CandidatesXicToXiPiPiMCD,
                                          o2::aod::CandidatesXicToXiPiPiMCP,
                                          o2::aod::JetTracksMCD,
                                          o2::aod::JetParticles,
                                          o2::aod::JDummys>;

template struct JetSubstructureMatchingMC<o2::soa::Join<o2::aod::DielectronChargedMCDetectorLevelJets, o2::aod::DielectronChargedMCDetectorLevelJetConstituents, o2::aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCParticleLevelJets>,
                                          o2::soa::Join<o2::aod::DielectronChargedMCParticleLevelJets, o2::aod::DielectronChargedMCParticleLevelJetConstituents, o2::aod::DielectronChargedMCParticleLevelJetsMatchedToDielectronChargedMCDetectorLevelJets>,
                                          o2::aod::DielectronChargedMCDetectorLevelSPsMatchedToDielectronChargedMCParticleLevelSPs,
                                          o2::aod::DielectronChargedMCParticleLevelSPsMatchedToDielectronChargedMCDetectorLevelSPs,
                                          o2::aod::DielectronChargedMCDetectorLevelPRsMatchedToDielectronChargedMCParticleLevelPRs,
                                          o2::aod::DielectronChargedMCParticleLevelPRsMatchedToDielectronChargedMCDetectorLevelPRs,
                                          o2::aod::DielectronChargedMCDetectorLevelSPs,
                                          o2::aod::DielectronChargedMCParticleLevelSPs,
                                          o2::aod::DielectronChargedMCDetectorLevelPRs,
                                          o2::aod::DielectronChargedMCParticleLevelPRs,
                                          o2::aod::CandidatesDielectronMCD,
                                          o2::aod::CandidatesDielectronMCP,
                                          o2::aod::JetTracksMCD,
                                          o2::aod::JetParticles,
                                          o2::aod::JDummys>;

template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::ChargedJets, o2::aod::ChargedJetConstituents, o2::aod::ChargedJetsMatchedToChargedEventWiseSubtractedJets>,
                                           o2::soa::Join<o2::aod::ChargedEventWiseSubtractedJets, o2::aod::ChargedEventWiseSubtractedJetConstituents, o2::aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets>,
                                           o2::aod::ChargedSPsMatchedToChargedEventWiseSubtractedSPs,
                                           o2::aod::ChargedEventWiseSubtractedSPsMatchedToChargedSPs,
                                           o2::aod::ChargedPRsMatchedToChargedEventWiseSubtractedPRs,
                                           o2::aod::ChargedEventWiseSubtractedPRsMatchedToChargedPRs,
                                           o2::aod::ChargedSPs,
                                           o2::aod::ChargedEventWiseSubtractedSPs,
                                           o2::aod::ChargedPRs,
                                           o2::aod::ChargedEventWiseSubtractedPRs,
                                           o2::aod::JCollisions,
                                           o2::aod::JetTracks,
                                           o2::aod::JetTracksSub,
                                           o2::aod::JDummys>;

template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::D0ChargedJets, o2::aod::D0ChargedJetConstituents, o2::aod::D0ChargedJetsMatchedToD0ChargedEventWiseSubtractedJets>,
                                           o2::soa::Join<o2::aod::D0ChargedEventWiseSubtractedJets, o2::aod::D0ChargedEventWiseSubtractedJetConstituents, o2::aod::D0ChargedEventWiseSubtractedJetsMatchedToD0ChargedJets>,
                                           o2::aod::D0ChargedSPsMatchedToD0ChargedEventWiseSubtractedSPs,
                                           o2::aod::D0ChargedEventWiseSubtractedSPsMatchedToD0ChargedSPs,
                                           o2::aod::D0ChargedPRsMatchedToD0ChargedEventWiseSubtractedPRs,
                                           o2::aod::D0ChargedEventWiseSubtractedPRsMatchedToD0ChargedPRs,
                                           o2::aod::D0ChargedSPs,
                                           o2::aod::D0ChargedEventWiseSubtractedSPs,
                                           o2::aod::D0ChargedPRs,
                                           o2::aod::D0ChargedEventWiseSubtractedPRs,
                                           o2::aod::CandidatesD0Data,
                                           o2::aod::JetTracks,
                                           o2::aod::JetTracksSubD0,
                                           o2::aod::JDummys>;

template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::DplusChargedJets, o2::aod::DplusChargedJetConstituents, o2::aod::DplusChargedJetsMatchedToDplusChargedEventWiseSubtractedJets>,
                                           o2::soa::Join<o2::aod::DplusChargedEventWiseSubtractedJets, o2::aod::DplusChargedEventWiseSubtractedJetConstituents, o2::aod::DplusChargedEventWiseSubtractedJetsMatchedToDplusChargedJets>,
                                           o2::aod::DplusChargedSPsMatchedToDplusChargedEventWiseSubtractedSPs,
                                           o2::aod::DplusChargedEventWiseSubtractedSPsMatchedToDplusChargedSPs,
                                           o2::aod::DplusChargedPRsMatchedToDplusChargedEventWiseSubtractedPRs,
                                           o2::aod::DplusChargedEventWiseSubtractedPRsMatchedToDplusChargedPRs,
                                           o2::aod::DplusChargedSPs,
                                           o2::aod::DplusChargedEventWiseSubtractedSPs,
                                           o2::aod::DplusChargedPRs,
                                           o2::aod::DplusChargedEventWiseSubtractedPRs,
                                           o2::aod::CandidatesDplusData,
                                           o2::aod::JetTracks,
                                           o2::aod::JetTracksSubDplus,
                                           o2::aod::JDummys>;

template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::DsChargedJets, o2::aod::DsChargedJetConstituents, o2::aod::DsChargedJetsMatchedToDsChargedEventWiseSubtractedJets>,
                                           o2::soa::Join<o2::aod::DsChargedEventWiseSubtractedJets, o2::aod::DsChargedEventWiseSubtractedJetConstituents, o2::aod::DsChargedEventWiseSubtractedJetsMatchedToDsChargedJets>,
                                           o2::aod::DsChargedSPsMatchedToDsChargedEventWiseSubtractedSPs,
                                           o2::aod::DsChargedEventWiseSubtractedSPsMatchedToDsChargedSPs,
                                           o2::aod::DsChargedPRsMatchedToDsChargedEventWiseSubtractedPRs,
                                           o2::aod::DsChargedEventWiseSubtractedPRsMatchedToDsChargedPRs,
                                           o2::aod::DsChargedSPs,
                                           o2::aod::DsChargedEventWiseSubtractedSPs,
                                           o2::aod::DsChargedPRs,
                                           o2::aod::DsChargedEventWiseSubtractedPRs,
                                           o2::aod::CandidatesDsData,
                                           o2::aod::JetTracks,
                                           o2::aod::JetTracksSubDs,
                                           o2::aod::JDummys>;

template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::DstarChargedJets, o2::aod::DstarChargedJetConstituents, o2::aod::DstarChargedJetsMatchedToDstarChargedEventWiseSubtractedJets>,
                                           o2::soa::Join<o2::aod::DstarChargedEventWiseSubtractedJets, o2::aod::DstarChargedEventWiseSubtractedJetConstituents, o2::aod::DstarChargedEventWiseSubtractedJetsMatchedToDstarChargedJets>,
                                           o2::aod::DstarChargedSPsMatchedToDstarChargedEventWiseSubtractedSPs,
                                           o2::aod::DstarChargedEventWiseSubtractedSPsMatchedToDstarChargedSPs,
                                           o2::aod::DstarChargedPRsMatchedToDstarChargedEventWiseSubtractedPRs,
                                           o2::aod::DstarChargedEventWiseSubtractedPRsMatchedToDstarChargedPRs,
                                           o2::aod::DstarChargedSPs,
                                           o2::aod::DstarChargedEventWiseSubtractedSPs,
                                           o2::aod::DstarChargedPRs,
                                           o2::aod::DstarChargedEventWiseSubtractedPRs,
                                           o2::aod::CandidatesDstarData,
                                           o2::aod::JetTracks,
                                           o2::aod::JetTracksSubDstar,
                                           o2::aod::JDummys>;

template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::LcChargedJets, o2::aod::LcChargedJetConstituents, o2::aod::LcChargedJetsMatchedToLcChargedEventWiseSubtractedJets>,
                                           o2::soa::Join<o2::aod::LcChargedEventWiseSubtractedJets, o2::aod::LcChargedEventWiseSubtractedJetConstituents, o2::aod::LcChargedEventWiseSubtractedJetsMatchedToLcChargedJets>,
                                           o2::aod::LcChargedSPsMatchedToLcChargedEventWiseSubtractedSPs,
                                           o2::aod::LcChargedEventWiseSubtractedSPsMatchedToLcChargedSPs,
                                           o2::aod::LcChargedPRsMatchedToLcChargedEventWiseSubtractedPRs,
                                           o2::aod::LcChargedEventWiseSubtractedPRsMatchedToLcChargedPRs,
                                           o2::aod::LcChargedSPs,
                                           o2::aod::LcChargedEventWiseSubtractedSPs,
                                           o2::aod::LcChargedPRs,
                                           o2::aod::LcChargedEventWiseSubtractedPRs,
                                           o2::aod::CandidatesLcData,
                                           o2::aod::JetTracks,
                                           o2::aod::JetTracksSubLc,
                                           o2::aod::JDummys>;

template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::B0ChargedJets, o2::aod::B0ChargedJetConstituents, o2::aod::B0ChargedJetsMatchedToB0ChargedEventWiseSubtractedJets>,
                                           o2::soa::Join<o2::aod::B0ChargedEventWiseSubtractedJets, o2::aod::B0ChargedEventWiseSubtractedJetConstituents, o2::aod::B0ChargedEventWiseSubtractedJetsMatchedToB0ChargedJets>,
                                           o2::aod::B0ChargedSPsMatchedToB0ChargedEventWiseSubtractedSPs,
                                           o2::aod::B0ChargedEventWiseSubtractedSPsMatchedToB0ChargedSPs,
                                           o2::aod::B0ChargedPRsMatchedToB0ChargedEventWiseSubtractedPRs,
                                           o2::aod::B0ChargedEventWiseSubtractedPRsMatchedToB0ChargedPRs,
                                           o2::aod::B0ChargedSPs,
                                           o2::aod::B0ChargedEventWiseSubtractedSPs,
                                           o2::aod::B0ChargedPRs,
                                           o2::aod::B0ChargedEventWiseSubtractedPRs,
                                           o2::aod::CandidatesB0Data,
                                           o2::aod::JetTracks,
                                           o2::aod::JetTracksSubB0,
                                           o2::aod::JDummys>;

template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::BplusChargedJets, o2::aod::BplusChargedJetConstituents, o2::aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets>,
                                           o2::soa::Join<o2::aod::BplusChargedEventWiseSubtractedJets, o2::aod::BplusChargedEventWiseSubtractedJetConstituents, o2::aod::BplusChargedEventWiseSubtractedJetsMatchedToBplusChargedJets>,
                                           o2::aod::BplusChargedSPsMatchedToBplusChargedEventWiseSubtractedSPs,
                                           o2::aod::BplusChargedEventWiseSubtractedSPsMatchedToBplusChargedSPs,
                                           o2::aod::BplusChargedPRsMatchedToBplusChargedEventWiseSubtractedPRs,
                                           o2::aod::BplusChargedEventWiseSubtractedPRsMatchedToBplusChargedPRs,
                                           o2::aod::BplusChargedSPs,
                                           o2::aod::BplusChargedEventWiseSubtractedSPs,
                                           o2::aod::BplusChargedPRs,
                                           o2::aod::BplusChargedEventWiseSubtractedPRs,
                                           o2::aod::CandidatesBplusData,
                                           o2::aod::JetTracks,
                                           o2::aod::JetTracksSubBplus,
                                           o2::aod::JDummys>;

template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::XicToXiPiPiChargedJets, o2::aod::XicToXiPiPiChargedJetConstituents, o2::aod::XicToXiPiPiChargedJetsMatchedToXicToXiPiPiChargedEventWiseSubtractedJets>,
                                           o2::soa::Join<o2::aod::XicToXiPiPiChargedEventWiseSubtractedJets, o2::aod::XicToXiPiPiChargedEventWiseSubtractedJetConstituents, o2::aod::XicToXiPiPiChargedEventWiseSubtractedJetsMatchedToXicToXiPiPiChargedJets>,
                                           o2::aod::XicToXiPiPiChargedSPsMatchedToXicToXiPiPiChargedEventWiseSubtractedSPs,
                                           o2::aod::XicToXiPiPiChargedEventWiseSubtractedSPsMatchedToXicToXiPiPiChargedSPs,
                                           o2::aod::XicToXiPiPiChargedPRsMatchedToXicToXiPiPiChargedEventWiseSubtractedPRs,
                                           o2::aod::XicToXiPiPiChargedEventWiseSubtractedPRsMatchedToXicToXiPiPiChargedPRs,
                                           o2::aod::XicToXiPiPiChargedSPs,
                                           o2::aod::XicToXiPiPiChargedEventWiseSubtractedSPs,
                                           o2::aod::XicToXiPiPiChargedPRs,
                                           o2::aod::XicToXiPiPiChargedEventWiseSubtractedPRs,
                                           o2::aod::CandidatesXicToXiPiPiData,
                                           o2::aod::JetTracks,
                                           o2::aod::JetTracksSubXicToXiPiPi,
                                           o2::aod::JDummys>;

template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::DielectronChargedJets, o2::aod::DielectronChargedJetConstituents, o2::aod::DielectronChargedJetsMatchedToDielectronChargedEventWiseSubtractedJets>,
                                           o2::soa::Join<o2::aod::DielectronChargedEventWiseSubtractedJets, o2::aod::DielectronChargedEventWiseSubtractedJetConstituents, o2::aod::DielectronChargedEventWiseSubtractedJetsMatchedToDielectronChargedJets>,
                                           o2::aod::DielectronChargedSPsMatchedToDielectronChargedEventWiseSubtractedSPs,
                                           o2::aod::DielectronChargedEventWiseSubtractedSPsMatchedToDielectronChargedSPs,
                                           o2::aod::DielectronChargedPRsMatchedToDielectronChargedEventWiseSubtractedPRs,
                                           o2::aod::DielectronChargedEventWiseSubtractedPRsMatchedToDielectronChargedPRs,
                                           o2::aod::DielectronChargedSPs,
                                           o2::aod::DielectronChargedEventWiseSubtractedSPs,
                                           o2::aod::DielectronChargedPRs,
                                           o2::aod::DielectronChargedEventWiseSubtractedPRs,
                                           o2::aod::CandidatesDielectronData,
                                           o2::aod::JetTracks,
                                           o2::aod::JetTracksSubDielectron,
                                           o2::aod::JDummys>;
