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

/// \file jetmatchingmc.cxx
/// \brief matching detector level and generator level jets
///
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Aimeric Lanodu <aimeric.landou@cern.ch>
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGMC_H_
#define PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGMC_H_

#include "PWGJE/Core/JetMatchingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>

#ifndef O2_NO_WORKFLOW_MAIN
#include <Framework/runDataProcessing.h> // IWYU pragma: export
#endif

#include <vector>

template <typename JetsBase, typename JetsTag, typename JetsBasetoTagMatchingTable, typename JetsTagtoBaseMatchingTable, typename CandidatesBase, typename CandidatesTag, typename ClustersBase>
struct JetMatchingMc {

  o2::framework::Configurable<bool> doMatchingGeo{"doMatchingGeo", true, "Enable geometric matching"};
  o2::framework::Configurable<bool> doMatchingPt{"doMatchingPt", true, "Enable pt matching"};
  o2::framework::Configurable<bool> doMatchingHf{"doMatchingHf", false, "Enable HF matching"};
  o2::framework::Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.24f, "Max matching distance"};
  o2::framework::Configurable<float> minPtFraction{"minPtFraction", 0.5f, "Minimum pt fraction for pt matching"};

  o2::framework::Produces<JetsBasetoTagMatchingTable> jetsBasetoTagMatchingTable;
  o2::framework::Produces<JetsTagtoBaseMatchingTable> jetsTagtoBaseMatchingTable;

  // preslicing jet collections, only for Mc-based collection
  static constexpr bool jetsBaseIsMc = o2::soa::relatedByIndex<o2::aod::JetMcCollisions, JetsBase>();
  static constexpr bool jetsTagIsMc = o2::soa::relatedByIndex<o2::aod::JetMcCollisions, JetsTag>();

  o2::framework::Preslice<JetsBase> baseJetsPerCollision = jetsBaseIsMc ? o2::aod::jet::mcCollisionId : o2::aod::jet::collisionId;
  o2::framework::Preslice<JetsTag> tagJetsPerCollision = jetsTagIsMc ? o2::aod::jet::mcCollisionId : o2::aod::jet::collisionId;

  o2::framework::PresliceUnsorted<o2::aod::JetCollisionsMCD> CollisionsPerMcCollision = o2::aod::jmccollisionlb::mcCollisionId;

  void init(o2::framework::InitContext const&)
  {
  }

  void processJets(o2::aod::JetMcCollisions const& mcCollisions, o2::aod::JetCollisionsMCD const& collisions,
                   JetsBase const& jetsBase, JetsTag const& jetsTag,
                   o2::aod::JetTracksMCD const& tracks,
                   ClustersBase const& clusters,
                   o2::aod::JetParticles const& particles,
                   CandidatesBase const& candidatesBase,
                   CandidatesTag const& candidatesTag)
  {
    // initialise objects used to store the matching index arrays (array in case a mcCollision is split) before filling the matching tables
    std::vector<std::vector<int>> jetsBasetoTagMatchingGeo, jetsBasetoTagMatchingPt, jetsBasetoTagMatchingHF;
    std::vector<std::vector<int>> jetsTagtoBaseMatchingGeo, jetsTagtoBaseMatchingPt, jetsTagtoBaseMatchingHF;
    //  waiting for framework fix to make sliced collection of same type as original collection:
    jetsBasetoTagMatchingGeo.assign(jetsBase.size(), {});
    jetsBasetoTagMatchingPt.assign(jetsBase.size(), {});
    jetsBasetoTagMatchingHF.assign(jetsBase.size(), {});
    jetsTagtoBaseMatchingGeo.assign(jetsTag.size(), {});
    jetsTagtoBaseMatchingPt.assign(jetsTag.size(), {});
    jetsTagtoBaseMatchingHF.assign(jetsTag.size(), {});

    for (const auto& mcCollision : mcCollisions) {

      const auto collisionsPerMcColl = collisions.sliceBy(CollisionsPerMcCollision, mcCollision.globalIndex());

      for (const auto& collision : collisionsPerMcColl) {

        const auto jetsBasePerColl = jetsBase.sliceBy(baseJetsPerCollision, jetsBaseIsMc ? mcCollision.globalIndex() : collision.globalIndex());
        const auto jetsTagPerColl = jetsTag.sliceBy(tagJetsPerCollision, jetsTagIsMc ? mcCollision.globalIndex() : collision.globalIndex());

        jetmatchingutilities::doAllMatching<jetsBaseIsMc, jetsTagIsMc>(jetsBasePerColl, jetsTagPerColl, jetsBasetoTagMatchingGeo, jetsBasetoTagMatchingPt, jetsBasetoTagMatchingHF, jetsTagtoBaseMatchingGeo, jetsTagtoBaseMatchingPt, jetsTagtoBaseMatchingHF, candidatesBase, tracks, clusters, candidatesTag, particles, particles, doMatchingGeo, doMatchingHf, doMatchingPt, maxMatchingDistance, minPtFraction);
      }
    }
    for (auto i = 0; i < jetsBase.size(); ++i) {
      jetsBasetoTagMatchingTable(jetsBasetoTagMatchingGeo[i], jetsBasetoTagMatchingPt[i], jetsBasetoTagMatchingHF[i]); // is (and needs to) be filled in order
    }
    for (auto i = 0; i < jetsTag.size(); i++) {
      jetsTagtoBaseMatchingTable(jetsTagtoBaseMatchingGeo[i], jetsTagtoBaseMatchingPt[i], jetsTagtoBaseMatchingHF[i]); // is (and needs to) be filled in order
    }
  }
  PROCESS_SWITCH(JetMatchingMc, processJets, "Perform jet matching", true);
};

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::ChargedMCDetectorLevelJets, o2::aod::ChargedMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::ChargedMCParticleLevelJets, o2::aod::ChargedMCParticleLevelJetConstituents>,
                                     o2::aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets,
                                     o2::aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets,
                                     o2::aod::JCollisions,
                                     o2::aod::JMcCollisions,
                                     o2::aod::JDummys>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::NeutralMCDetectorLevelJets, o2::aod::NeutralMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::NeutralMCParticleLevelJets, o2::aod::NeutralMCParticleLevelJetConstituents>,
                                     o2::aod::NeutralMCDetectorLevelJetsMatchedToNeutralMCParticleLevelJets,
                                     o2::aod::NeutralMCParticleLevelJetsMatchedToNeutralMCDetectorLevelJets,
                                     o2::aod::JCollisions,
                                     o2::aod::JMcCollisions,
                                     o2::aod::JetClustersMCD>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::FullMCDetectorLevelJets, o2::aod::FullMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::FullMCParticleLevelJets, o2::aod::FullMCParticleLevelJetConstituents>,
                                     o2::aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets,
                                     o2::aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets,
                                     o2::aod::JCollisions,
                                     o2::aod::JMcCollisions,
                                     o2::aod::JetClustersMCD>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::D0ChargedMCDetectorLevelJets, o2::aod::D0ChargedMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::D0ChargedMCParticleLevelJets, o2::aod::D0ChargedMCParticleLevelJetConstituents>,
                                     o2::aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets,
                                     o2::aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets,
                                     o2::aod::CandidatesD0MCD,
                                     o2::aod::CandidatesD0MCP,
                                     o2::aod::JDummys>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::DplusChargedMCDetectorLevelJets, o2::aod::DplusChargedMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::DplusChargedMCParticleLevelJets, o2::aod::DplusChargedMCParticleLevelJetConstituents>,
                                     o2::aod::DplusChargedMCDetectorLevelJetsMatchedToDplusChargedMCParticleLevelJets,
                                     o2::aod::DplusChargedMCParticleLevelJetsMatchedToDplusChargedMCDetectorLevelJets,
                                     o2::aod::CandidatesDplusMCD,
                                     o2::aod::CandidatesDplusMCP,
                                     o2::aod::JDummys>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::DsChargedMCDetectorLevelJets, o2::aod::DsChargedMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::DsChargedMCParticleLevelJets, o2::aod::DsChargedMCParticleLevelJetConstituents>,
                                     o2::aod::DsChargedMCDetectorLevelJetsMatchedToDsChargedMCParticleLevelJets,
                                     o2::aod::DsChargedMCParticleLevelJetsMatchedToDsChargedMCDetectorLevelJets,
                                     o2::aod::CandidatesDsMCD,
                                     o2::aod::CandidatesDsMCP,
                                     o2::aod::JDummys>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::DstarChargedMCDetectorLevelJets, o2::aod::DstarChargedMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::DstarChargedMCParticleLevelJets, o2::aod::DstarChargedMCParticleLevelJetConstituents>,
                                     o2::aod::DstarChargedMCDetectorLevelJetsMatchedToDstarChargedMCParticleLevelJets,
                                     o2::aod::DstarChargedMCParticleLevelJetsMatchedToDstarChargedMCDetectorLevelJets,
                                     o2::aod::CandidatesDstarMCD,
                                     o2::aod::CandidatesDstarMCP,
                                     o2::aod::JDummys>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::LcChargedMCDetectorLevelJets, o2::aod::LcChargedMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::LcChargedMCParticleLevelJets, o2::aod::LcChargedMCParticleLevelJetConstituents>,
                                     o2::aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets,
                                     o2::aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets,
                                     o2::aod::CandidatesLcMCD,
                                     o2::aod::CandidatesLcMCP,
                                     o2::aod::JDummys>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::B0ChargedMCDetectorLevelJets, o2::aod::B0ChargedMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::B0ChargedMCParticleLevelJets, o2::aod::B0ChargedMCParticleLevelJetConstituents>,
                                     o2::aod::B0ChargedMCDetectorLevelJetsMatchedToB0ChargedMCParticleLevelJets,
                                     o2::aod::B0ChargedMCParticleLevelJetsMatchedToB0ChargedMCDetectorLevelJets,
                                     o2::aod::CandidatesB0MCD,
                                     o2::aod::CandidatesB0MCP,
                                     o2::aod::JDummys>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::BplusChargedMCDetectorLevelJets, o2::aod::BplusChargedMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::BplusChargedMCParticleLevelJets, o2::aod::BplusChargedMCParticleLevelJetConstituents>,
                                     o2::aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets,
                                     o2::aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets,
                                     o2::aod::CandidatesBplusMCD,
                                     o2::aod::CandidatesBplusMCP,
                                     o2::aod::JDummys>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::XicToXiPiPiChargedMCDetectorLevelJets, o2::aod::XicToXiPiPiChargedMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::XicToXiPiPiChargedMCParticleLevelJets, o2::aod::XicToXiPiPiChargedMCParticleLevelJetConstituents>,
                                     o2::aod::XicToXiPiPiChargedMCDetectorLevelJetsMatchedToXicToXiPiPiChargedMCParticleLevelJets,
                                     o2::aod::XicToXiPiPiChargedMCParticleLevelJetsMatchedToXicToXiPiPiChargedMCDetectorLevelJets,
                                     o2::aod::CandidatesXicToXiPiPiMCD,
                                     o2::aod::CandidatesXicToXiPiPiMCP,
                                     o2::aod::JDummys>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::DielectronChargedMCDetectorLevelJets, o2::aod::DielectronChargedMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::DielectronChargedMCParticleLevelJets, o2::aod::DielectronChargedMCParticleLevelJetConstituents>,
                                     o2::aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCParticleLevelJets,
                                     o2::aod::DielectronChargedMCParticleLevelJetsMatchedToDielectronChargedMCDetectorLevelJets,
                                     o2::aod::CandidatesDielectronMCD,
                                     o2::aod::CandidatesDielectronMCP,
                                     o2::aod::JDummys>;

extern template struct JetMatchingMc<o2::soa::Join<o2::aod::V0ChargedMCDetectorLevelJets, o2::aod::V0ChargedMCDetectorLevelJetConstituents>,
                                     o2::soa::Join<o2::aod::V0ChargedMCParticleLevelJets, o2::aod::V0ChargedMCParticleLevelJetConstituents>,
                                     o2::aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets,
                                     o2::aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets,
                                     o2::aod::CandidatesV0MCD,
                                     o2::aod::CandidatesV0MCP,
                                     o2::aod::JDummys>;

#endif // PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGMC_H_
