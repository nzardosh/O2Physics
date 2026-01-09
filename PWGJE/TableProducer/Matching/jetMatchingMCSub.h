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

/// \file jetmatchingmcsub.cxx
/// \brief matching event-wise constituent subtracted detector level and unsubtracted generated level jets (this is usseful as a template for embedding  matching)
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGMCSUB_H_
#define PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGMCSUB_H_

#include "PWGJE/Core/JetMatchingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>

#ifndef O2_NO_WORKFLOW_MAIN
#include <Framework/runDataProcessing.h> // IWYU pragma: export
#endif

#include <vector>

template <typename JetsBase, typename JetsTag, typename JetsBasetoTagMatchingTable, typename JetsTagtoBaseMatchingTable, typename Candidates>
struct JetMatchingMcSub {

  o2::framework::Configurable<bool> doMatchingGeo{"doMatchingGeo", true, "Enable geometric matching"};
  o2::framework::Configurable<bool> doMatchingPt{"doMatchingPt", true, "Enable pt matching"};
  o2::framework::Configurable<bool> doMatchingHf{"doMatchingHf", false, "Enable HF matching"};
  o2::framework::Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.24f, "Max matching distance"};
  o2::framework::Configurable<float> minPtFraction{"minPtFraction", 0.5f, "Minimum pt fraction for pt matching"};

  o2::framework::Produces<JetsBasetoTagMatchingTable> jetsBasetoTagMatchingTable;
  o2::framework::Produces<JetsTagtoBaseMatchingTable> jetsTagtoBaseMatchingTable;

  // preslicing jet collections, only for Mc-based collection
  static constexpr bool jetsBaseIsMc = false;
  static constexpr bool jetsTagIsMc = false;

  o2::framework::Preslice<JetsBase> baseJetsPerCollision = o2::aod::jet::collisionId;
  o2::framework::Preslice<JetsTag> tagJetsPerCollision = o2::aod::jet::collisionId;

  void init(o2::framework::InitContext const&)
  {
  }

  void processJets(o2::aod::JetCollisions const& collisions,
                   JetsBase const& jetsBase, JetsTag const& jetsTag,
                   o2::aod::JetTracks const& tracks,
                   o2::aod::JetTracksSub const& tracksSub,
                   Candidates const& candidates)
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

    for (const auto& collision : collisions) {

      const auto jetsBasePerColl = jetsBase.sliceBy(baseJetsPerCollision, collision.globalIndex());
      const auto jetsTagPerColl = jetsTag.sliceBy(tagJetsPerCollision, collision.globalIndex());

      jetmatchingutilities::doAllMatching<jetsBaseIsMc, jetsTagIsMc>(jetsBasePerColl, jetsTagPerColl, jetsBasetoTagMatchingGeo, jetsBasetoTagMatchingPt, jetsBasetoTagMatchingHF, jetsTagtoBaseMatchingGeo, jetsTagtoBaseMatchingPt, jetsTagtoBaseMatchingHF, candidates, tracks, tracks, candidates, tracksSub, tracksSub, doMatchingGeo, doMatchingHf, doMatchingPt, maxMatchingDistance, minPtFraction);
    }

    for (auto i = 0; i < jetsBase.size(); ++i) {
      jetsBasetoTagMatchingTable(jetsBasetoTagMatchingGeo[i], jetsBasetoTagMatchingPt[i], jetsBasetoTagMatchingHF[i]); // is (and needs to) be filled in order
    }
    for (auto i = 0; i < jetsTag.size(); i++) {
      jetsTagtoBaseMatchingTable(jetsTagtoBaseMatchingGeo[i], jetsTagtoBaseMatchingPt[i], jetsTagtoBaseMatchingHF[i]); // is (and needs to) be filled in order
    }
  }
  PROCESS_SWITCH(JetMatchingMcSub, processJets, "Perform jet matching", true);
};

extern template struct JetMatchingMcSub<o2::soa::Join<o2::aod::ChargedMCDetectorLevelJets, o2::aod::ChargedMCDetectorLevelJetConstituents>,
                                        o2::soa::Join<o2::aod::ChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                        o2::aod::ChargedMCDetectorLevelJetsMatchedToChargedMCDetectorLevelEventWiseSubtractedJets,
                                        o2::aod::ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToChargedMCDetectorLevelJets,
                                        o2::aod::JDummys>;

extern template struct JetMatchingMcSub<o2::soa::Join<o2::aod::D0ChargedMCDetectorLevelJets, o2::aod::D0ChargedMCDetectorLevelJetConstituents>,
                                        o2::soa::Join<o2::aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                        o2::aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCDetectorLevelEventWiseSubtractedJets,
                                        o2::aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToD0ChargedMCDetectorLevelJets,
                                        o2::aod::CandidatesD0MCD>;

extern template struct JetMatchingMcSub<o2::soa::Join<o2::aod::DplusChargedMCDetectorLevelJets, o2::aod::DplusChargedMCDetectorLevelJetConstituents>,
                                        o2::soa::Join<o2::aod::DplusChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DplusChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                        o2::aod::DplusChargedMCDetectorLevelJetsMatchedToDplusChargedMCDetectorLevelEventWiseSubtractedJets,
                                        o2::aod::DplusChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToDplusChargedMCDetectorLevelJets,
                                        o2::aod::CandidatesDplusMCD>;

extern template struct JetMatchingMcSub<o2::soa::Join<o2::aod::DsChargedMCDetectorLevelJets, o2::aod::DsChargedMCDetectorLevelJetConstituents>,
                                        o2::soa::Join<o2::aod::DsChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DsChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                        o2::aod::DsChargedMCDetectorLevelJetsMatchedToDsChargedMCDetectorLevelEventWiseSubtractedJets,
                                        o2::aod::DsChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToDsChargedMCDetectorLevelJets,
                                        o2::aod::CandidatesDsMCD>;

extern template struct JetMatchingMcSub<o2::soa::Join<o2::aod::DstarChargedMCDetectorLevelJets, o2::aod::DstarChargedMCDetectorLevelJetConstituents>,
                                        o2::soa::Join<o2::aod::DstarChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DstarChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                        o2::aod::DstarChargedMCDetectorLevelJetsMatchedToDstarChargedMCDetectorLevelEventWiseSubtractedJets,
                                        o2::aod::DstarChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToDstarChargedMCDetectorLevelJets,
                                        o2::aod::CandidatesDstarMCD>;

extern template struct JetMatchingMcSub<o2::soa::Join<o2::aod::LcChargedMCDetectorLevelJets, o2::aod::LcChargedMCDetectorLevelJetConstituents>,
                                        o2::soa::Join<o2::aod::LcChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::LcChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                        o2::aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCDetectorLevelEventWiseSubtractedJets,
                                        o2::aod::LcChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToLcChargedMCDetectorLevelJets,
                                        o2::aod::CandidatesLcMCD>;

extern template struct JetMatchingMcSub<o2::soa::Join<o2::aod::B0ChargedMCDetectorLevelJets, o2::aod::B0ChargedMCDetectorLevelJetConstituents>,
                                        o2::soa::Join<o2::aod::B0ChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::B0ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                        o2::aod::B0ChargedMCDetectorLevelJetsMatchedToB0ChargedMCDetectorLevelEventWiseSubtractedJets,
                                        o2::aod::B0ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToB0ChargedMCDetectorLevelJets,
                                        o2::aod::CandidatesB0MCD>;

extern template struct JetMatchingMcSub<o2::soa::Join<o2::aod::BplusChargedMCDetectorLevelJets, o2::aod::BplusChargedMCDetectorLevelJetConstituents>,
                                        o2::soa::Join<o2::aod::BplusChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::BplusChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                        o2::aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCDetectorLevelEventWiseSubtractedJets,
                                        o2::aod::BplusChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToBplusChargedMCDetectorLevelJets,
                                        o2::aod::CandidatesBplusMCD>;

extern template struct JetMatchingMcSub<o2::soa::Join<o2::aod::XicToXiPiPiChargedMCDetectorLevelJets, o2::aod::XicToXiPiPiChargedMCDetectorLevelJetConstituents>,
                                        o2::soa::Join<o2::aod::XicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::XicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                        o2::aod::XicToXiPiPiChargedMCDetectorLevelJetsMatchedToXicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJets,
                                        o2::aod::XicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToXicToXiPiPiChargedMCDetectorLevelJets,
                                        o2::aod::CandidatesXicToXiPiPiMCD>;

extern template struct JetMatchingMcSub<o2::soa::Join<o2::aod::DielectronChargedMCDetectorLevelJets, o2::aod::DielectronChargedMCDetectorLevelJetConstituents>,
                                        o2::soa::Join<o2::aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJets, o2::aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                        o2::aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCDetectorLevelEventWiseSubtractedJets,
                                        o2::aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToDielectronChargedMCDetectorLevelJets,
                                        o2::aod::CandidatesDielectronMCD>;

#endif // PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGMCSUB_H_
