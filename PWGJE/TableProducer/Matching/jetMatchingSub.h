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

/// \file jetmatching.cxx
/// \brief matching event-wise constituent subtracted data jets and unsubtracted data jets
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGSUB_H_
#define PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGSUB_H_

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

template <typename JetsBase, typename JetsTag, typename JetsBasetoTagMatchingTable, typename JetsTagtoBaseMatchingTable, typename TracksTag, typename Candidates>
struct JetMatchingSub {

  o2::framework::Configurable<bool> doMatchingGeo{"doMatchingGeo", true, "Enable geometric matching"};
  o2::framework::Configurable<bool> doMatchingPt{"doMatchingPt", true, "Enable pt matching"};
  o2::framework::Configurable<bool> doMatchingHf{"doMatchingHf", false, "Enable HF matching"};
  o2::framework::Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.24f, "Max matching distance"};
  o2::framework::Configurable<float> minPtFraction{"minPtFraction", 0.5f, "Minimum pt fraction for pt matching"};

  o2::framework::Produces<JetsBasetoTagMatchingTable> jetsBasetoTagMatchingTable;
  o2::framework::Produces<JetsTagtoBaseMatchingTable> jetsTagtoBaseMatchingTable;

  // preslicing jet collections, only for Mc-based collection
  static constexpr bool jetsBaseIsMc = o2::soa::relatedByIndex<o2::aod::JMcCollisions, JetsBase>();
  static constexpr bool jetsTagIsMc = o2::soa::relatedByIndex<o2::aod::JMcCollisions, JetsTag>();

  o2::framework::Preslice<JetsBase> baseJetsPerCollision = jetsBaseIsMc ? o2::aod::jet::mcCollisionId : o2::aod::jet::collisionId;
  o2::framework::Preslice<JetsTag> tagJetsPerCollision = jetsTagIsMc ? o2::aod::jet::mcCollisionId : o2::aod::jet::collisionId;

  void init(o2::framework::InitContext const&)
  {
  }

  void processJets(o2::aod::JetCollisions const& collisions,
                   JetsBase const& jetsBase, JetsTag const& jetsTag,
                   o2::aod::JetTracks const& tracks, TracksTag const& tracksSub, Candidates const& candidates)
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
  PROCESS_SWITCH(JetMatchingSub, processJets, "Perform jet matching", true);
};

extern template struct JetMatchingSub<o2::soa::Join<o2::aod::ChargedJets, o2::aod::ChargedJetConstituents>,
                                      o2::soa::Join<o2::aod::ChargedEventWiseSubtractedJets, o2::aod::ChargedEventWiseSubtractedJetConstituents>,
                                      o2::aod::ChargedJetsMatchedToChargedEventWiseSubtractedJets,
                                      o2::aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets,
                                      o2::aod::JTrackSubs,
                                      o2::aod::JDummys>;

extern template struct JetMatchingSub<o2::soa::Join<o2::aod::D0ChargedJets, o2::aod::D0ChargedJetConstituents>,
                                      o2::soa::Join<o2::aod::D0ChargedEventWiseSubtractedJets, o2::aod::D0ChargedEventWiseSubtractedJetConstituents>,
                                      o2::aod::D0ChargedJetsMatchedToD0ChargedEventWiseSubtractedJets,
                                      o2::aod::D0ChargedEventWiseSubtractedJetsMatchedToD0ChargedJets,
                                      o2::aod::JTrackD0Subs,
                                      o2::aod::CandidatesD0Data>;

extern template struct JetMatchingSub<o2::soa::Join<o2::aod::DplusChargedJets, o2::aod::DplusChargedJetConstituents>,
                                      o2::soa::Join<o2::aod::DplusChargedEventWiseSubtractedJets, o2::aod::DplusChargedEventWiseSubtractedJetConstituents>,
                                      o2::aod::DplusChargedJetsMatchedToDplusChargedEventWiseSubtractedJets,
                                      o2::aod::DplusChargedEventWiseSubtractedJetsMatchedToDplusChargedJets,
                                      o2::aod::JTrackDplusSubs,
                                      o2::aod::CandidatesDplusData>;

extern template struct JetMatchingSub<o2::soa::Join<o2::aod::DsChargedJets, o2::aod::DsChargedJetConstituents>,
                                      o2::soa::Join<o2::aod::DsChargedEventWiseSubtractedJets, o2::aod::DsChargedEventWiseSubtractedJetConstituents>,
                                      o2::aod::DsChargedJetsMatchedToDsChargedEventWiseSubtractedJets,
                                      o2::aod::DsChargedEventWiseSubtractedJetsMatchedToDsChargedJets,
                                      o2::aod::JTrackDsSubs,
                                      o2::aod::CandidatesDsData>;

extern template struct JetMatchingSub<o2::soa::Join<o2::aod::DstarChargedJets, o2::aod::DstarChargedJetConstituents>,
                                      o2::soa::Join<o2::aod::DstarChargedEventWiseSubtractedJets, o2::aod::DstarChargedEventWiseSubtractedJetConstituents>,
                                      o2::aod::DstarChargedJetsMatchedToDstarChargedEventWiseSubtractedJets,
                                      o2::aod::DstarChargedEventWiseSubtractedJetsMatchedToDstarChargedJets,
                                      o2::aod::JTrackDstarSubs,
                                      o2::aod::CandidatesDstarData>;

extern template struct JetMatchingSub<o2::soa::Join<o2::aod::LcChargedJets, o2::aod::LcChargedJetConstituents>,
                                      o2::soa::Join<o2::aod::LcChargedEventWiseSubtractedJets, o2::aod::LcChargedEventWiseSubtractedJetConstituents>,
                                      o2::aod::LcChargedJetsMatchedToLcChargedEventWiseSubtractedJets,
                                      o2::aod::LcChargedEventWiseSubtractedJetsMatchedToLcChargedJets,
                                      o2::aod::JTrackLcSubs,
                                      o2::aod::CandidatesLcData>;

extern template struct JetMatchingSub<o2::soa::Join<o2::aod::B0ChargedJets, o2::aod::B0ChargedJetConstituents>,
                                      o2::soa::Join<o2::aod::B0ChargedEventWiseSubtractedJets, o2::aod::B0ChargedEventWiseSubtractedJetConstituents>,
                                      o2::aod::B0ChargedJetsMatchedToB0ChargedEventWiseSubtractedJets,
                                      o2::aod::B0ChargedEventWiseSubtractedJetsMatchedToB0ChargedJets,
                                      o2::aod::JTrackB0Subs,
                                      o2::aod::CandidatesB0Data>;

extern template struct JetMatchingSub<o2::soa::Join<o2::aod::BplusChargedJets, o2::aod::BplusChargedJetConstituents>,
                                      o2::soa::Join<o2::aod::BplusChargedEventWiseSubtractedJets, o2::aod::BplusChargedEventWiseSubtractedJetConstituents>,
                                      o2::aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets,
                                      o2::aod::BplusChargedEventWiseSubtractedJetsMatchedToBplusChargedJets,
                                      o2::aod::JTrackBplusSubs,
                                      o2::aod::CandidatesBplusData>;

extern template struct JetMatchingSub<o2::soa::Join<o2::aod::XicToXiPiPiChargedJets, o2::aod::XicToXiPiPiChargedJetConstituents>,
                                      o2::soa::Join<o2::aod::XicToXiPiPiChargedEventWiseSubtractedJets, o2::aod::XicToXiPiPiChargedEventWiseSubtractedJetConstituents>,
                                      o2::aod::XicToXiPiPiChargedJetsMatchedToXicToXiPiPiChargedEventWiseSubtractedJets,
                                      o2::aod::XicToXiPiPiChargedEventWiseSubtractedJetsMatchedToXicToXiPiPiChargedJets,
                                      o2::aod::JTrackXicToXiPiPiSubs,
                                      o2::aod::CandidatesXicToXiPiPiData>;

extern template struct JetMatchingSub<o2::soa::Join<o2::aod::DielectronChargedJets, o2::aod::DielectronChargedJetConstituents>,
                                      o2::soa::Join<o2::aod::DielectronChargedEventWiseSubtractedJets, o2::aod::DielectronChargedEventWiseSubtractedJetConstituents>,
                                      o2::aod::DielectronChargedJetsMatchedToDielectronChargedEventWiseSubtractedJets,
                                      o2::aod::DielectronChargedEventWiseSubtractedJetsMatchedToDielectronChargedJets,
                                      o2::aod::JTrackDielectronSubs,
                                      o2::aod::CandidatesDielectronData>;

#endif // PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGSUB_H_
