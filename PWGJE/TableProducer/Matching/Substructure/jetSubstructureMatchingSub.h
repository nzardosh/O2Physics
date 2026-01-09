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

// jet analysis tasks (subscribing to jet finder task)
// this file should be removed once we can pass coloumns as templates!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#ifndef PWGJE_TABLEPRODUCER_MATCHING_SUBSTRUCTURE_JETSUBSTRUCTUREMATCHINGSUB_H_
#define PWGJE_TABLEPRODUCER_MATCHING_SUBSTRUCTURE_JETSUBSTRUCTUREMATCHINGSUB_H_

#include "PWGJE/Core/JetMatchingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>

#ifndef O2_NO_WORKFLOW_MAIN
#include <Framework/runDataProcessing.h> // IWYU pragma: export
#endif

#include <vector>

template <typename JetsBase, typename JetsTag, typename SplittingsBasetoTagMatchingTable, typename SplittingsTagtoBaseMatchingTable, typename PairsBasetoTagMatchingTable, typename PairsTagtoBaseMatchingTable, typename SplittingsBase, typename SplittingsTag, typename PairsBase, typename PairsTag, typename Candidates, typename TracksBase, typename TracksTag, typename ClustersBase>
struct JetSubstructureMatchingSub {

  o2::framework::Produces<SplittingsBasetoTagMatchingTable> splittingsBasetoTagMatchingTable;
  o2::framework::Produces<SplittingsTagtoBaseMatchingTable> splittingsTagtoBaseMatchingTable;

  o2::framework::Produces<PairsBasetoTagMatchingTable> pairsBasetoTagMatchingTable;
  o2::framework::Produces<PairsTagtoBaseMatchingTable> pairsTagtoBaseMatchingTable;

  o2::framework::Configurable<bool> doMatchingGeo{"doMatchingGeo", true, "Enable geometric matching"};
  o2::framework::Configurable<bool> doMatchingPt{"doMatchingPt", true, "Enable pt matching"};
  o2::framework::Configurable<bool> doMatchingHf{"doMatchingHf", false, "Enable HF matching"};
  o2::framework::Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.24f, "Max matching distance"};
  o2::framework::Configurable<float> minPtFraction{"minPtFraction", 0.5f, "Minimum pt fraction for pt matching"};
  o2::framework::Configurable<bool> requireGeoMatchedJets{"requireGeoMatchedJets", false, "require jets are geo matched as well"};
  o2::framework::Configurable<bool> requirePtMatchedJets{"requirePtMatchedJets", false, "require jets are pT matched as well"};
  o2::framework::Configurable<bool> requireHFMatchedJets{"requireHFMatchedJets", false, "require jets are HF matched as well"};

  static constexpr bool jetsBaseIsMc = o2::soa::relatedByIndex<o2::aod::JetMcCollisions, JetsBase>();
  static constexpr bool jetsTagIsMc = o2::soa::relatedByIndex<o2::aod::JetMcCollisions, JetsTag>();

  void init(o2::framework::InitContext const&)
  {
  }

  o2::framework::PresliceOptional<o2::aod::ChargedSPs> BaseSplittingsPerBaseJetInclusive = o2::aod::chargedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::ChargedEventWiseSubtractedSPs> TagSplittingsPerTagJetInclusive = o2::aod::chargedeventwisesubtractedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::D0ChargedSPs> BaseSplittingsPerBaseJetD0 = o2::aod::d0chargedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::D0ChargedEventWiseSubtractedSPs> TagSplittingsPerTagJetD0 = o2::aod::d0chargedeventwisesubtractedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DplusChargedSPs> BaseSplittingsPerBaseJetDplus = o2::aod::dpluschargedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DplusChargedEventWiseSubtractedSPs> TagSplittingsPerTagJetDplus = o2::aod::dpluschargedeventwisesubtractedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DsChargedSPs> BaseSplittingsPerBaseJetDs = o2::aod::dschargedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DsChargedEventWiseSubtractedSPs> TagSplittingsPerTagJetDs = o2::aod::dschargedeventwisesubtractedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DstarChargedSPs> BaseSplittingsPerBaseJetDstar = o2::aod::dstarchargedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DstarChargedEventWiseSubtractedSPs> TagSplittingsPerTagJetDstar = o2::aod::dstarchargedeventwisesubtractedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::LcChargedSPs> BaseSplittingsPerBaseJetLc = o2::aod::lcchargedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::LcChargedEventWiseSubtractedSPs> TagSplittingsPerTagJetLc = o2::aod::lcchargedeventwisesubtractedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::B0ChargedSPs> BaseSplittingsPerBaseJetB0 = o2::aod::b0chargedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::B0ChargedEventWiseSubtractedSPs> TagSplittingsPerTagJetB0 = o2::aod::b0chargedeventwisesubtractedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::BplusChargedSPs> BaseSplittingsPerBaseJetBplus = o2::aod::bpluschargedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::BplusChargedEventWiseSubtractedSPs> TagSplittingsPerTagJetBplus = o2::aod::bpluschargedeventwisesubtractedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::XicToXiPiPiChargedSPs> BaseSplittingsPerBaseJetXicToXiPiPi = o2::aod::xictoxipipichargedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::XicToXiPiPiChargedEventWiseSubtractedSPs> TagSplittingsPerTagJetXicToXiPiPi = o2::aod::xictoxipipichargedeventwisesubtractedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DielectronChargedSPs> BaseSplittingsPerBaseJetDielectron = o2::aod::dielectronchargedsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DielectronChargedEventWiseSubtractedSPs> TagSplittingsPerTagJetDielectron = o2::aod::dielectronchargedeventwisesubtractedsplitting::jetId;

  o2::framework::PresliceOptional<o2::aod::ChargedPRs> BasePairsPerBaseJetInclusive = o2::aod::chargedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::ChargedEventWiseSubtractedPRs> TagPairsPerTagJetInclusive = o2::aod::chargedeventwisesubtractedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::D0ChargedPRs> BasePairsPerBaseJetD0 = o2::aod::d0chargedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::D0ChargedEventWiseSubtractedPRs> TagPairsPerTagJetD0 = o2::aod::d0chargedeventwisesubtractedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DplusChargedPRs> BasePairsPerBaseJetDplus = o2::aod::dpluschargedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DplusChargedEventWiseSubtractedPRs> TagPairsPerTagJetDplus = o2::aod::dpluschargedeventwisesubtractedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DsChargedPRs> BasePairsPerBaseJetDs = o2::aod::dschargedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DsChargedEventWiseSubtractedPRs> TagPairsPerTagJetDs = o2::aod::dschargedeventwisesubtractedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DstarChargedPRs> BasePairsPerBaseJetDstar = o2::aod::dstarchargedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DstarChargedEventWiseSubtractedPRs> TagPairsPerTagJetDstar = o2::aod::dstarchargedeventwisesubtractedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::LcChargedPRs> BasePairsPerBaseJetLc = o2::aod::lcchargedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::LcChargedEventWiseSubtractedPRs> TagPairsPerTagJetLc = o2::aod::lcchargedeventwisesubtractedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::B0ChargedPRs> BasePairsPerBaseJetB0 = o2::aod::b0chargedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::B0ChargedEventWiseSubtractedPRs> TagPairsPerTagJetB0 = o2::aod::b0chargedeventwisesubtractedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::BplusChargedPRs> BasePairsPerBaseJetBplus = o2::aod::bpluschargedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::BplusChargedEventWiseSubtractedPRs> TagPairsPerTagJetBplus = o2::aod::bpluschargedeventwisesubtractedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::XicToXiPiPiChargedPRs> BasePairsPerBaseJetXicToXiPiPi = o2::aod::xictoxipipichargedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::XicToXiPiPiChargedEventWiseSubtractedPRs> TagPairsPerTagJetXicToXiPiPi = o2::aod::xictoxipipichargedeventwisesubtractedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DielectronChargedPRs> BasePairsPerBaseJetDielectron = o2::aod::dielectronchargedpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DielectronChargedEventWiseSubtractedPRs> TagPairsPerTagJetDielectron = o2::aod::dielectronchargedeventwisesubtractedpair::jetId;

  // workaround till binding nodes can be passed as template arguments
  template <typename CandidateTable, typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename Q, typename R, typename S, typename A, typename B>
  auto slicedPerJetForMatching(T const& table, U const& jet, V const& perIncluisveJet, M const& perD0Jet, N const& perDplusJet, O const& perDsJet, P const& perDstarJet, Q const& perLcJet, R const& perB0Jet, S const& perBplusJet, A const& perXicToXiPiPiJet, B const& perDielectronJet)
  {
    if constexpr (jethfutilities::isHFTable<CandidateTable>() || jethfutilities::isHFMcTable<CandidateTable>()) {
      return jethfutilities::slicedPerHFJet<CandidateTable>(table, jet, perD0Jet, perDplusJet, perDsJet, perDstarJet, perLcJet, perB0Jet, perBplusJet, perXicToXiPiPiJet);
    } else if constexpr (jetdqutilities::isDielectronTable<CandidateTable>() || jetdqutilities::isDielectronMcTable<CandidateTable>()) {
      return jetdqutilities::slicedPerDielectronJet<CandidateTable>(table, jet, perDielectronJet);
    } else {
      return table.sliceBy(perIncluisveJet, jet.globalIndex());
    }
  }

  template <typename T>
  auto defaultMatchedJets(T const& jetTag)
  {
    if constexpr (jetcandidateutilities::isCandidateTable<Candidates>()) {
      return jetTag.template matchedJetCand_as<JetsBase>();
    } else {
      return jetTag.template matchedJetGeo_as<JetsBase>();
    }
  }

  void processData(JetsTag const& jetsTag,
                   JetsBase const&,
                   SplittingsBase const& jetsBaseSplittings,
                   SplittingsTag const& jetsTagSplittings,
                   PairsBase const& jetsBasePairs,
                   PairsTag const& jetsTagPairs,
                   Candidates const& candidates,
                   ClustersBase const& clustersBase,
                   TracksBase const& tracksBase, TracksTag const& tracksTag)
  {

    std::vector<std::vector<int>> jetsBasetoTagSplittingsMatchingGeo, jetsBasetoTagSplittingsMatchingPt, jetsBasetoTagSplittingsMatchingHF;
    jetsBasetoTagSplittingsMatchingGeo.assign(jetsBaseSplittings.size(), {});
    jetsBasetoTagSplittingsMatchingPt.assign(jetsBaseSplittings.size(), {});
    jetsBasetoTagSplittingsMatchingHF.assign(jetsBaseSplittings.size(), {});

    std::vector<std::vector<int>> jetsTagtoBaseSplittingsMatchingGeo, jetsTagtoBaseSplittingsMatchingPt, jetsTagtoBaseSplittingsMatchingHF;
    jetsTagtoBaseSplittingsMatchingGeo.assign(jetsTagSplittings.size(), {});
    jetsTagtoBaseSplittingsMatchingPt.assign(jetsTagSplittings.size(), {});
    jetsTagtoBaseSplittingsMatchingHF.assign(jetsTagSplittings.size(), {});

    std::vector<int> jetTagSplittingsMap;
    std::vector<int> jetBaseSplittingsMap;
    jetTagSplittingsMap.resize(jetsTagSplittings.size(), -1);
    jetBaseSplittingsMap.resize(jetsBaseSplittings.size(), -1);

    std::vector<std::vector<int>> jetsBasetoTagPairsMatching;
    std::vector<std::vector<int>> jetsTagtoBasePairsMatching;
    jetsBasetoTagPairsMatching.assign(jetsBasePairs.size(), {});
    jetsTagtoBasePairsMatching.assign(jetsTagPairs.size(), {});

    std::vector<int> jetTagPairsMap;
    std::vector<int> jetBasePairsMap;
    jetTagPairsMap.resize(jetsTagPairs.size(), -1);
    jetBasePairsMap.resize(jetsBasePairs.size(), -1);

    for (auto jetTag : jetsTag) {
      bool hasMatchedJet = false;
      if constexpr (jetcandidateutilities::isCandidateTable<Candidates>()) {
        hasMatchedJet = jetTag.has_matchedJetCand();
      } else {
        hasMatchedJet = jetTag.has_matchedJetGeo();
      }
      if (hasMatchedJet) {
        // auto const& jetTagSplittings = jetsTagSplittings.sliceBy(TagSplittingsPerTagJet, jetTag.globalIndex());
        auto const& jetTagSplittings = slicedPerJetForMatching<Candidates>(jetsTagSplittings, jetTag, TagSplittingsPerTagJetInclusive, TagSplittingsPerTagJetD0, TagSplittingsPerTagJetDplus, TagSplittingsPerTagJetDs, TagSplittingsPerTagJetDstar, TagSplittingsPerTagJetLc, TagSplittingsPerTagJetB0, TagSplittingsPerTagJetBplus, TagSplittingsPerTagJetXicToXiPiPi, TagSplittingsPerTagJetDielectron);
        int tagSplittingIndex = 0;
        for (auto const& jetTagSplitting : jetTagSplittings) {
          jetTagSplittingsMap[jetTagSplitting.globalIndex()] = tagSplittingIndex;
          tagSplittingIndex++;
        }
        // auto const& jetTagPairs = jetsTagPairs.sliceBy(TagPairsPerTagJet, jetTag.globalIndex());
        auto const& jetTagPairs = slicedPerJetForMatching<Candidates>(jetsTagPairs, jetTag, TagPairsPerTagJetInclusive, TagPairsPerTagJetD0, TagPairsPerTagJetDplus, TagPairsPerTagJetDs, TagPairsPerTagJetDstar, TagPairsPerTagJetLc, TagPairsPerTagJetB0, TagPairsPerTagJetBplus, TagPairsPerTagJetXicToXiPiPi, TagPairsPerTagJetDielectron);
        int tagPairIndex = 0;
        for (auto const& jetTagPair : jetTagPairs) {
          jetTagPairsMap[jetTagPair.globalIndex()] = tagPairIndex;
          tagPairIndex++;
        }
        for (auto& jetBase : defaultMatchedJets(jetTag)) {
          if (requireGeoMatchedJets) {
            bool jetsMatchedWithGeo = false;
            for (auto& jetBaseForMatchGeo : jetTag.template matchedJetGeo_as<JetsBase>()) {
              if (jetBaseForMatchGeo.globalIndex() == jetBase.globalIndex()) {
                jetsMatchedWithGeo = true;
              }
            }
            if (!jetsMatchedWithGeo) {
              continue;
            }
          }
          if (requirePtMatchedJets) {
            bool jetsMatchedWithPt = false;
            for (auto& jetBaseForMatchPt : jetTag.template matchedJetPt_as<JetsBase>()) {
              if (jetBaseForMatchPt.globalIndex() == jetBase.globalIndex()) {
                jetsMatchedWithPt = true;
              }
            }
            if (!jetsMatchedWithPt) {
              continue;
            }
          }
          if (requireHFMatchedJets) {
            bool jetsMatchedWithHF = false;
            for (auto& jetBaseForMatchHF : jetTag.template matchedJetCand_as<JetsBase>()) {
              if (jetBaseForMatchHF.globalIndex() == jetBase.globalIndex()) {
                jetsMatchedWithHF = true;
              }
            }
            if (!jetsMatchedWithHF) {
              continue;
            }
          }
          // auto const& jetBaseSplittings = jetsBaseSplittings.sliceBy(BaseSplittingsPerBaseJet, jetBase.globalIndex());
          auto const& jetBaseSplittings = slicedPerJetForMatching<Candidates>(jetsBaseSplittings, jetBase, BaseSplittingsPerBaseJetInclusive, BaseSplittingsPerBaseJetD0, BaseSplittingsPerBaseJetDplus, BaseSplittingsPerBaseJetDs, BaseSplittingsPerBaseJetDstar, BaseSplittingsPerBaseJetLc, BaseSplittingsPerBaseJetB0, BaseSplittingsPerBaseJetBplus, BaseSplittingsPerBaseJetXicToXiPiPi, BaseSplittingsPerBaseJetDielectron);
          int baseSplittingIndex = 0;
          for (auto const& jetBaseSplitting : jetBaseSplittings) {
            jetBaseSplittingsMap[jetBaseSplitting.globalIndex()] = baseSplittingIndex;
            baseSplittingIndex++;
          }
          jetmatchingutilities::doAllMatching<jetsBaseIsMc, jetsTagIsMc>(jetBaseSplittings, jetTagSplittings, jetsBasetoTagSplittingsMatchingGeo, jetsBasetoTagSplittingsMatchingPt, jetsBasetoTagSplittingsMatchingHF, jetsTagtoBaseSplittingsMatchingGeo, jetsTagtoBaseSplittingsMatchingPt, jetsTagtoBaseSplittingsMatchingHF, candidates, tracksBase, clustersBase, candidates, tracksTag, tracksTag, doMatchingGeo, doMatchingHf, doMatchingPt, maxMatchingDistance, minPtFraction);
          // auto const& jetBasePairs = jetsBasePairs.sliceBy(BasePairsPerBaseJet, jetBase.globalIndex());
          auto const& jetBasePairs = slicedPerJetForMatching<Candidates>(jetsBasePairs, jetBase, BasePairsPerBaseJetInclusive, BasePairsPerBaseJetD0, BasePairsPerBaseJetDplus, BasePairsPerBaseJetDs, BasePairsPerBaseJetDstar, BasePairsPerBaseJetLc, BasePairsPerBaseJetB0, BasePairsPerBaseJetBplus, BasePairsPerBaseJetXicToXiPiPi, BasePairsPerBaseJetDielectron);
          int basePairIndex = 0;
          for (auto const& jetBasePair : jetBasePairs) {
            jetBasePairsMap[jetBasePair.globalIndex()] = basePairIndex;
            basePairIndex++;
          }
          jetmatchingutilities::doPairMatching<jetsBaseIsMc, jetsTagIsMc>(jetBasePairs, jetTagPairs, jetsBasetoTagPairsMatching, jetsTagtoBasePairsMatching, candidates, tracksBase, candidates, tracksTag);
        }
      }
    }

    for (auto jetsTagSplitting : jetsTagSplittings) {
      std::vector<int> tagToBaseMatchingGeoIndex;
      std::vector<int> tagToBaseMatchingPtIndex;
      std::vector<int> tagToBaseMatchingHFIndex;
      for (auto jetBaseSplittingIndex : jetsTagtoBaseSplittingsMatchingGeo[jetsTagSplitting.globalIndex()]) {
        tagToBaseMatchingGeoIndex.push_back(jetBaseSplittingsMap[jetBaseSplittingIndex]);
      }
      for (auto jetBaseSplittingIndex : jetsTagtoBaseSplittingsMatchingPt[jetsTagSplitting.globalIndex()]) {
        tagToBaseMatchingPtIndex.push_back(jetBaseSplittingsMap[jetBaseSplittingIndex]);
      }
      for (auto jetBaseSplittingIndex : jetsTagtoBaseSplittingsMatchingHF[jetsTagSplitting.globalIndex()]) {
        tagToBaseMatchingHFIndex.push_back(jetBaseSplittingsMap[jetBaseSplittingIndex]);
      }
      splittingsTagtoBaseMatchingTable(tagToBaseMatchingGeoIndex, tagToBaseMatchingPtIndex, tagToBaseMatchingHFIndex);
    }
    for (auto jetsBaseSplitting : jetsBaseSplittings) {
      std::vector<int> baseToTagMatchingGeoIndex;
      std::vector<int> baseToTagMatchingPtIndex;
      std::vector<int> baseToTagMatchingHFIndex;
      for (auto jetTagSplittingIndex : jetsBasetoTagSplittingsMatchingGeo[jetsBaseSplitting.globalIndex()]) {
        baseToTagMatchingGeoIndex.push_back(jetTagSplittingsMap[jetTagSplittingIndex]);
      }
      for (auto jetTagSplittingIndex : jetsBasetoTagSplittingsMatchingPt[jetsBaseSplitting.globalIndex()]) {
        baseToTagMatchingPtIndex.push_back(jetTagSplittingsMap[jetTagSplittingIndex]);
      }
      for (auto jetTagSplittingIndex : jetsBasetoTagSplittingsMatchingHF[jetsBaseSplitting.globalIndex()]) {
        baseToTagMatchingHFIndex.push_back(jetTagSplittingsMap[jetTagSplittingIndex]);
      }
      splittingsBasetoTagMatchingTable(baseToTagMatchingGeoIndex, baseToTagMatchingPtIndex, baseToTagMatchingHFIndex);
    }
    for (auto jetsTagPair : jetsTagPairs) {
      std::vector<int> tagToBaseMatchingIndex;
      for (auto jetBasePairIndex : jetsTagtoBasePairsMatching[jetsTagPair.globalIndex()]) {
        tagToBaseMatchingIndex.push_back(jetBasePairsMap[jetBasePairIndex]);
      }
      pairsTagtoBaseMatchingTable(tagToBaseMatchingIndex);
    }
    for (auto jetsBasePair : jetsBasePairs) {
      std::vector<int> baseToTagMatchingIndex;
      for (auto jetTagPairIndex : jetsBasetoTagPairsMatching[jetsBasePair.globalIndex()]) {
        baseToTagMatchingIndex.push_back(jetTagPairsMap[jetTagPairIndex]);
      }
      pairsBasetoTagMatchingTable(baseToTagMatchingIndex);
    }
  }
  PROCESS_SWITCH(JetSubstructureMatchingSub, processData, "charged jet substructure", true);
};

extern template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::ChargedJets, o2::aod::ChargedJetConstituents, o2::aod::ChargedJetsMatchedToChargedEventWiseSubtractedJets>,
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

extern template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::D0ChargedJets, o2::aod::D0ChargedJetConstituents, o2::aod::D0ChargedJetsMatchedToD0ChargedEventWiseSubtractedJets>,
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

extern template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::DplusChargedJets, o2::aod::DplusChargedJetConstituents, o2::aod::DplusChargedJetsMatchedToDplusChargedEventWiseSubtractedJets>,
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

extern template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::DsChargedJets, o2::aod::DsChargedJetConstituents, o2::aod::DsChargedJetsMatchedToDsChargedEventWiseSubtractedJets>,
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

extern template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::DstarChargedJets, o2::aod::DstarChargedJetConstituents, o2::aod::DstarChargedJetsMatchedToDstarChargedEventWiseSubtractedJets>,
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

extern template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::LcChargedJets, o2::aod::LcChargedJetConstituents, o2::aod::LcChargedJetsMatchedToLcChargedEventWiseSubtractedJets>,
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

extern template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::B0ChargedJets, o2::aod::B0ChargedJetConstituents, o2::aod::B0ChargedJetsMatchedToB0ChargedEventWiseSubtractedJets>,
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

extern template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::BplusChargedJets, o2::aod::BplusChargedJetConstituents, o2::aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets>,
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

extern template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::XicToXiPiPiChargedJets, o2::aod::XicToXiPiPiChargedJetConstituents, o2::aod::XicToXiPiPiChargedJetsMatchedToXicToXiPiPiChargedEventWiseSubtractedJets>,
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

extern template struct JetSubstructureMatchingSub<o2::soa::Join<o2::aod::DielectronChargedJets, o2::aod::DielectronChargedJetConstituents, o2::aod::DielectronChargedJetsMatchedToDielectronChargedEventWiseSubtractedJets>,
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

#endif // PWGJE_TABLEPRODUCER_MATCHING_SUBSTRUCTURE_JETSUBSTRUCTUREMATCHINGSUB_H_
