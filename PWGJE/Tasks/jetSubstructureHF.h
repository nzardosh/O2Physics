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

// heavy-flavour jet substructure task (subscribing to jet finder hf task)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#ifndef PWGJE_TASKS_JETSUBSTRUCTUREHF_H_
#define PWGJE_TASKS_JETSUBSTRUCTUREHF_H_

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDQUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetSubstructureUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>

#include <TMath.h>

#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <math.h>

// NB: runDataProcessing.h must be included after customize!

template <typename JetTableData, typename JetTableMCD, typename JetTableMCP, typename JetTableDataSub, typename CandidateTable, typename CandidateTableMCP, typename SubstructureTableData, typename RecoilTableData, typename SplittingsTableData, typename PairsTableData, typename SubstructureTableMCD, typename RecoilTableMCD, typename SplittingsTableMCD, typename PairsTableMCD, typename SubstructureTableMCP, typename RecoilTableMCP, typename SplittingsTableMCP, typename PairsTableMCP, typename SubstructureTableDataSub, typename SplittingsTableDataSub, typename PairsTableDataSub, typename TracksSub>
struct JetSubstructureHFTask {
  o2::framework::Produces<SubstructureTableData> jetSubstructureDataTable;
  o2::framework::Produces<SubstructureTableMCD> jetSubstructureMCDTable;
  o2::framework::Produces<SubstructureTableMCP> jetSubstructureMCPTable;
  o2::framework::Produces<SubstructureTableDataSub> jetSubstructureDataSubTable;

  o2::framework::Produces<SplittingsTableData> jetSplittingsDataTable;
  o2::framework::Produces<SplittingsTableMCD> jetSplittingsMCDTable;
  o2::framework::Produces<SplittingsTableMCP> jetSplittingsMCPTable;
  o2::framework::Produces<SplittingsTableDataSub> jetSplittingsDataSubTable;

  o2::framework::Produces<PairsTableData> jetPairsDataTable;
  o2::framework::Produces<PairsTableMCD> jetPairsMCDTable;
  o2::framework::Produces<PairsTableMCP> jetPairsMCPTable;
  o2::framework::Produces<PairsTableDataSub> jetPairsDataSubTable;

  o2::framework::Produces<RecoilTableData> jetRecoilDataTable;
  o2::framework::Produces<RecoilTableMCD> jetRecoilMCDTable;
  o2::framework::Produces<RecoilTableMCP> jetRecoilMCPTable;

  // Jet level configurables
  o2::framework::Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  o2::framework::Configurable<float> beta{"beta", 0.0, "soft drop beta"};
  o2::framework::Configurable<float> kappa{"kappa", 1.0, "angularity kappa"};
  o2::framework::Configurable<float> alpha{"alpha", 1.0, "angularity alpha"};
  o2::framework::Configurable<bool> doPairBkg{"doPairBkg", true, "save bkg pairs"};
  o2::framework::Configurable<float> pairConstituentPtMin{"pairConstituentPtMin", 1.0, "pt cut off for constituents going into pairs"};
  o2::framework::Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  o2::framework::Configurable<bool> applyTrackingEfficiency{"applyTrackingEfficiency", {false}, "configurable to decide whether to apply artificial tracking efficiency (discarding tracks) in jet finding"};
  o2::framework::Configurable<std::vector<double>> trackingEfficiencyPtBinning{"trackingEfficiencyPtBinning", {0., 10, 999.}, "pt binning of tracking efficiency array if applyTrackingEfficiency is true"};
  o2::framework::Configurable<std::vector<double>> trackingEfficiency{"trackingEfficiency", {1.0, 1.0}, "tracking efficiency array applied to jet finding if applyTrackingEfficiency is true"};
  o2::framework::Configurable<float> recoilRegion{"recoilRegion", 0.6, "recoil acceptance in phi"};

  o2::framework::Service<o2::framework::O2DatabasePDG> pdg;
  float candMass;

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  std::vector<float> energyMotherVec;
  std::vector<float> ptLeadingVec;
  std::vector<float> ptSubLeadingVec;
  std::vector<float> thetaVec;
  std::vector<float> nSub;
  std::vector<float> pairJetPtVec;
  std::vector<float> pairJetEnergyVec;
  std::vector<float> pairJetThetaVec;
  std::vector<float> pairJetPerpCone1PtVec;
  std::vector<float> pairJetPerpCone1EnergyVec;
  std::vector<float> pairJetPerpCone1ThetaVec;
  std::vector<float> pairPerpCone1PerpCone1PtVec;
  std::vector<float> pairPerpCone1PerpCone1EnergyVec;
  std::vector<float> pairPerpCone1PerpCone1ThetaVec;
  std::vector<float> pairPerpCone1PerpCone2PtVec;
  std::vector<float> pairPerpCone1PerpCone2EnergyVec;
  std::vector<float> pairPerpCone1PerpCone2ThetaVec;
  float angularity;
  float leadingConstituentPt;
  float perpConeRho;

  o2::framework::HistogramRegistry registry;

  int trackSelection = -1;

  void init(o2::framework::InitContext const&)
  {
    registry.add("h2_jet_pt_jet_zg", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {o2::framework::HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {o2::framework::HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_nsd", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {o2::framework::HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_part_jet_zg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{z}_{g}^{part}", {o2::framework::HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_rg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{R}_{g}^{part}", {o2::framework::HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_nsd_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{n}_{SD}^{part}", {o2::framework::HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {o2::framework::HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {o2::framework::HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {o2::framework::HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
    jetReclusterer.ghostRepeatN = 0;

    candMass = jetcandidateutilities::getTablePDGMass<CandidateTable>();

    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    if (applyTrackingEfficiency) {
      if (trackingEfficiencyPtBinning->size() < 2) {
        LOGP(fatal, "jetFinder workflow: trackingEfficiencyPtBinning configurable should have at least two bin edges");
      }
      if (trackingEfficiency->size() + 1 != trackingEfficiencyPtBinning->size()) {
        LOGP(fatal, "jetFinder workflow: trackingEfficiency configurable should have exactly one less entry than the number of bin edges set in trackingEfficiencyPtBinning configurable");
      }
    }
  }

  o2::framework::Preslice<o2::aod::JetTracks> TracksPerCollision = o2::aod::jtrack::collisionId;
  o2::framework::PresliceOptional<o2::aod::JTrackD0Subs> TracksPerD0DataSub = o2::aod::bkgd0::candidateId;
  o2::framework::PresliceOptional<o2::aod::JTrackDplusSubs> TracksPerDplusDataSub = o2::aod::bkgdplus::candidateId;
  o2::framework::PresliceOptional<o2::aod::JTrackDsSubs> TracksPerDsDataSub = o2::aod::bkgds::candidateId;
  o2::framework::PresliceOptional<o2::aod::JTrackDstarSubs> TracksPerDstarDataSub = o2::aod::bkgdstar::candidateId;
  o2::framework::PresliceOptional<o2::aod::JTrackLcSubs> TracksPerLcDataSub = o2::aod::bkglc::candidateId;
  o2::framework::PresliceOptional<o2::aod::JTrackB0Subs> TracksPerB0DataSub = o2::aod::bkgb0::candidateId;
  o2::framework::PresliceOptional<o2::aod::JTrackBplusSubs> TracksPerBplusDataSub = o2::aod::bkgbplus::candidateId;
  o2::framework::PresliceOptional<o2::aod::JTrackXicToXiPiPiSubs> TracksPerXicToXiPiPiDataSub = o2::aod::bkgxictoxipipi::candidateId;
  o2::framework::PresliceOptional<o2::aod::JTrackDielectronSubs> TracksPerDielectronDataSub = o2::aod::bkgdielectron::candidateId;
  o2::framework::Preslice<o2::aod::JetParticles> ParticlesPerMcCollision = o2::aod::jmcparticle::mcCollisionId;

  template <typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename Q, typename R>
  auto selectSlicer(T const& D0Slicer, U const& DplusSlicer, V const& DsSlicer, M const& DstarSlicer, N const& LcSlicer, O const& B0Slicer, P const& BplusSlicer, Q const& XicToXiPiPiSlicer, R const& DielectronSlicer)
  {
    if constexpr (jethfutilities::isD0Table<CandidateTable>()) {
      return D0Slicer;
    } else if constexpr (jethfutilities::isDplusTable<CandidateTable>()) {
      return DplusSlicer;
    } else if constexpr (jethfutilities::isDsTable<CandidateTable>()) {
      return DsSlicer;
    } else if constexpr (jethfutilities::isDstarTable<CandidateTable>()) {
      return DstarSlicer;
    } else if constexpr (jethfutilities::isLcTable<CandidateTable>()) {
      return LcSlicer;
    } else if constexpr (jethfutilities::isB0Table<CandidateTable>()) {
      return B0Slicer;
    } else if constexpr (jethfutilities::isBplusTable<CandidateTable>()) {
      return BplusSlicer;
    } else if constexpr (jethfutilities::isXicToXiPiPiTable<CandidateTable>()) {
      return XicToXiPiPiSlicer;
    } else if constexpr (jetdqutilities::isDielectronTable<CandidateTable>()) {
      return DielectronSlicer;
    } else {
      return D0Slicer;
    }
  }

  template <bool isMCP, bool isSubtracted, typename T, typename U>
  void jetReclustering(T const& jet, U& splittingTable, int nHFCandidates)
  {
    energyMotherVec.clear();
    ptLeadingVec.clear();
    ptSubLeadingVec.clear();
    thetaVec.clear();
    jetReclustered.clear();
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);
    fastjet::PseudoJet daughterSubJet = jetReclustered[0];
    fastjet::PseudoJet parentSubJet1;
    fastjet::PseudoJet parentSubJet2;
    bool softDropped = false;
    auto nsd = 0.0;
    while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {

      int nHFInSubjet1 = 0;
      for (auto& subjet1Constituent : parentSubJet1.constituents()) {
        if (subjet1Constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidate)) {
          nHFInSubjet1++;
        }
      }
      if (nHFCandidates == 1 && nHFInSubjet1 == 0) {
        std::swap(parentSubJet1, parentSubJet2);
      }
      if (nHFCandidates == 2) {
        if (nHFInSubjet1 == 2) {
          daughterSubJet = parentSubJet1;
          continue;
        }
        if (nHFInSubjet1 == 0) {
          daughterSubJet = parentSubJet2;
          continue;
        }
      }

      std::vector<int32_t> tracks;
      std::vector<int32_t> candidates;
      std::vector<int32_t> clusters;
      for (const auto& constituent : sorted_by_pt(parentSubJet2.constituents())) {
        if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::track)) {
          tracks.push_back(constituent.template user_info<fastjetutilities::fastjet_user_info>().getIndex());
        }
        if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidate)) {
          candidates.push_back(constituent.template user_info<fastjetutilities::fastjet_user_info>().getIndex());
        }
      }
      splittingTable(jet.globalIndex(), tracks, clusters, candidates, parentSubJet2.perp(), parentSubJet2.eta(), parentSubJet2.phi(), 0);
      auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
      auto theta = parentSubJet1.delta_R(parentSubJet2);
      energyMotherVec.push_back(daughterSubJet.e());
      ptLeadingVec.push_back(parentSubJet1.pt());
      ptSubLeadingVec.push_back(parentSubJet2.pt());
      thetaVec.push_back(theta);
      if (z >= zCut * TMath::Power(theta / (jet.r() / 100.f), beta)) {
        if (!softDropped) {
          auto zg = z;
          auto rg = theta;
          if constexpr (!isSubtracted && !isMCP) {
            registry.fill(HIST("h2_jet_pt_jet_zg"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_jet_rg"), jet.pt(), rg);
          }
          if constexpr (!isSubtracted && isMCP) {
            registry.fill(HIST("h2_jet_pt_part_jet_zg_part"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_part_jet_rg_part"), jet.pt(), rg);
          }
          if constexpr (isSubtracted && !isMCP) {
            registry.fill(HIST("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted"), jet.pt(), rg);
          }
          softDropped = true;
        }
        nsd++;
      }
      daughterSubJet = parentSubJet1;
      if (nHFCandidates == 2 && nHFInSubjet1 == 1) {
        break;
      }
    }
    if constexpr (!isSubtracted && !isMCP) {
      registry.fill(HIST("h2_jet_pt_jet_nsd"), jet.pt(), nsd);
    }
    if constexpr (!isSubtracted && isMCP) {
      registry.fill(HIST("h2_jet_pt_part_jet_nsd_part"), jet.pt(), nsd);
    }
    if constexpr (isSubtracted && !isMCP) {
      registry.fill(HIST("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted"), jet.pt(), nsd);
    }
  }

  template <bool isMC, bool isSubtracted, typename T, typename U, typename V, typename M, typename N>
  void jetPairing(T const& jet, U const& tracks, V const& /*candidates*/, M const& slicer, N& pairTable)
  {
    pairJetPtVec.clear();
    pairJetEnergyVec.clear();
    pairJetThetaVec.clear();
    std::vector<std::decay_t<typename U::iterator>> tracksVec;
    std::vector<std::decay_t<typename V::iterator>> candidatesVec;
    std::vector<int32_t> tracksVecIds;
    std::vector<int32_t> candidatesVecIds;
    for (auto& constituent : jet.template tracks_as<U>()) {
      if (constituent.pt() >= pairConstituentPtMin) {
        tracksVec.push_back(constituent);
        tracksVecIds.push_back(constituent.globalIndex());
      }
    }
    for (auto& candidate : jet.template candidates_as<V>()) {
      candidatesVec.push_back(candidate);
      candidatesVecIds.push_back(candidate.globalIndex());
    }
    if (tracksVec.size() >= 1) {
      for (typename std::vector<std::decay_t<typename U::iterator>>::size_type track1Index = 0; track1Index < tracksVec.size(); track1Index++) {
        for (typename std::vector<std::decay_t<typename U::iterator>>::size_type track2Index = track1Index + 1; track2Index < tracksVec.size(); track2Index++) {
          pairJetPtVec.push_back(tracksVec.at(track1Index).pt() * tracksVec.at(track2Index).pt());
          pairJetEnergyVec.push_back(2.0 * tracksVec.at(track1Index).energy() * tracksVec.at(track2Index).energy());
          pairJetThetaVec.push_back(jetutilities::deltaR(tracksVec.at(track1Index), tracksVec.at(track2Index)));
          pairTable(jet.globalIndex(), tracksVecIds.at(track1Index), tracksVecIds.at(track2Index), -1, -1);
        }
      }
    }
    if (candidatesVec.size() >= 1) {
      for (typename std::vector<std::decay_t<typename V::iterator>>::size_type candidate1Index = 0; candidate1Index < candidatesVec.size(); candidate1Index++) {
        for (typename std::vector<std::decay_t<typename V::iterator>>::size_type candidate2Index = candidate1Index + 1; candidate2Index < candidatesVec.size(); candidate2Index++) {
          pairJetPtVec.push_back(candidatesVec.at(candidate1Index).pt() * candidatesVec.at(candidate2Index).pt());
          auto candidate1Energy = std::sqrt((candidatesVec.at(candidate1Index).p() * candidatesVec.at(candidate1Index).p()) + (candMass * candMass));
          auto candidate2Energy = std::sqrt((candidatesVec.at(candidate2Index).p() * candidatesVec.at(candidate2Index).p()) + (candMass * candMass));
          pairJetEnergyVec.push_back(2.0 * candidate1Energy * candidate2Energy);
          pairJetThetaVec.push_back(jetutilities::deltaR(candidatesVec.at(candidate1Index), candidatesVec.at(candidate2Index)));
          pairTable(jet.globalIndex(), -1, -1, candidatesVecIds.at(candidate1Index), candidatesVecIds.at(candidate2Index));
        }
      }
    }
    if (candidatesVec.size() >= 1 && tracksVec.size() >= 1) {
      for (typename std::vector<std::decay_t<typename V::iterator>>::size_type candidateIndex = 0; candidateIndex < candidatesVec.size(); candidateIndex++) { // could just directly get the candidate and tracks here but keeping it consistent with above
        for (typename std::vector<std::decay_t<typename U::iterator>>::size_type trackIndex = 0; trackIndex < tracksVec.size(); trackIndex++) {
          pairJetPtVec.push_back(candidatesVec.at(candidateIndex).pt() * tracksVec.at(trackIndex).pt());
          auto candidateEnergy = std::sqrt((candidatesVec.at(candidateIndex).p() * candidatesVec.at(candidateIndex).p()) + (candMass * candMass));
          pairJetEnergyVec.push_back(2.0 * candidateEnergy * tracksVec.at(trackIndex).energy());
          pairJetThetaVec.push_back(jetutilities::deltaR(candidatesVec.at(candidateIndex), tracksVec.at(trackIndex)));
          pairTable(jet.globalIndex(), tracksVecIds.at(trackIndex), -1, candidatesVecIds.at(candidateIndex), -1);
        }
      }
    }

    pairJetPerpCone1PtVec.clear();
    pairJetPerpCone1EnergyVec.clear();
    pairJetPerpCone1ThetaVec.clear();
    pairPerpCone1PerpCone1PtVec.clear();
    pairPerpCone1PerpCone1EnergyVec.clear();
    pairPerpCone1PerpCone1ThetaVec.clear();
    pairPerpCone1PerpCone2PtVec.clear();
    pairPerpCone1PerpCone2EnergyVec.clear();
    pairPerpCone1PerpCone2ThetaVec.clear();
    int32_t slicerId = -1;
    if constexpr (isSubtracted) {
      auto const& candidate = jet.template candidates_first_as<V>();
      slicerId = candidate.globalIndex();
    } else {
      if constexpr (!isMC) {
        slicerId = jet.collisionId();
      } else {
        slicerId = jet.mcCollisionId();
      }
    }
    auto tracksPerCollision = tracks.sliceBy(slicer, slicerId);

    float perpCone1Phi = RecoDecay::constrainAngle<float, float>(jet.phi() + (M_PI / 2.));
    float perpCone2Phi = RecoDecay::constrainAngle<float, float>(jet.phi() - (M_PI / 2.));
    float perpCone1Pt = 0.0;
    float perpCone2Pt = 0.0;
    std::vector<typename U::iterator> tracksPerpCone1Vec;
    std::vector<typename U::iterator> tracksPerpCone2Vec;
    for (auto const& track : tracksPerCollision) {
      if (!jetderiveddatautilities::applyTrackKinematics(track)) {
        continue;
      }

      if constexpr (!std::is_same_v<std::decay_t<U>, o2::aod::JetParticles>) {
        if (!jetfindingutilities::isTrackSelected<typename U::iterator, typename U::iterator>(track, trackSelection, applyTrackingEfficiency, trackingEfficiency, trackingEfficiencyPtBinning)) {
          continue;
        }
      }
      bool isdaughterTrack = false;
      for (auto& candidate : jet.template candidates_as<V>()) {
        if (jetcandidateutilities::isDaughterTrack(track, candidate)) {
          isdaughterTrack = true;
          break;
        }
      }
      if (isdaughterTrack) {
        continue;
      }
      float deltaPhi1 = track.phi() - perpCone1Phi;
      deltaPhi1 = RecoDecay::constrainAngle<float, float>(deltaPhi1, -M_PI);
      float deltaPhi2 = track.phi() - perpCone2Phi;
      deltaPhi2 = RecoDecay::constrainAngle<float, float>(deltaPhi2, -M_PI);
      float deltaEta = jet.eta() - track.eta();

      if (TMath::Sqrt((deltaPhi1 * deltaPhi1) + (deltaEta * deltaEta)) <= jet.r() / 100.0) {
        if (track.pt() >= pairConstituentPtMin) {
          tracksPerpCone1Vec.push_back(track);
        }
        perpCone1Pt += track.pt();
      }
      if (TMath::Sqrt((deltaPhi2 * deltaPhi2) + (deltaEta * deltaEta)) <= jet.r() / 100.0) {
        if (track.pt() >= pairConstituentPtMin) {
          tracksPerpCone2Vec.push_back(track);
        }
        perpCone2Pt += track.pt();
      }
    }
    perpConeRho = (perpCone1Pt + perpCone2Pt) / (2 * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)); // currently done per jet - could be better to do for leading jet if pushing to very low pT
    if (doPairBkg) {
      if (tracksVec.size() >= 1 && tracksPerpCone1Vec.size() >= 1) {
        for (typename std::vector<typename U::iterator>::size_type track1Index = 0; track1Index < tracksVec.size(); track1Index++) {
          for (typename std::vector<typename U::iterator>::size_type track2Index = 0; track2Index < tracksPerpCone1Vec.size(); track2Index++) {
            pairJetPerpCone1PtVec.push_back(tracksVec.at(track1Index).pt() * tracksPerpCone1Vec.at(track2Index).pt());
            pairJetPerpCone1EnergyVec.push_back(2.0 * tracksVec.at(track1Index).energy() * tracksPerpCone1Vec.at(track2Index).energy());
            float dPhi = RecoDecay::constrainAngle(tracksVec.at(track1Index).phi() - (tracksPerpCone1Vec.at(track2Index).phi() - (M_PI / 2.)), -M_PI);
            float dEta = tracksVec.at(track1Index).eta() - tracksPerpCone1Vec.at(track2Index).eta();
            pairJetPerpCone1ThetaVec.push_back(std::sqrt(dEta * dEta + dPhi * dPhi));
          }
        }
      }

      if (candidatesVec.size() >= 1 && tracksPerpCone1Vec.size() >= 1) {
        for (typename std::vector<std::decay_t<typename V::iterator>>::size_type candidate1Index = 0; candidate1Index < candidatesVec.size(); candidate1Index++) {
          for (typename std::vector<typename U::iterator>::size_type track2Index = 0; track2Index < tracksPerpCone1Vec.size(); track2Index++) {
            pairJetPerpCone1PtVec.push_back(candidatesVec.at(candidate1Index).pt() * tracksPerpCone1Vec.at(track2Index).pt());
            auto candidate1Energy = std::sqrt((candidatesVec.at(candidate1Index).p() * candidatesVec.at(candidate1Index).p()) + (candMass * candMass));
            pairJetPerpCone1EnergyVec.push_back(2.0 * candidate1Energy * tracksPerpCone1Vec.at(track2Index).energy());
            float dPhi = RecoDecay::constrainAngle(candidatesVec.at(candidate1Index).phi() - (tracksPerpCone1Vec.at(track2Index).phi() - (M_PI / 2.)), -M_PI);
            float dEta = candidatesVec.at(candidate1Index).eta() - tracksPerpCone1Vec.at(track2Index).eta();
            pairJetPerpCone1ThetaVec.push_back(std::sqrt(dEta * dEta + dPhi * dPhi));
          }
        }
      }

      if (tracksPerpCone1Vec.size() >= 1) {
        for (typename std::vector<typename U::iterator>::size_type track1Index = 0; track1Index < tracksPerpCone1Vec.size(); track1Index++) {
          for (typename std::vector<typename U::iterator>::size_type track2Index = track1Index + 1; track2Index < tracksPerpCone1Vec.size(); track2Index++) {
            pairPerpCone1PerpCone1PtVec.push_back(tracksPerpCone1Vec.at(track1Index).pt() * tracksPerpCone1Vec.at(track2Index).pt());
            pairPerpCone1PerpCone1EnergyVec.push_back(2.0 * tracksPerpCone1Vec.at(track1Index).energy() * tracksPerpCone1Vec.at(track2Index).energy());
            pairPerpCone1PerpCone1ThetaVec.push_back(jetutilities::deltaR(tracksPerpCone1Vec.at(track1Index), tracksPerpCone1Vec.at(track2Index)));
          }
        }
      }

      if (tracksPerpCone1Vec.size() >= 1 && tracksPerpCone2Vec.size() >= 1) {
        for (typename std::vector<typename U::iterator>::size_type track1Index = 0; track1Index < tracksPerpCone1Vec.size(); track1Index++) {
          for (typename std::vector<typename U::iterator>::size_type track2Index = 0; track2Index < tracksPerpCone2Vec.size(); track2Index++) {
            pairPerpCone1PerpCone2PtVec.push_back(tracksPerpCone1Vec.at(track1Index).pt() * tracksPerpCone2Vec.at(track2Index).pt());
            pairPerpCone1PerpCone2EnergyVec.push_back(2.0 * tracksPerpCone1Vec.at(track1Index).energy() * tracksPerpCone2Vec.at(track2Index).energy());
            float dPhi = RecoDecay::constrainAngle((tracksPerpCone1Vec.at(track1Index).phi() - (M_PI / 2.)) - (tracksPerpCone2Vec.at(track2Index).phi() + (M_PI / 2.)), -M_PI);
            float dEta = tracksPerpCone1Vec.at(track1Index).eta() - tracksPerpCone2Vec.at(track2Index).eta();
            pairPerpCone1PerpCone2ThetaVec.push_back(std::sqrt(dEta * dEta + dPhi * dPhi));
          }
        }
      }
    }
  }

  template <typename T, typename U, typename V>
  void jetSubstructureSimple(T const& jet, U const& /*tracks*/, V const& /*candidates*/)
  {
    angularity = 0.0;
    leadingConstituentPt = 0.0;
    for (auto& candidate : jet.template candidates_as<V>()) {
      if (candidate.pt() >= leadingConstituentPt) {
        leadingConstituentPt = candidate.pt();
      }
      angularity += std::pow(candidate.pt(), kappa) * std::pow(jetutilities::deltaR(jet, candidate), alpha);
    }
    for (auto& constituent : jet.template tracks_as<U>()) {
      if (constituent.pt() >= leadingConstituentPt) {
        leadingConstituentPt = constituent.pt();
      }
      angularity += std::pow(constituent.pt(), kappa) * std::pow(jetutilities::deltaR(jet, constituent), alpha);
    }
    angularity /= (jet.pt() * (jet.r() / 100.f));
  }

  template <bool isSubtracted, typename T, typename U, typename V, typename M, typename N, typename O, typename P>
  void analyseCharged(T const& jet, U const& tracks, V const& candidates, M const& trackSlicer, N& outputTable, O& splittingTable, P& pairTable)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<U>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    int nHFCandidates = 0;
    for (auto& jetHFCandidate : jet.template candidates_as<V>()) {
      fastjetutilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidate), candMass);
      nHFCandidates++;
    }
    nSub = jetsubstructureutilities::getNSubjettiness(jet, tracks, tracks, candidates, 2, fastjet::contrib::CA_Axes(), true, zCut, beta);
    jetReclustering<false, isSubtracted>(jet, splittingTable, nHFCandidates);
    jetPairing<false, isSubtracted>(jet, tracks, candidates, trackSlicer, pairTable);
    jetSubstructureSimple(jet, tracks, candidates);
    outputTable(energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, nSub[0], nSub[1], nSub[2], pairJetPtVec, pairJetEnergyVec, pairJetThetaVec, pairJetPerpCone1PtVec, pairJetPerpCone1EnergyVec, pairJetPerpCone1ThetaVec, pairPerpCone1PerpCone1PtVec, pairPerpCone1PerpCone1EnergyVec, pairPerpCone1PerpCone1ThetaVec, pairPerpCone1PerpCone2PtVec, pairPerpCone1PerpCone2EnergyVec, pairPerpCone1PerpCone2ThetaVec, angularity, leadingConstituentPt, perpConeRho);
  }

  template <typename T, typename U, typename V, typename M>
  void analyseRecoilCharged(T const& jets, U const& /*tracks*/, V const& candidates, M& outputRecoilTable)
  {
    for (auto const& candidate : candidates) {
      for (auto const& jet : jets) {
        if (std::abs(RecoDecay::constrainAngle(jet.phi() - candidate.phi(), -M_PI)) > (M_PI - recoilRegion)) {
          outputRecoilTable(jet.globalIndex(), candidate.globalIndex(), jet.pt(), jet.phi(), jet.eta(), jet.r(), jet.tracksIds().size()); // add variables for recoil jet tagging
          break;
        }
      }
    }
  }

  void processChargedJetsData(typename JetTableData::iterator const& jet,
                              CandidateTable const& candidates,
                              o2::aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks, candidates, TracksPerCollision, jetSubstructureDataTable, jetSplittingsDataTable, jetPairsDataTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsData, "HF jet substructure on data", false);

  void processChargedJetsDataSub(typename JetTableDataSub::iterator const& jet,
                                 CandidateTable const& candidates,
                                 TracksSub const& tracks)
  {
    analyseCharged<true>(jet, tracks, candidates, selectSlicer(TracksPerD0DataSub, TracksPerDplusDataSub, TracksPerDsDataSub, TracksPerDstarDataSub, TracksPerLcDataSub, TracksPerB0DataSub, TracksPerBplusDataSub, TracksPerXicToXiPiPiDataSub, TracksPerDielectronDataSub), jetSubstructureDataSubTable, jetSplittingsDataSubTable, jetPairsDataSubTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsDataSub, "HF jet substructure on data", false);

  void processChargedJetsMCD(typename JetTableMCD::iterator const& jet,
                             CandidateTable const& candidates,
                             o2::aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks, candidates, TracksPerCollision, jetSubstructureMCDTable, jetSplittingsMCDTable, jetPairsMCDTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsMCD, "HF jet substructure on data", false);

  void processChargedJetsMCP(typename JetTableMCP::iterator const& jet,
                             o2::aod::JetParticles const& particles,
                             CandidateTableMCP const& candidates)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<o2::aod::JetParticles>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
    }
    int nHFCandidates = 0;
    for (auto& jetHFCandidate : jet.template candidates_as<CandidateTableMCP>()) {
      fastjetutilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidate), candMass);
      nHFCandidates++;
    }
    nSub = jetsubstructureutilities::getNSubjettiness(jet, particles, particles, candidates, 2, fastjet::contrib::CA_Axes(), true, zCut, beta);
    jetReclustering<true, false>(jet, jetSplittingsMCPTable, nHFCandidates);
    jetPairing<true, false>(jet, particles, candidates, ParticlesPerMcCollision, jetPairsMCPTable);
    jetSubstructureSimple(jet, particles, candidates);
    jetSubstructureMCPTable(energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, nSub[0], nSub[1], nSub[2], pairJetPtVec, pairJetEnergyVec, pairJetThetaVec, pairJetPerpCone1PtVec, pairJetPerpCone1EnergyVec, pairJetPerpCone1ThetaVec, pairPerpCone1PerpCone1PtVec, pairPerpCone1PerpCone1EnergyVec, pairPerpCone1PerpCone1ThetaVec, pairPerpCone1PerpCone2PtVec, pairPerpCone1PerpCone2EnergyVec, pairPerpCone1PerpCone2ThetaVec, angularity, leadingConstituentPt, perpConeRho);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsMCP, "HF jet substructure on MC particle level", false);

  void processChargedRecoilJetsData(o2::aod::JetCollision const& /*collision*/,
                                    o2::soa::Join<o2::aod::ChargedJets, o2::aod::ChargedJetConstituents> const& jets,
                                    CandidateTable const& candidates,
                                    o2::aod::JetTracks const& tracks)
  {
    analyseRecoilCharged(jets, tracks, candidates, jetRecoilDataTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedRecoilJetsData, "HF recoil jet data", false);

  void processChargedRecoilJetsMCD(o2::aod::JetCollision const& /*collision*/,
                                   o2::soa::Join<o2::aod::ChargedMCDetectorLevelJets, o2::aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                                   CandidateTable const& candidates,
                                   o2::aod::JetTracks const& tracks)
  {
    analyseRecoilCharged(jets, tracks, candidates, jetRecoilMCDTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedRecoilJetsMCD, "HF recoil jet mcd", false);

  void processChargedRecoilJetsMCP(o2::aod::JetCollision const& /*collision*/,
                                   o2::soa::Join<o2::aod::ChargedMCParticleLevelJets, o2::aod::ChargedMCParticleLevelJetConstituents> const& jets,
                                   CandidateTableMCP const& candidates,
                                   o2::aod::JetParticles const& tracks)
  {
    analyseRecoilCharged(jets, tracks, candidates, jetRecoilMCPTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedRecoilJetsMCP, "HF recoil jet mcp", false);
};

extern template struct JetSubstructureHFTask<o2::soa::Join<o2::aod::D0ChargedJets, o2::aod::D0ChargedJetConstituents>, o2::soa::Join<o2::aod::D0ChargedMCDetectorLevelJets, o2::aod::D0ChargedMCDetectorLevelJetConstituents>, o2::soa::Join<o2::aod::D0ChargedMCParticleLevelJets, o2::aod::D0ChargedMCParticleLevelJetConstituents>, o2::soa::Join<o2::aod::D0ChargedEventWiseSubtractedJets, o2::aod::D0ChargedEventWiseSubtractedJetConstituents>, o2::aod::CandidatesD0Data, o2::aod::CandidatesD0MCP, o2::aod::D0CJetSSs, o2::aod::D0CJetRs, o2::aod::D0ChargedSPs, o2::aod::D0ChargedPRs, o2::aod::D0CMCDJetSSs, o2::aod::D0CMCDJetRs, o2::aod::D0ChargedMCDetectorLevelSPs, o2::aod::D0ChargedMCDetectorLevelPRs, o2::aod::D0CMCPJetSSs, o2::aod::D0CMCPJetRs, o2::aod::D0ChargedMCParticleLevelSPs, o2::aod::D0ChargedMCParticleLevelPRs, o2::aod::D0CEWSJetSSs, o2::aod::D0ChargedEventWiseSubtractedSPs, o2::aod::D0ChargedEventWiseSubtractedPRs, o2::aod::JTrackD0Subs>;

extern template struct JetSubstructureHFTask<o2::soa::Join<o2::aod::DplusChargedJets, o2::aod::DplusChargedJetConstituents>, o2::soa::Join<o2::aod::DplusChargedMCDetectorLevelJets, o2::aod::DplusChargedMCDetectorLevelJetConstituents>, o2::soa::Join<o2::aod::DplusChargedMCParticleLevelJets, o2::aod::DplusChargedMCParticleLevelJetConstituents>, o2::soa::Join<o2::aod::DplusChargedEventWiseSubtractedJets, o2::aod::DplusChargedEventWiseSubtractedJetConstituents>, o2::aod::CandidatesDplusData, o2::aod::CandidatesDplusMCP, o2::aod::DplusCJetSSs, o2::aod::DplusCJetRs, o2::aod::DplusChargedSPs, o2::aod::DplusChargedPRs, o2::aod::DplusCMCDJetSSs, o2::aod::DplusCMCDJetRs, o2::aod::DplusChargedMCDetectorLevelSPs, o2::aod::DplusChargedMCDetectorLevelPRs, o2::aod::DplusCMCPJetSSs, o2::aod::DplusCMCPJetRs, o2::aod::DplusChargedMCParticleLevelSPs, o2::aod::DplusChargedMCParticleLevelPRs, o2::aod::DplusCEWSJetSSs, o2::aod::DplusChargedEventWiseSubtractedSPs, o2::aod::DplusChargedEventWiseSubtractedPRs, o2::aod::JTrackDplusSubs>;

extern template struct JetSubstructureHFTask<o2::soa::Join<o2::aod::DsChargedJets, o2::aod::DsChargedJetConstituents>, o2::soa::Join<o2::aod::DsChargedMCDetectorLevelJets, o2::aod::DsChargedMCDetectorLevelJetConstituents>, o2::soa::Join<o2::aod::DsChargedMCParticleLevelJets, o2::aod::DsChargedMCParticleLevelJetConstituents>, o2::soa::Join<o2::aod::DsChargedEventWiseSubtractedJets, o2::aod::DsChargedEventWiseSubtractedJetConstituents>, o2::aod::CandidatesDsData, o2::aod::CandidatesDsMCP, o2::aod::DsCJetSSs, o2::aod::DsCJetRs, o2::aod::DsChargedSPs, o2::aod::DsChargedPRs, o2::aod::DsCMCDJetSSs, o2::aod::DsCMCDJetRs, o2::aod::DsChargedMCDetectorLevelSPs, o2::aod::DsChargedMCDetectorLevelPRs, o2::aod::DsCMCPJetSSs, o2::aod::DsCMCPJetRs, o2::aod::DsChargedMCParticleLevelSPs, o2::aod::DsChargedMCParticleLevelPRs, o2::aod::DsCEWSJetSSs, o2::aod::DsChargedEventWiseSubtractedSPs, o2::aod::DsChargedEventWiseSubtractedPRs, o2::aod::JTrackDsSubs>;

extern template struct JetSubstructureHFTask<o2::soa::Join<o2::aod::DstarChargedJets, o2::aod::DstarChargedJetConstituents>, o2::soa::Join<o2::aod::DstarChargedMCDetectorLevelJets, o2::aod::DstarChargedMCDetectorLevelJetConstituents>, o2::soa::Join<o2::aod::DstarChargedMCParticleLevelJets, o2::aod::DstarChargedMCParticleLevelJetConstituents>, o2::soa::Join<o2::aod::DstarChargedEventWiseSubtractedJets, o2::aod::DstarChargedEventWiseSubtractedJetConstituents>, o2::aod::CandidatesDstarData, o2::aod::CandidatesDstarMCP, o2::aod::DstarCJetSSs, o2::aod::DstarCJetRs, o2::aod::DstarChargedSPs, o2::aod::DstarChargedPRs, o2::aod::DstarCMCDJetSSs, o2::aod::DstarCMCDJetRs, o2::aod::DstarChargedMCDetectorLevelSPs, o2::aod::DstarChargedMCDetectorLevelPRs, o2::aod::DstarCMCPJetSSs, o2::aod::DstarCMCPJetRs, o2::aod::DstarChargedMCParticleLevelSPs, o2::aod::DstarChargedMCParticleLevelPRs, o2::aod::DstarCEWSJetSSs, o2::aod::DstarChargedEventWiseSubtractedSPs, o2::aod::DstarChargedEventWiseSubtractedPRs, o2::aod::JTrackDstarSubs>;

extern template struct JetSubstructureHFTask<o2::soa::Join<o2::aod::LcChargedJets, o2::aod::LcChargedJetConstituents>, o2::soa::Join<o2::aod::LcChargedMCDetectorLevelJets, o2::aod::LcChargedMCDetectorLevelJetConstituents>, o2::soa::Join<o2::aod::LcChargedMCParticleLevelJets, o2::aod::LcChargedMCParticleLevelJetConstituents>, o2::soa::Join<o2::aod::LcChargedEventWiseSubtractedJets, o2::aod::LcChargedEventWiseSubtractedJetConstituents>, o2::aod::CandidatesLcData, o2::aod::CandidatesLcMCP, o2::aod::LcCJetSSs, o2::aod::LcCJetRs, o2::aod::LcChargedSPs, o2::aod::LcChargedPRs, o2::aod::LcCMCDJetSSs, o2::aod::LcCMCDJetRs, o2::aod::LcChargedMCDetectorLevelSPs, o2::aod::LcChargedMCDetectorLevelPRs, o2::aod::LcCMCPJetSSs, o2::aod::LcCMCPJetRs, o2::aod::LcChargedMCParticleLevelSPs, o2::aod::LcChargedMCParticleLevelPRs, o2::aod::LcCEWSJetSSs, o2::aod::LcChargedEventWiseSubtractedSPs, o2::aod::LcChargedEventWiseSubtractedPRs, o2::aod::JTrackLcSubs>;

extern template struct JetSubstructureHFTask<o2::soa::Join<o2::aod::B0ChargedJets, o2::aod::B0ChargedJetConstituents>, o2::soa::Join<o2::aod::B0ChargedMCDetectorLevelJets, o2::aod::B0ChargedMCDetectorLevelJetConstituents>, o2::soa::Join<o2::aod::B0ChargedMCParticleLevelJets, o2::aod::B0ChargedMCParticleLevelJetConstituents>, o2::soa::Join<o2::aod::B0ChargedEventWiseSubtractedJets, o2::aod::B0ChargedEventWiseSubtractedJetConstituents>, o2::aod::CandidatesB0Data, o2::aod::CandidatesB0MCP, o2::aod::B0CJetSSs, o2::aod::B0CJetRs, o2::aod::B0ChargedSPs, o2::aod::B0ChargedPRs, o2::aod::B0CMCDJetSSs, o2::aod::B0CMCDJetRs, o2::aod::B0ChargedMCDetectorLevelSPs, o2::aod::B0ChargedMCDetectorLevelPRs, o2::aod::B0CMCPJetSSs, o2::aod::B0CMCPJetRs, o2::aod::B0ChargedMCParticleLevelSPs, o2::aod::B0ChargedMCParticleLevelPRs, o2::aod::B0CEWSJetSSs, o2::aod::B0ChargedEventWiseSubtractedSPs, o2::aod::B0ChargedEventWiseSubtractedPRs, o2::aod::JTrackB0Subs>;

extern template struct JetSubstructureHFTask<o2::soa::Join<o2::aod::BplusChargedJets, o2::aod::BplusChargedJetConstituents>, o2::soa::Join<o2::aod::BplusChargedMCDetectorLevelJets, o2::aod::BplusChargedMCDetectorLevelJetConstituents>, o2::soa::Join<o2::aod::BplusChargedMCParticleLevelJets, o2::aod::BplusChargedMCParticleLevelJetConstituents>, o2::soa::Join<o2::aod::BplusChargedEventWiseSubtractedJets, o2::aod::BplusChargedEventWiseSubtractedJetConstituents>, o2::aod::CandidatesBplusData, o2::aod::CandidatesBplusMCP, o2::aod::BplusCJetSSs, o2::aod::BplusCJetRs, o2::aod::BplusChargedSPs, o2::aod::BplusChargedPRs, o2::aod::BplusCMCDJetSSs, o2::aod::BplusCMCDJetRs, o2::aod::BplusChargedMCDetectorLevelSPs, o2::aod::BplusChargedMCDetectorLevelPRs, o2::aod::BplusCMCPJetSSs, o2::aod::BplusCMCPJetRs, o2::aod::BplusChargedMCParticleLevelSPs, o2::aod::BplusChargedMCParticleLevelPRs, o2::aod::BplusCEWSJetSSs, o2::aod::BplusChargedEventWiseSubtractedSPs, o2::aod::BplusChargedEventWiseSubtractedPRs, o2::aod::JTrackBplusSubs>;

extern template struct JetSubstructureHFTask<o2::soa::Join<o2::aod::XicToXiPiPiChargedJets, o2::aod::XicToXiPiPiChargedJetConstituents>, o2::soa::Join<o2::aod::XicToXiPiPiChargedMCDetectorLevelJets, o2::aod::XicToXiPiPiChargedMCDetectorLevelJetConstituents>, o2::soa::Join<o2::aod::XicToXiPiPiChargedMCParticleLevelJets, o2::aod::XicToXiPiPiChargedMCParticleLevelJetConstituents>, o2::soa::Join<o2::aod::XicToXiPiPiChargedEventWiseSubtractedJets, o2::aod::XicToXiPiPiChargedEventWiseSubtractedJetConstituents>, o2::aod::CandidatesXicToXiPiPiData, o2::aod::CandidatesXicToXiPiPiMCP, o2::aod::XicToXiPiPiCJetSSs, o2::aod::XicToXiPiPiCJetRs, o2::aod::XicToXiPiPiChargedSPs, o2::aod::XicToXiPiPiChargedPRs, o2::aod::XicToXiPiPiCMCDJetSSs, o2::aod::XicToXiPiPiCMCDJetRs, o2::aod::XicToXiPiPiChargedMCDetectorLevelSPs, o2::aod::XicToXiPiPiChargedMCDetectorLevelPRs, o2::aod::XicToXiPiPiCMCPJetSSs, o2::aod::XicToXiPiPiCMCPJetRs, o2::aod::XicToXiPiPiChargedMCParticleLevelSPs, o2::aod::XicToXiPiPiChargedMCParticleLevelPRs, o2::aod::XicToXiPiPiCEWSJetSSs, o2::aod::XicToXiPiPiChargedEventWiseSubtractedSPs, o2::aod::XicToXiPiPiChargedEventWiseSubtractedPRs, o2::aod::JTrackXicToXiPiPiSubs>;

extern template struct JetSubstructureHFTask<o2::soa::Join<o2::aod::DielectronChargedJets, o2::aod::DielectronChargedJetConstituents>, o2::soa::Join<o2::aod::DielectronChargedMCDetectorLevelJets, o2::aod::DielectronChargedMCDetectorLevelJetConstituents>, o2::soa::Join<o2::aod::DielectronChargedMCParticleLevelJets, o2::aod::DielectronChargedMCParticleLevelJetConstituents>, o2::soa::Join<o2::aod::DielectronChargedEventWiseSubtractedJets, o2::aod::DielectronChargedEventWiseSubtractedJetConstituents>, o2::aod::CandidatesDielectronData, o2::aod::CandidatesDielectronMCP, o2::aod::DielectronCJetSSs, o2::aod::DielectronCJetRs, o2::aod::DielectronChargedSPs, o2::aod::DielectronChargedPRs, o2::aod::DielectronCMCDJetSSs, o2::aod::DielectronCMCDJetRs, o2::aod::DielectronChargedMCDetectorLevelSPs, o2::aod::DielectronChargedMCDetectorLevelPRs, o2::aod::DielectronCMCPJetSSs, o2::aod::DielectronCMCPJetRs, o2::aod::DielectronChargedMCParticleLevelSPs, o2::aod::DielectronChargedMCParticleLevelPRs, o2::aod::DielectronCEWSJetSSs, o2::aod::DielectronChargedEventWiseSubtractedSPs, o2::aod::DielectronChargedEventWiseSubtractedPRs, o2::aod::JTrackDielectronSubs>;

#endif // PWGJE_TASKS_JETSUBSTRUCTUREHF_H_
