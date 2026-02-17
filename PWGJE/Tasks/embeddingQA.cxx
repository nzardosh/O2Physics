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

// jet tutorial task for hands on tutorial session (09/11/2023)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct EmbeddingQATask {
  HistogramRegistry registry{"registry",
                             {{"h_tracks_all_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_tracks_embedded_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_tracks_sub_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_particles_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jets_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jets_sub_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jets_mcd_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jets_mcp_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_detector_response", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_jet_response", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_jet_background_response", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}}}};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  std::vector<int> eventSelection;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  void processTracks(aod::JetCollision const&, aod::JetTracks const& tracks, aod::JetTracksSub const& subTracks)
  {
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        registry.fill(HIST("h_tracks_all_pt"), track.pt());
      }
      if (jetderiveddatautilities::selectTrack(track, trackSelection, true)) {
        registry.fill(HIST("h_tracks_embedded_pt"), track.pt());
      }
    }
    for (auto const& track : subTracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        registry.fill(HIST("h_tracks_sub_pt"), track.pt());
      }
    }
  }
  PROCESS_SWITCH(EmbeddingQATask, processTracks, "track QA", true);

  void processParticles(aod::JetMcCollision const&, aod::JetParticles const& particles)
  {
    for (auto const& particle : particles) {
      registry.fill(HIST("h_particles_pt"), particle.pt());
    }
  }
  PROCESS_SWITCH(EmbeddingQATask, processParticles, "track QA", true);

  void processJets(aod::ChargedJets const& jets, aod::ChargedEventWiseSubtractedJets const& subJets, aod::ChargedMCDetectorLevelJets const& mcdJets, aod::ChargedMCParticleLevelJets const& mcpJets)
  {
    for (auto const& jet : jets) {
      registry.fill(HIST("h_jets_pt"), jet.pt());
    }
    for (auto const& jet : subJets) {
      registry.fill(HIST("h_jets_sub_pt"), jet.pt());
    }
    for (auto const& jet : mcdJets) {
      registry.fill(HIST("h_jets_mcd_pt"), jet.pt());
    }
    for (auto const& jet : mcpJets) {
      registry.fill(HIST("h_jets_mcp_pt"), jet.pt());
    }
  }
  PROCESS_SWITCH(EmbeddingQATask, processJets, "jet QA", true);

  void processResponse(aod::ChargedJets const& jets, soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets> const& subJets, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedEventWiseSubtractedJets> const& mcdJets, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets> const& mcpJets)
  {

    for (auto const& mcpJet : mcpJets) {
      for (auto const& mcdJet : mcpJet.matchedJetGeo_as<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedEventWiseSubtractedJets>>()) {
        registry.fill(HIST("h_jet_detector_response"), mcpJet.pt(), mcdJet.pt());
        for (auto const& subJet : mcdJet.matchedJetGeo_as<soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets>>()) {
          registry.fill(HIST("h_jet_response"), mcpJet.pt(), subJet.pt());
          for (auto const& jet : subJet.matchedJetGeo_as<aod::ChargedJets>()) {
            registry.fill(HIST("h_jet_background_response"), subJet.pt(), jet.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(EmbeddingQATask, processResponse, "jet response", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EmbeddingQATask>(cfgc, TaskName{"embedding-qa"})}; }
