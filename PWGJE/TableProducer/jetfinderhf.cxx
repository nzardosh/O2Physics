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

// jet finder task
//
// Authors: Nima Zardoshti, Jochen Klein

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec hfjetMode = {
    "hfjetMode",
    VariantType::String,
    "",
    {"HF jet finder mode."},
  };
  workflowOptions.push_back(hfjetMode);
}

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename JetTable, typename ConstituentTable>
struct JetFinderHFTask {
  Produces<JetTable> jetsTable;
  Produces<ConstituentTable> constituentsTable;
  OutputObj<TH1F> hJetPt{"h_jet_pt"};
  OutputObj<TH1F> hJetPhi{"h_jet_phi"};
  OutputObj<TH1F> hJetEta{"h_jet_eta"};
  OutputObj<TH1F> hJetNTracks{"h_jet_ntracks"};
  OutputObj<TH1F> hD0Pt{"h_D0_pt"};

  Service<O2DatabasePDG> pdg;
  TrackSelection globalTracks;

  std::vector<fastjet::PseudoJet> jets;
  std::vector<fastjet::PseudoJet> inputParticles;
  JetFinder jetFinder;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.8, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.8, "maximum track eta"};

  // HF candidate level configurables
  Configurable<float> candPtMin{"candPtMin", 0.0, "minimum candidate pT"};
  Configurable<float> candPtMax{"candPtMax", 100.0, "maximum candidate pT"};
  Configurable<float> candYMin{"candYMin", -0.8, "minimum candidate eta"};
  Configurable<float> candYMax{"candYMax", 0.8, "maximum candidate eta"};
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};

  // jet level configurables
  Configurable<std::vector<double>> jetR{"jetR", {0.4}, "jet resolution parameters"};
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<int> jetType{"jetType", 1, "Type of stored jets. 0 = full, 1 = charged, 2 = neutral"};
  Configurable<int> jetAlgorithm{"jetAlgorithm", 2, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
  Configurable<int> jetRecombScheme{"jetRecombScheme", 0, "jet recombination scheme. 0 = E-scheme, 1 = pT-scheme, 2 = pT2-scheme"};
  Configurable<float> jetGhostArea{"jetGhostArea", 0.005, "jet ghost area"};

  void init(InitContext const&)
  {
    // set up global tracks and adjust as necessary
    globalTracks = getGlobalTrackSelection();
    globalTracks.SetEtaRange(trackEtaMin, trackEtaMax);

    hJetPt.setObject(new TH1F("h_jet_pt", "jet p_{T};p_{T} (GeV/#it{c})",
                              100, 0., 100.));
    hJetPhi.setObject(new TH1F("h_jet_phi", "jet #phi; #phi",
                               140, -7.0, 7.0));
    hJetEta.setObject(new TH1F("h_jet_eta", "jet #eta; #eta",
                               30, -1.5, 1.5));
    hJetNTracks.setObject(new TH1F("h_jet_ntracks", "jet N tracks ; N tracks",
                                   150, -0.5, 99.5));
    hD0Pt.setObject(new TH1F("h_D0_pt", "jet p_{T,D};p_{T,D} (GeV/#it{c})",
                             100, 0., 100.));

    jetFinder.etaMin = trackEtaMin;
    jetFinder.etaMax = trackEtaMax;
    jetFinder.jetPtMin = jetPtMin;
    jetFinder.jetPtMax = jetPtMax;
    jetFinder.algorithm = static_cast<fastjet::JetAlgorithm>(static_cast<int>(jetAlgorithm));
    jetFinder.recombScheme = static_cast<fastjet::RecombinationScheme>(static_cast<int>(jetRecombScheme));
    jetFinder.ghostArea = jetGhostArea;
  }

  //need enum as configurable
  enum pdgCode { pdgD0 = 421 };

  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter trackCuts = (aod::track::pt >= trackPtMin && aod::track::pt < trackPtMax && aod::track::eta > trackEtaMin && aod::track::eta < trackEtaMax);
  Filter partCuts = (aod::mcparticle::pt >= trackPtMin && aod::mcparticle::pt < trackPtMax);
  Filter candCuts = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar);

  template <typename T>
  void jetFinding(T const& collision)
  {
    auto jetRValues = static_cast<std::vector<double>>(jetR);
    auto candidatepT = 0.0;
    for (auto R : jetRValues) {
      jetFinder.jetR = R;
      jets.clear();
      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));
      for (const auto& jet : jets) {
        bool isHFJet = false;
        // if (jet.eta() < trackEtaMin + jetR || jet.eta() > trackEtaMax - jetR) { //should be done in jetFinder already
        //  continue;
        // }
        // if (jet.perp() < jetPtMin || jet.perp() >= jetPtMax) { //should be done in jetFinder already
        //  continue;
        // }
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidateHF)) {
            isHFJet = true;
            break;
          }
        }
        std::vector<int> trackconst;
        std::vector<int> candconst;

        if (isHFJet) {
          jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                    jet.E(), jet.m(), jet.area(), std::round(R * 100.0) / 100.0);
          // const auto& constituents = sorted_by_pt(jet.constituents());
          for (const auto& constituent : sorted_by_pt(jet.constituents())) {
            // need to add seperate thing for constituent subtraction

            if (constituent.user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::track)) {
              trackconst.push_back(constituent.user_info<FastJetUtilities::fastjet_user_info>().getIndex());
            }
            if (constituent.user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidateHF)) {
              candconst.push_back(constituent.user_info<FastJetUtilities::fastjet_user_info>().getIndex());
              candidatepT = constituent.pt();
            }
          }
          // candconst.push_back(candidate.globalIndex()); // is this grouped per collision too?
          constituentsTable(jetsTable.lastIndex(), trackconst, std::vector<int>(), candconst);
          hJetPt->Fill(jet.pt());
          hJetPhi->Fill(jet.phi());
          hJetEta->Fill(jet.rap());
          hJetNTracks->Fill(jet.constituents().size());
          hD0Pt->Fill(candidatepT);
          break;
        }
      }
    }
  }

  using JetTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                   JetTracks const& tracks,
                   soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> const& candidates)
  {
    // TODO: retrieve pion mass from somewhere
    if (!collision.sel8())
      return;

    //this loop should be made more efficient
    for (auto& candidate : candidates) {
      if (yD0(candidate) < candYMin || yD0(candidate) > candYMax) {
        continue;
      }
      if (candidate.pt() < candPtMin || candidate.pt() >= candPtMax) {
        continue;
      }
      inputParticles.clear();
      for (auto& track : tracks) {
        if (!globalTracks.IsSelected(track)) {
          continue;
        }
        if (candidate.prong0_as<JetTracks>().globalIndex() == track.globalIndex() || candidate.prong1_as<JetTracks>().globalIndex() == track.globalIndex()) {
          continue;
        }
        FastJetUtilities::fillTracks(track, inputParticles, track.globalIndex(), static_cast<int>(JetConstituentStatus::track));
      }
      FastJetUtilities::fillTracks(candidate, inputParticles, candidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(pdgD0));
      jetFinding(collision);
      /*
      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));
      for (const auto& jet : jets) {
        isHFJet = false;
        if (jet.eta() < trackEtaMin + jetR || jet.eta() > trackEtaMax - jetR) { //should be done in jetFinder already
          continue;
        }
        if (jet.perp() < jetPtMin || jet.perp() >= jetPtMax) { //should be done in jetFinder already
          continue;
        }
        std::vector<int> trackconst;
        std::vector<int> candconst;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_index() == -1) {
            isHFJet = true;
            break;
          }
        }

        if (isHFJet) {
          jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                    jet.E(), jet.m(), jet.area(), jetFinder.jetR);
          const auto& constituents = sorted_by_pt(jet.constituents());
          for (const auto& constituent : constituents) {
            if (constituent.user_index() != -1) {
              // auto track = tracks.rawIteratorAt(constituent.user_index() - tracks.offset());
              trackconst.push_back(constituent.user_index());
            }
          }
          candconst.push_back(candidate.globalIndex()); // is this grouped per collision too?
          constituentsTable(jetsTable.lastIndex(), trackconst, std::vector<int>(), candconst);
          hJetPt->Fill(jet.pt());
          hJetPhi->Fill(jet.phi());
          hJetEta->Fill(jet.rap());
          hJetNTracks->Fill(jet.constituents().size());
          hD0Pt->Fill(candidate.pt());
          break;
        }
      }
      */
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processData, "HF jet finding on data", true);

  void processMCD(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                  JetTracks const& tracks,
                  soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> const& candidates)
  {
    // TODO: retrieve pion mass from somewhere
    // this loop should be made more efficient
    // TODO: should probably refine the candidate selection
    for (auto& candidate : candidates) {
      if (yD0(candidate) < candYMin || yD0(candidate) > candYMax) {
        continue;
      }
      if (candidate.pt() < candPtMin || candidate.pt() >= candPtMax) {
        continue;
      }
      // are the next two ifs needed?
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (!(std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      inputParticles.clear();
      for (auto& track : tracks) {
        if (!globalTracks.IsSelected(track)) {
          continue;
        }
        if (candidate.prong0_as<JetTracks>().globalIndex() == track.globalIndex() || candidate.prong1_as<JetTracks>().globalIndex() == track.globalIndex()) {
          continue;
        }
        FastJetUtilities::fillTracks(track, inputParticles, track.globalIndex(), static_cast<int>(JetConstituentStatus::track));
      }
      FastJetUtilities::fillTracks(candidate, inputParticles, candidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(pdgD0));
      jetFinding(collision);
      /*
      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));
      for (const auto& jet : jets) {
        isHFJet = false;
        if (jet.eta() < trackEtaMin + jetR || jet.eta() > trackEtaMax - jetR) { // is this needed?
          continue;
        }
        if (jet.perp() < jetPtMin || jet.perp() >= jetPtMax) {
          continue;
        }
        std::vector<int> trackconst;
        std::vector<int> candconst;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_index() == -1) {
            isHFJet = true;
            break;
          }
        }
        if (isHFJet) {
          jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                    jet.E(), jet.m(), jet.area(), jetFinder.jetR);
          const auto& constituents = sorted_by_pt(jet.constituents());
          for (const auto& constituent : constituents) {
            if (constituent.user_index() != -1) {
              // auto track = tracks.rawIteratorAt(constituent.user_index() - tracks.offset());
              trackconst.push_back(constituent.user_index());
            }
          }
          candconst.push_back(candidate.globalIndex()); // check if its correct
          constituentsTable(jetsTable.lastIndex(), trackconst, std::vector<int>(), candconst);
          hJetPt->Fill(jet.pt());
          hJetPhi->Fill(jet.phi());
          hJetEta->Fill(jet.rap());
          hJetNTracks->Fill(jet.constituents().size());
          hD0Pt->Fill(candidate.pt());
          break;
        }
      }
      */
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processMCD, "HF jet finding on MC detector level", false);

  void processMCP(aod::McCollision const& collision,
                  soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>> const& particles)
  {
    LOG(debug) << "Per Event MCP";
    // TODO: retrieve pion mass from somewhere
    // TODO: probably should do this as a filter
    std::vector<soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>::iterator> candidates;
    for (auto const& particle : particles) {
      // TODO: generalise to any D0
      if (std::abs(particle.flagMcMatchGen()) & (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        auto particleY = RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        if (particleY < candYMin || particleY > candYMax) {
          continue;
        }
        if (particle.pt() < candPtMin || particle.pt() >= candPtMax) {
          continue;
        }
        candidates.push_back(particle);
      }
    }

    //this loop should be made more efficient
    for (auto& candidate : candidates) {
      inputParticles.clear();
      for (auto& particle : particles) {
        // exclude neutral particles
        // TODO: can we do this through the filter?
        if (particle.eta() < trackEtaMin || particle.eta() > trackEtaMax) {
          continue;
        }
        if (particle.getGenStatusCode() != 1) {
          continue;
        }
        auto pdgParticle = pdg->GetParticle(particle.pdgCode());
        auto pdgCharge = pdgParticle ? std::abs(pdgParticle->Charge()) : -1.0;
        if (jetType == static_cast<int>(JetType::charged) && pdgCharge < 3.0) {
          continue;
        }
        if (jetType == static_cast<int>(JetType::neutral) && pdgCharge != 0.0) {
          continue;
        }
        // TODO: check what mass to use?
        const auto daughters = candidate.daughtersIds();
        if (std::find(std::begin(daughters), std::end(daughters), particle.globalIndex()) != std::end(daughters)) {
          continue;
        }
        FastJetUtilities::fillTracks(particle, inputParticles, particle.globalIndex(), static_cast<int>(JetConstituentStatus::track), RecoDecay::getMassPDG(particle.pdgCode()));
      }
      FastJetUtilities::fillTracks(candidate, inputParticles, candidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candidate.pdgCode()));
      jetFinding(collision);
      /*
      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));

      for (const auto& jet : jets) {
        if (jet.eta() < trackEtaMin + jetR || jet.eta() > trackEtaMax - jetR) {
          continue;
        }
        if (jet.perp() < jetPtMin || jet.perp() >= jetPtMax) {
          continue;
        }
        isHFJet = false;
        std::vector<int> trackconst;
        std::vector<int> candconst;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_index() == -1) {
            isHFJet = true;
            break;
          }
        }
        if (isHFJet) {
          jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                    jet.E(), jet.m(), jet.area(), jetFinder.jetR);
          for (const auto& constituent : jet.constituents()) {
            if (constituent.user_index() == -1)
              continue;
            trackconst.push_back(constituent.user_index());
          }

          candconst.push_back(candidate.globalIndex());
          constituentsTable(jetsTable.lastIndex(), trackconst, std::vector<int>(), candconst);
          hJetPt->Fill(jet.pt());
          hJetPhi->Fill(jet.phi());
          hJetEta->Fill(jet.rap());
          hJetNTracks->Fill(jet.constituents().size());
          hD0Pt->Fill(candidate.pt());
          break;
        }
      }
      */
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processMCP, "HF jet finding on MC particle level", false);
};

using JetFinderHF = JetFinderHFTask<o2::aod::HFJets, o2::aod::HFJetConstituents>;
using MCParticleLevelJetFinderHF = JetFinderHFTask<o2::aod::MCParticleLevelHFJets, o2::aod::MCParticleLevelHFJetConstituents>;
using MCDetectorLevelJetFinderHF = JetFinderHFTask<o2::aod::MCDetectorLevelHFJets, o2::aod::MCDetectorLevelHFJetConstituents>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  auto hfjetMode = cfgc.options().get<std::string>("hfjetMode");

  if (hfjetMode.find("data") != std::string::npos || hfjetMode.empty())
    tasks.emplace_back(adaptAnalysisTask<JetFinderHF>(cfgc,
                                                      SetDefaultProcesses{{{"processData", true}}},
                                                      TaskName{"jet-finder-hf-data"}));

  if (hfjetMode.find("mcp") != std::string::npos || hfjetMode.empty())
    tasks.emplace_back(adaptAnalysisTask<MCParticleLevelJetFinderHF>(cfgc,
                                                                     SetDefaultProcesses{{{"processData", false},{"processMCP", true}}},
                                                                     TaskName{"jet-finder-hf-mcp"}));

  if (hfjetMode.find("mcd") != std::string::npos || hfjetMode.empty())
    tasks.emplace_back(adaptAnalysisTask<MCDetectorLevelJetFinderHF>(cfgc,
                                                                     SetDefaultProcesses{{{"processData", false},{"processMCD", true}}},
                                                                     TaskName{"jet-finder-hf-mcd"}));

  return WorkflowSpec{tasks};
}
