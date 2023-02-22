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
// Author: Jochen Klein, Nima Zardoshti, Raymond Ehlers

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/RecoDecay.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
/*
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec jetData = {
    "jetInputData",
    VariantType::String,
    "",
    {"Jet input data type. Options include Data, MCParticleLevel, MCDetectorLevel, and HybridIntermediate."},
  };
  workflowOptions.push_back(jetData);
  ConfigParamSpec jetType = {
    "jetType",
    VariantType::Int,
    1,
    {"Type of jet. Options include full = 0, charged = 1, neutral = 2."},
  };
  workflowOptions.push_back(jetType);
}
*/
/*
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec jetMode = {
    "jetMode",
    VariantType::String,
    "",
    {"jet finder mode."},
  };
  workflowOptions.push_back(jetMode);
}
*/

#include "Framework/runDataProcessing.h"

template <typename JetTable, typename ConstituentTable, typename ConstituentSubTable>
struct JetFinderTask {
  Produces<JetTable> jetsTable;
  Produces<ConstituentTable> constituentsTable;
  Produces<ConstituentSubTable> constituentsSubTable;
  OutputObj<TH2F> hJetPt{"h_jet_pt"};
  OutputObj<TH2F> hJetPhi{"h_jet_phi"};
  OutputObj<TH2F> hJetEta{"h_jet_eta"};
  OutputObj<TH2F> hJetN{"h_jet_n"};

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

  // cluster level configurables
  Configurable<int> mClusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};

  // jet level configurables
  Configurable<std::vector<double>> jetR{"jetR", {0.4}, "jet resolution parameters"};
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<int> jetTypeParticleLevel{"jetTypeParticleLevel", 1, "Type of stored jets at MC particle leevel. 0 = full, 1 = charged, 2 = neutral"};
  Configurable<int> jetAlgorithm{"jetAlgorithm", 2, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
  Configurable<int> jetRecombScheme{"jetRecombScheme", 0, "jet recombination scheme. 0 = E-scheme, 1 = pT-scheme, 2 = pT2-scheme"};
  Configurable<float> jetGhostArea{"jetGhostArea", 0.005, "jet ghost area"};
  Configurable<bool> DoRhoAreaSub{"DoRhoAreaSub", false, "do rho area subtraction"};
  Configurable<bool> DoConstSub{"DoConstSub", false, "do constituent subtraction"};

  void init(InitContext const&)
  {

    globalTracks = getGlobalTrackSelection();
    globalTracks.SetEtaRange(trackEtaMin, trackEtaMax);

    hJetPt.setObject(new TH2F("h_jet_pt", "jet p_{T};p_{T} (GeV/#it{c})",
                              100, 0., 100., 10, 0.05, 1.05));
    hJetPhi.setObject(new TH2F("h_jet_phi", "jet #phi;#phi",
                               80, -1., 7., 10, 0.05, 1.05));
    hJetEta.setObject(new TH2F("h_jet_eta", "jet #eta;#eta",
                               70, -0.7, 0.7, 10, 0.05, 1.05));
    hJetN.setObject(new TH2F("h_jet_n", "jet n;n constituents",
                             30, 0., 30., 10, 0.05, 1.05));
    if (DoRhoAreaSub) {
      jetFinder.setBkgSubMode(JetFinder::BkgSubMode::rhoAreaSub);
    }
    if (DoConstSub) {
      jetFinder.setBkgSubMode(JetFinder::BkgSubMode::constSub);
    }

    jetFinder.etaMin = trackEtaMin;
    jetFinder.etaMax = trackEtaMax;
    jetFinder.jetPtMin = jetPtMin;
    jetFinder.jetPtMax = jetPtMax;
    jetFinder.algorithm = static_cast<fastjet::JetAlgorithm>(static_cast<int>(jetAlgorithm));
    jetFinder.recombScheme = static_cast<fastjet::RecombinationScheme>(static_cast<int>(jetRecombScheme));
    jetFinder.ghostArea = jetGhostArea;
  }


  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter trackCuts = (aod::track::pt >= trackPtMin && aod::track::pt < trackPtMax && aod::track::eta > trackEtaMin && aod::track::eta < trackEtaMax); // do we need eta cut both here and in globalselection?
  Filter partCuts = (aod::mcparticle::pt >= trackPtMin && aod::mcparticle::pt < trackPtMax);
  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition);

  using JetTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;
  using JetClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;

  template <typename T>
  void processTracks(T const& tracks)
  {
    for (auto& track : tracks) {
      if (!globalTracks.IsSelected(track)) {
        continue;
      }
      FastJetUtilities::fillTracks(track, inputParticles, track.globalIndex());
    }
  }

  template <typename T>
  void processClusters(T const& clusters)
  {
    for (auto& cluster : *clusters) {
      // add cluster selections
      FastJetUtilities::fillClusters(cluster, inputParticles, cluster.globalIndex());
    }
  }

  template <typename T>
  void jetFinding(T const& collision)
  {
    LOG(debug) << "Process Implementation";
    // NOTE: Can't just iterate directly - we have to cast first
    auto jetRValues = static_cast<std::vector<double>>(jetR);
    for (auto R : jetRValues) {
      // Update jet finder R and find jets
      jetFinder.jetR = R;
      jets.clear();
      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));

      for (const auto& jet : jets) {
        jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                  jet.E(), jet.m(), jet.area(), std::round(R * 100.0) / 100.0);

        std::vector<int> trackconst; // FIXME:initialise as empty
        std::vector<int> clusterconst;

        for (const auto& constituent : sorted_by_pt(jet.constituents())) { // event or jetwise
          if (DoConstSub) {                                                // FIXME: needs to be addressed in Haadi's PR
            constituentsSubTable(jetsTable.lastIndex(), constituent.pt(), constituent.eta(), constituent.phi(),
                                 constituent.E(), constituent.m(), constituent.user_index());
          }

          if (constituent.user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::track)) {
            trackconst.push_back(constituent.user_info<FastJetUtilities::fastjet_user_info>().getIndex());
          }
          if (constituent.user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::cluster)) {
            clusterconst.push_back(constituent.user_info<FastJetUtilities::fastjet_user_info>().getIndex());
          }
          constituentsTable(jetsTable.lastIndex(), trackconst, clusterconst, std::vector<int>());
        }

        hJetPt->Fill(jet.pt(), R);
        hJetPhi->Fill(jet.phi(), R);
        hJetEta->Fill(jet.eta(), R);
        hJetN->Fill(jet.constituents().size(), R);
      }
    }
  }

  void processParticleLevel(aod::McCollision const& collision, aod::McParticles const& particles)
  {
    // Initialziation and event selection
    // TODO: MC event selection?
    inputParticles.clear();
    for (auto& particle : particles) {
      if (particle.getGenStatusCode() != 1) {
        continue;
      }
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      auto pdgCharge = pdgParticle ? std::abs(pdgParticle->Charge()) : -1.0;
      if (jetTypeParticleLevel == static_cast<int>(JetType::charged) && pdgCharge < 3.0) { // pdg charge is in increments of 1/3
        continue;
      }
      if (jetTypeParticleLevel == static_cast<int>(JetType::neutral) && pdgCharge != 0.0) {
        continue;
      }
      FastJetUtilities::fillTracks(particle, inputParticles, particle.globalIndex(), static_cast<int>(JetConstituentStatus::track), RecoDecay::getMassPDG(particle.pdgCode()));
    }
    jetFinding(collision);
  }

  PROCESS_SWITCH(JetFinderTask, processParticleLevel, "Particle level jet finding", false);

  void processDataCharged(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                          JetTracks const& tracks)
  {
    LOG(debug) << "Process data charged!";
    processTracks(tracks);
    jetFinding(collision);
  }

  PROCESS_SWITCH(JetFinderTask, processDataCharged, "Data jet finding for charged jets", true);

  void processDataNeutral(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                          JetClusters const& clusters)
  {
    LOG(debug) << "Process data neutral!";
    processClusters(&clusters);
    jetFinding(collision);
  }
  PROCESS_SWITCH(JetFinderTask, processDataNeutral, "Data jet finding for neutral jets", false);

  void processDataFull(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                       JetTracks const& tracks,
                       JetClusters const& clusters)
  {
    LOG(debug) << "Process data full!";
    processTracks(tracks);
    processClusters(&clusters);
    jetFinding(collision);
  }

  PROCESS_SWITCH(JetFinderTask, processDataFull, "Data jet finding for full and neutral jets", false);
};

using JetFinderDataCharged = JetFinderTask<o2::aod::Jets, o2::aod::JetConstituents, o2::aod::JetConstituentsSub>;
using JetFinderDataNeutral = JetFinderTask<o2::aod::Jets, o2::aod::JetConstituents, o2::aod::JetConstituentsSub>;
using JetFinderDataFull = JetFinderTask<o2::aod::Jets, o2::aod::JetConstituents, o2::aod::JetConstituentsSub>;
using JetFinderMCParticleLevel = JetFinderTask<o2::aod::MCParticleLevelJets, o2::aod::MCParticleLevelJetConstituents, o2::aod::MCParticleLevelJetConstituentsSub>;
using JetFinderMCDetectorLevelCharged = JetFinderTask<o2::aod::MCDetectorLevelJets, o2::aod::MCDetectorLevelJetConstituents, o2::aod::MCDetectorLevelJetConstituentsSub>;
using JetFinderMCDetectorLevelNeutral = JetFinderTask<o2::aod::MCDetectorLevelJets, o2::aod::MCDetectorLevelJetConstituents, o2::aod::MCDetectorLevelJetConstituentsSub>;
using JetFinderMCDetectorLevelFull = JetFinderTask<o2::aod::MCDetectorLevelJets, o2::aod::MCDetectorLevelJetConstituents, o2::aod::MCDetectorLevelJetConstituentsSub>;
using JetFinderHybridIntermediate = JetFinderTask<o2::aod::HybridIntermediateJets, o2::aod::HybridIntermediateJetConstituents, o2::aod::HybridIntermediateJetConstituentsSub>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  std::vector<o2::framework::DataProcessorSpec> tasks;

  // auto jetMode = cfgc.options().get<std::string>("hfjetMode");

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderDataCharged>(cfgc,
                                            SetDefaultProcesses{{{"processDataCharged", true}}}, TaskName{"jet-finder-data"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderDataNeutral>(cfgc,
                                            SetDefaultProcesses{{{"processDataCharged", false}, {"processDataNeutral", true}}}, TaskName{"jet-finder-data-neutral"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderDataFull>(cfgc,
                                         SetDefaultProcesses{{{"processDataCharged", false}, {"processDataFull", true}}}, TaskName{"jet-finder-data-full"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderMCDetectorLevelCharged>(cfgc,
                                                       SetDefaultProcesses{{{"processDataCharged", true}}}, TaskName{"jet-finder-MC-detector-level"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderMCDetectorLevelCharged>(cfgc,
                                                       SetDefaultProcesses{{{"processDataCharged", false}, {"processDataNeutral", true}}}, TaskName{"jet-finder-MC-detector-level-neutral"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderMCDetectorLevelCharged>(cfgc,
                                                       SetDefaultProcesses{{{"processDataCharged", false}, {"processDataFull", true}}}, TaskName{"jet-finder-MC-detector-level-full"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderHybridIntermediate>(cfgc,
                                                   SetDefaultProcesses{{{"processDataFull", true}}}, TaskName{"jet-finder-hybrid-intermedaite-full"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderMCParticleLevel>(cfgc,
                                                SetDefaultProcesses{{{"processParticleLevel", true}, {"processDataCharged", false}}}, TaskName{"jet-finder-MC"}));

  return WorkflowSpec{tasks};

  /*


    // Determine the jet input data types.
    // NOTE: It's possible to create multiple jet finders by passing multiple jet input data types in the string.
    //       Each type is separated by a comma.
    // FIXME: String validation and normalization. There's too much room for error with the strings. (Is there a better way to do this?)
    auto jetInputData = cfgc.options().get<std::string>("jetInputData");
    const std::map<std::string, JetInputData_t> jetInputDataTypes = {
      {"Data", JetInputData_t::Data},
      {"MCParticleLevel", JetInputData_t::MCParticleLevel},
      {"MCDetectorLevel", JetInputData_t::MCDetectorLevel},
      {"HybridIntermediate", JetInputData_t::HybridIntermediate},
      {"", JetInputData_t::Data}, // Default to data
    };
    // Tokenize using stringstream
    std::vector<std::string> jetInputDataOptions;
    std::stringstream ss;
    ss << jetInputData;
    while (ss.good()) {
      std::string substring;
      getline(ss, substring, ',');
      jetInputDataOptions.push_back(substring);
    }
    // Jet type
    auto jetType = static_cast<JetType_t>(cfgc.options().get<int>("jetType"));
    std::vector<o2::framework::DataProcessorSpec> tasks;
    for (auto opt : jetInputDataOptions) {
      auto jetData = jetInputDataTypes.at(opt);
      switch (jetData) {
        case JetInputData_t::MCParticleLevel:
          // We don't need a switch on jet type here because the arguments for particle level jets are the same
          // for charged and full jets
          tasks.emplace_back(
            adaptAnalysisTask<JetFinderMCParticleLevel>(cfgc,
                                                        SetDefaultProcesses{{{"processParticleLevel", true}, {"processDataCharged", false}, {"processDataFull", false}}}, TaskName{"jet-finder-MC"}));
          break;
        case JetInputData_t::MCDetectorLevel:
          if (jetType == JetType_t::full || jetType == JetType_t::neutral) {
            tasks.emplace_back(
              adaptAnalysisTask<JetFinderMCDetectorLevel>(cfgc,
                                                          SetDefaultProcesses{{{"processParticleLevel", false}, {"processDataCharged", false}, {"processDataFull", true}}}, TaskName{"jet-finder-MC-detector-level"}));
          } else {
            tasks.emplace_back(
              adaptAnalysisTask<JetFinderMCDetectorLevel>(cfgc, TaskName{"jet-finder-MC-detector-level"}));
          }
          break;
        case JetInputData_t::HybridIntermediate:
          if (jetType == JetType_t::full || jetType == JetType_t::neutral) {
            tasks.emplace_back(
              adaptAnalysisTask<JetFinderHybridIntermediate>(cfgc,
                                                             SetDefaultProcesses{{{"processParticleLevel", false}, {"processDataCharged", false}, {"processDataFull", true}}}, TaskName{"jet-finder-hybrid-intermedaite"}));
          } else {
            tasks.emplace_back(
              adaptAnalysisTask<JetFinderHybridIntermediate>(cfgc, TaskName{"jet-finder-hybrid-intermedaite"}));
          }
          break;
        case JetInputData_t::Data: // Intentionally fall through to the default for data.
        default:
          if (jetType == JetType_t::full || jetType == JetType_t::neutral) {
            tasks.emplace_back(
              adaptAnalysisTask<JetFinderData>(cfgc,
                                               SetDefaultProcesses{{{"processParticleLevel", false}, {"processDataCharged", false}, {"processDataFull", true}}}, TaskName{"jet-finder-data"}));
          } else {
            tasks.emplace_back(
              adaptAnalysisTask<JetFinderData>(cfgc, TaskName{"jet-finder-data"}));
          }
          break;
      }
    }
    return WorkflowSpec{tasks};*/
}
