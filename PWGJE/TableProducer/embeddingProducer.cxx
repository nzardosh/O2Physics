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

// task to produce embedded tables
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/Vertex.h"
#include <CommonConstants/MathConstants.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
O2ORIGIN("EMB");
}

struct JetEmbeddingProducerTask {

  Produces<aod::StoredJBCs> storedJBCsTable;
  Produces<aod::StoredJCollisions> storedJCollisionsTable;
  Produces<aod::StoredJCollisionUPCs> storedJCollisionUPCsTable;
  Produces<aod::StoredJCollisionMcInfos> storedJCollisionMcInfosTable;
  Produces<aod::StoredJMcCollisionLbs> storedJMcCollisionsLabelTable;
  Produces<aod::StoredJTracks> storedJTracksTable;
  Produces<aod::StoredJMcTrackLbs> storedJMcTracksLabelTable;

  Produces<aod::StoredJMcCollisions> storedJMcCollisionsTable;
  Produces<aod::StoredJMcParticles> storedJMcParticlesTable;

  Preslice<aod::JetCollisions> perBCCollisions = aod::jcollision::bcId;
  Preslice<aod::JetTracks> perCollisionTracks = aod::jtrack::collisionId;
  Preslice<aod::JetParticles> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;
  Preslice<aod::JTracksFrom<o2::aod::Hash<"EMB"_h>>> perCollisionTargetTracks = aod::jtrack::collisionId;

  void init(InitContext const&)
  {
  }

  std::vector<uint32_t> collisionMapping;
  std::vector<int32_t> mcCollisionMapping;
  std::vector<int32_t> particleMapping;

  void processCollisionMaps(aod::JCollisionsFrom<o2::aod::Hash<"EMB"_h>>& targetCollisions, aod::JBCsFrom<o2::aod::Hash<"EMB"_h>> const& targetBCs, aod::JCollisions& signalCollisions, aod::JBCs const& signalBCs)
  {
    //targetCollisions.bindExternalIndices(&targetBCs);
    //signalCollisions.bindExternalIndices(&signalBCs);
    collisionMapping.clear();
    collisionMapping.resize(targetCollisions.size(), 4294967290); // magic number, but -1 will crash later

    for (auto const& targetCollision : targetCollisions) {
      auto globalBC = targetCollision.bc_as<aod::JBCsFrom<o2::aod::Hash<"EMB"_h>>>().globalBC();
      for (auto const& signalBC : signalBCs) {
        if (signalBC.globalBC() == globalBC) {
          const auto bcCollisions = signalCollisions.sliceBy(perBCCollisions, signalBC.globalIndex());
          for (auto const& signalCollision : bcCollisions) {
            collisionMapping[targetCollision.globalIndex()] = signalCollision.globalIndex();
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetEmbeddingProducerTask, processCollisionMaps, "produces collision table mapping", true);

  void processEmbeddingMCP(aod::JetMcCollisions const& mcCollisions, aod::JetParticles& particles)
  {
   // particles.bindExternalIndices(&mcCollisions);
    mcCollisionMapping.clear();
    mcCollisionMapping.resize(mcCollisions.size(), -1);
    particleMapping.clear();
    std::cout << "particles size: " << particles.size() << std::endl;
    particleMapping.resize(particles.size(), -1);
    int particleTableIndex = 0;
    for (auto const& mcCollision : mcCollisions) {
      storedJMcCollisionsTable(mcCollision.bcId(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.multFV0A(), mcCollision.multFT0A(), mcCollision.multFT0C(), mcCollision.centFT0M(), mcCollision.weight(), mcCollision.accepted(), mcCollision.attempted(), mcCollision.xsectGen(), mcCollision.xsectErr(), mcCollision.ptHard(), mcCollision.eventSel(), mcCollision.rct_raw(), mcCollision.getGeneratorId(), mcCollision.getSubGeneratorId(), mcCollision.getSourceId(), mcCollision.impactParameter(), mcCollision.eventPlaneAngle());
      mcCollisionMapping[mcCollision.globalIndex()] = storedJMcCollisionsTable.lastIndex();

      const auto particlesPerMcCollision = particles.sliceBy(ParticlesPerMcCollision, mcCollision.globalIndex());

      for (auto particle : particlesPerMcCollision) {
        particleMapping[particle.globalIndex()] = particleTableIndex;
        particleTableIndex++;
      }
      for (auto particle : particlesPerMcCollision) {

        std::vector<int32_t> mothersIds;
        int daughtersIds[2] = {-1, -1};
        if (particle.has_mothers()) {
          auto mothersIdTemps = particle.mothersIds();
          for (auto mothersIdTemp : mothersIdTemps) {
            mothersIds.push_back(particleMapping[mothersIdTemp]);
          }
        }
        if (particle.has_daughters()) {
          auto i = 0;
          for (auto daughterId : particle.daughtersIds()) {
            if (i > 1) {
              break;
            }
            daughtersIds[i] = particleMapping[daughterId];
            i++;
          }
        }

        storedJMcParticlesTable(mcCollisionMapping[mcCollision.globalIndex()], particle.pt(), particle.eta(), particle.phi(), particle.y(), particle.e(), particle.pdgCode(), particle.statusCode(), particle.flags(), mothersIds, daughtersIds);
      }
    }
  }
  PROCESS_SWITCH(JetEmbeddingProducerTask, processEmbeddingMCP, "produces mcp embedding table", true);

  // void processEmbedding(aod::JCollisionsFrom<o2::aod::Hash<"EMB"_h>>::iterator const& targetCollision, aod::JBCsFrom<o2::aod::Hash<"EMB"_h>> const&, aod::JetCollisionsMCD const& signalCollisions, aod::JTracksFrom<o2::aod::Hash<"EMB"_h>> const& targetTracks, aod::JetTracksMCD const& signalTracks)
  void processEmbedding(aod::JCollisionsFrom<o2::aod::Hash<"EMB"_h>>& targetCollisions, aod::JBCsFrom<o2::aod::Hash<"EMB"_h>> const& targetBCs, aod::JetCollisionsMCD const& signalCollisions, aod::JetTracksMCD const& signalTracks, aod::JTracksFrom<o2::aod::Hash<"EMB"_h>>& targetTracks)
  {
    //targetCollisions.bindExternalIndices(&targetBCs);
    //targetTracks.bindExternalIndices(&targetCollisions);
    std::cout << "here1" << std::endl;
    for (auto const& targetCollision : targetCollisions) {
      if (collisionMapping[targetCollision.globalIndex()] == 4294967290) {
        continue;
      }
      auto signalCollision = signalCollisions.iteratorAt(collisionMapping[targetCollision.globalIndex()]);
      std::cout << collisionMapping[targetCollision.globalIndex()] << "  " << signalCollisions.size() << std::endl;
      const auto signalCollisionTracks = signalTracks.sliceBy(perCollisionTracks, signalCollision.globalIndex());
      const auto bc = targetCollision.bc_as<aod::JBCsFrom<o2::aod::Hash<"EMB"_h>>>();
      storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.triggerMask(), bc.timestamp(), bc.alias_raw(), bc.selection_raw(), bc.rct_raw());
      storedJCollisionsTable(storedJBCsTable.lastIndex(), targetCollision.posX(), targetCollision.posY(), targetCollision.posZ(), targetCollision.collisionTime(), targetCollision.multFV0A(), targetCollision.multFV0C(), targetCollision.multFT0A(), targetCollision.multFT0C(), targetCollision.centFV0A(), targetCollision.centFV0M(), targetCollision.centFT0A(), targetCollision.centFT0C(), targetCollision.centFT0M(), targetCollision.centFT0CVariant1(), targetCollision.hadronicRate(), targetCollision.trackOccupancyInTimeRange(), targetCollision.alias_raw(), targetCollision.eventSel(), targetCollision.rct_raw(), targetCollision.triggerSel());
      storedJCollisionMcInfosTable(signalCollision.weight(), jetderiveddatautilities::JCollisionSubGeneratorId::none);
      if (signalCollision.has_mcCollision()) {
        storedJMcCollisionsLabelTable(mcCollisionMapping[signalCollision.mcCollisionId()]);
      } else {
        storedJMcCollisionsLabelTable(-1);
      }
      auto const& targetTracksSliced = targetTracks.sliceBy(perCollisionTargetTracks, targetCollision.globalIndex());
      for (auto const& targetTrack : targetTracksSliced) {
        storedJTracksTable(storedJCollisionsTable.lastIndex(), targetTrack.pt(), targetTrack.eta(), targetTrack.phi(), targetTrack.trackSel());
        storedJMcTracksLabelTable(-1);
      }
      for (auto const& signalTrack : signalCollisionTracks) {
        storedJTracksTable(storedJCollisionsTable.lastIndex(), signalTrack.pt(), signalTrack.eta(), signalTrack.phi(), signalTrack.trackSel());
        if (signalTrack.has_mcParticle()) {
          storedJMcTracksLabelTable(particleMapping[signalTrack.mcParticleId()]);
        } else {
          storedJMcTracksLabelTable(-1);
        }
      }
    }
  }
  PROCESS_SWITCH(JetEmbeddingProducerTask, processEmbedding, "produces embedding table", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetEmbeddingProducerTask>(cfgc, TaskName{"jet-embedding-producer"})};
}
