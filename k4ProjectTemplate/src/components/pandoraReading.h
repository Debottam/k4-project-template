#pragma once

#include <k4FWCore/Transformer.h>

#include <edm4hep/ReconstructedParticleCollection.h>

// marlin
#include <marlinutil/HelixClass_double.h>

// pandora
#include <Objects/Helix.h>

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Volumes.h"
#include "GaudiKernel/ServiceHandle.h"
#include "k4Interface/IGeoSvc.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/TrackCollection.h"
#include <podio/LinkCollection.h>

using MCTruthRecoLinkCollection =
    podio::LinkCollection<edm4hep::ReconstructedParticle, edm4hep::MCParticle>;


class pandoraReading :
    public k4FWCore::Transformer<
        edm4hep::ReconstructedParticleCollection(
            const edm4hep::ReconstructedParticleCollection&,
            const edm4hep::MCParticleCollection&,
            const MCTruthRecoLinkCollection&,
            const edm4hep::ClusterCollection&,
            const edm4hep::TrackCollection& )> {

public:

    pandoraReading(const std::string& name,
                    ISvcLocator* svc);

    edm4hep::ReconstructedParticleCollection operator()(
        const edm4hep::ReconstructedParticleCollection&,
        const edm4hep::MCParticleCollection&,
        const MCTruthRecoLinkCollection&,
        const edm4hep::ClusterCollection&,
        const edm4hep::TrackCollection&) const override;

private:
     ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Geometry service"};
};
