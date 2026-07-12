#pragma once

#include <k4FWCore/Transformer.h>

#include <edm4hep/ReconstructedParticleCollection.h>

class pandoraReading :
    public k4FWCore::Transformer<
        edm4hep::ReconstructedParticleCollection(
            const edm4hep::ReconstructedParticleCollection&)> {

public:

    pandoraReading(const std::string& name,
                    ISvcLocator* svc);

    edm4hep::ReconstructedParticleCollection operator()(
        const edm4hep::ReconstructedParticleCollection&) const override;
};
