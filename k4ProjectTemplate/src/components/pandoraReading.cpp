#include "pandoraReading.h"

#include <GaudiKernel/MsgStream.h>

DECLARE_COMPONENT(pandoraReading)

pandoraReading::pandoraReading( const std::string& name,ISvcLocator* svc)

    : Transformer(name,svc,{
              KeyValue{"InputPFOs",
                       "PandoraPFOs"}
          },

          {
              KeyValue{"OutputPFOs",
                       "RecoveredPFOs"}
          })
{
std::cout << "pandoraReading constructor called" << std::endl;
}

edm4hep::ReconstructedParticleCollection 
pandoraReading::operator()(const edm4hep::ReconstructedParticleCollection& input) const {

    //throw std::runtime_error("operator() reached");
    edm4hep::ReconstructedParticleCollection output;
   
    std::cout << "==============Number of PFOs ============ "<< input.size()<< std::endl;

    for (const auto& pfo : input) {

        auto newPFO = output.create();

        //------------------------------------------------------
        // copy basic information
        //------------------------------------------------------

        //newPFO.setType(pfo.getType());

        newPFO.setEnergy(pfo.getEnergy());

        newPFO.setCharge(pfo.getCharge());

        newPFO.setMass(pfo.getMass());

        newPFO.setMomentum(pfo.getMomentum());

        newPFO.setReferencePoint(
            pfo.getReferencePoint());

        //------------------------------------------------------
        // copy relations
        //------------------------------------------------------

        for (auto track : pfo.getTracks())
           // newPFO.addToTracks(track);
           info() << "Track ID: collectionID = "
           << track.getObjectID().collectionID
           << ", index = "
           << track.getObjectID().index
           << endmsg;
        
        for (auto cluster : pfo.getClusters())
           // newPFO.addToClusters(cluster);
            info() << "Cluster ID: collectionID = "
           << cluster.getObjectID().collectionID
           << ", index = "
           << cluster.getObjectID().index
           << endmsg;


        for (auto particle : pfo.getParticles())
            newPFO.addToParticles(particle);

        info()
            << "Energy = "
            << pfo.getEnergy()
            << " Tracks = "
            << pfo.getTracks().size()
            << " Clusters = "
            << pfo.getClusters().size()
            << endmsg;
    }

    return output;
}
