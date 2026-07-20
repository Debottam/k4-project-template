#include "pandoraReading.h"

#include <GaudiKernel/MsgStream.h>

DECLARE_COMPONENT(pandoraReading)



double getFieldFromCompact() {
    dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
    const double position[3] = {0, 0, 0};      // position to calculate magnetic field at (the origin in this case)
    double magneticFieldVector[3] = {0, 0, 0}; // initialise object to hold magnetic field
    mainDetector.field().magneticField(position, magneticFieldVector); // get the magnetic field vector from DD4hep
    return magneticFieldVector[2] / dd4hep::tesla;                     // z component at (0,0,0)
  }



pandoraReading::pandoraReading( const std::string& name,ISvcLocator* svc)

    : Transformer(name,svc,{
              KeyValue{"InputPFOs", "PandoraPFOs"},
              KeyValue{"InputMCs", "MCParticles"},
              KeyValue{"RecoMCLinks", "MCTruthRecoLink"},
              KeyValue{"InputClusters", "PandoraClusters"},
              KeyValue{"InputTracks", "SiTracks_Refitted"}
          },

          {
              KeyValue{"OutputPFOs",
                       "RecoveredPFOs"}
          })
{
std::cout << "pandoraReading constructor called" << std::endl;
}

edm4hep::ReconstructedParticleCollection 
pandoraReading::operator()(const edm4hep::ReconstructedParticleCollection& input,
                           const edm4hep::MCParticleCollection& mcs,
                           const MCTruthRecoLinkCollection& links,
                           const edm4hep::ClusterCollection& clusters,
                           const edm4hep::TrackCollection& tracks) const {

    //==============DetectorRelatedOutput===============//

    auto detector = m_geoSvc->getDetector();
    std::cout << "Number of detectors = "
              << detector->detectors().size()
              << std::endl;

    for (const auto& det : detector->detectors()) {
        std::cout << "Detector: " << det.first << std::endl;
    }


    double eCalBarrelInnerR = detector->constantAsDouble("ECalBarrel_inner_radius")/dd4hep::mm;
     std::cout << "eCalBarrelInnerR: " << eCalBarrelInnerR << std::endl;

    double eCalEndCapInnerR = detector->constantAsDouble("ECalEndcap_inner_radius")/dd4hep::mm;
     std::cout << "eCalEndCapInnerR: " << eCalEndCapInnerR << std::endl;

    double eCalEndCapOuterR = detector->constantAsDouble("ECalEndcap_outer_radius")/dd4hep::mm;
     std::cout << "eCalEndCapOuterR: " << eCalEndCapOuterR << std::endl;

    double eCalEndCapInnerZ =detector->constantAsDouble("ECalEndcap_min_z")/dd4hep::mm;
     std::cout << "eCalEndCapInnerZ: " <<eCalEndCapInnerZ << std::endl; 

    double m_Bz = getFieldFromCompact();
    std::cout << "B field (T) is : " << m_Bz << std::endl;

    //double radius =
    //detector->constantAsDouble("ECalBarrel_inner_radius");

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
        std::cout <<"pfo.getPDG(): "<<pfo.getPDG() << '\n';
        std::cout <<"pfo.getObjectID(): "<<pfo.getObjectID().index << '\n';

        for (const auto& link : links) {
            auto reco = link.getFrom();
            auto mc   = link.getTo();
            if(reco.getObjectID().index==pfo.getObjectID().index){
              std::cout <<"mc.getPDG(): "<< mc.getPDG() << '\n';
              break;
            }
        }
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
