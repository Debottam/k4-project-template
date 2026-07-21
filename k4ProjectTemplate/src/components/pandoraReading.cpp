#include "pandoraReading.h"

#include <GaudiKernel/MsgStream.h>

DECLARE_COMPONENT(pandoraReading)



edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection,
                                                  const HelixClass_double& helixAtLastHit, double Bz) {

  edm4hep::TrackState trackState_AtCalorimeter = edm4hep::TrackState{};

  double posAtCalorimeter[] = {ecalProjection.GetX(), ecalProjection.GetY(), ecalProjection.GetZ()};
  int debug_lvl = 0;
  if (debug_lvl > 0) {
    std::cout << "Projection at calo: x, y, z, r = " << posAtCalorimeter[0] << " " << posAtCalorimeter[1] << " "
              << posAtCalorimeter[2] << " "
              << sqrt(posAtCalorimeter[0] * posAtCalorimeter[0] + posAtCalorimeter[1] * posAtCalorimeter[1])
              << std::endl;
  }

  // get extrapolated momentum from the helix with ref point at last hit
  double momAtCalorimeter[] = {0., 0., 0.};
  helixAtLastHit.getExtrapolatedMomentum(posAtCalorimeter, momAtCalorimeter);

  // produce new helix at calorimeter position
  auto helixAtCalorimeter = HelixClass_double();
  helixAtCalorimeter.Initialize_VP(posAtCalorimeter, momAtCalorimeter, helixAtLastHit.getCharge(), Bz);

  // fill the TrackState parameters
  trackState_AtCalorimeter.location = edm4hep::TrackState::AtCalorimeter;
  trackState_AtCalorimeter.D0 = helixAtCalorimeter.getD0();
  trackState_AtCalorimeter.phi = std::atan2(momAtCalorimeter[1], momAtCalorimeter[0]);
  trackState_AtCalorimeter.omega = helixAtCalorimeter.getOmega();
  trackState_AtCalorimeter.Z0 = helixAtCalorimeter.getZ0();
  trackState_AtCalorimeter.tanLambda = helixAtCalorimeter.getTanLambda();
  trackState_AtCalorimeter.referencePoint =
      edm4hep::Vector3f(posAtCalorimeter[0], posAtCalorimeter[1], posAtCalorimeter[2]);
  return trackState_AtCalorimeter;
}


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

    double eCalBarrelMaxZ = detector->constantAsDouble("ECalBarrel_half_length")/dd4hep::mm;
     std::cout << "eCalBarrelMaxZ: " << eCalBarrelMaxZ << std::endl;

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
        for (auto track : pfo.getTracks()){
           // newPFO.addToTracks(track);
           info() << "Track ID: collectionID = "
           << track.getObjectID().collectionID
           << ", index = "
           << track.getObjectID().index
           << endmsg;
            for (const auto& state : track.getTrackStates()) {

                 if (state.location == edm4hep::TrackState::AtLastHit) {
                          
                    double phi       = state.phi;
                    double omega     = state.omega;
                    double tanLambda = state.tanLambda;
                    double pt = 0.2998 * m_Bz / std::abs(omega);

                    double px = pt * std::cos(phi);
                    double py = pt * std::sin(phi);
                    double pz = pt * tanLambda;
                    auto ref = state.referencePoint;
                    double posAtLastHit[] = {ref.x, ref.y, ref.z};
                    double momAtLastHit[] = {px,py,pz};                    
                    int charge = (omega > 0.f) ? 1 : -1;
                    // produce new helix at last hit position
                    auto helixAtLastHit = HelixClass_double();
                    helixAtLastHit.Initialize_VP(posAtLastHit, momAtLastHit,charge, m_Bz);
                     // TrackState at Calorimeter
                    if (eCalBarrelInnerR > 0. || eCalEndCapInnerR > 0.) {
                        pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
                        pandora::CartesianVector secondBestECalProjection(0.f, 0.f, 0.f);
                        float minGenericTime(std::numeric_limits<float>::max());

                        // create helix to project
                        // rather than using parameters at production, better to use those from
                        // last hit
                        pandora::CartesianVector pos(ref.x, ref.y, ref.z);
                        pandora::CartesianVector mom(px,py,pz);
                        int charge = (omega > 0.f) ? 1 : -1;
                        const pandora::Helix helix(pos, mom, charge, m_Bz);
                        const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());
                        const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);
              
                        // First project to endcap
                        pandora::CartesianVector endCapProjection(0.f, 0.f, 0.f);
                        bool hasEndCapProjection(false);
                        if (eCalEndCapInnerR > 0) {
                            float genericTime(std::numeric_limits<float>::max());
                            const pandora::StatusCode statusCode(helix.GetPointInZ(static_cast<float>(signPz) * eCalEndCapInnerZ,
                                                                   referencePoint, endCapProjection, genericTime));
                            float x = endCapProjection.GetX();
                            float y = endCapProjection.GetY();
                            float r = std::sqrt(x * x + y * y);
                            if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime) &&
                                   (r >= eCalEndCapInnerR) && (r <= eCalEndCapOuterR)) {
                                minGenericTime = genericTime;
                                bestECalProjection = endCapProjection;
                                hasEndCapProjection = true;
                            }
                         }
                         // Then project to barrel surface(s), and keep projection
                         // if extrapolation is within the z acceptance of the detector
                         pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
                         bool hasBarrelProjection = false;
                         if (eCalBarrelInnerR > 0) {
                            float genericTime(std::numeric_limits<float>::max());
                            const pandora::StatusCode statusCode(
                            helix.GetPointOnCircle(eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));
                            if ((pandora::STATUS_CODE_SUCCESS == statusCode) &&
                                (std::fabs(barrelProjection.GetZ()) <= eCalBarrelMaxZ)) {
                                hasBarrelProjection = true;
                                if (genericTime < minGenericTime) {
                                   minGenericTime = genericTime;
                                   secondBestECalProjection = bestECalProjection;
                                   bestECalProjection = barrelProjection;
                                } 
                                else {
                                   secondBestECalProjection = barrelProjection;
                                }
                            }
                         }
                        edm4hep::TrackState trackState_AtCalorimeter = getExtrapolationAtCalorimeter(bestECalProjection, helixAtLastHit, m_Bz);
                    }

                 }
            }

        }
        
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
