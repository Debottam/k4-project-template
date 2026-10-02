#include "ElectronDresser.h"

#include <GaudiKernel/MsgStream.h>

#include <algorithm>

#include <tuple>

#include <cmath>

DECLARE_COMPONENT(ElectronDresser)

// struct for storing cluster and assocated track from PFO
struct ClusterTrackMap {
    edm4hep::Cluster cluster;
    std::vector<edm4hep::Track> tracks;
};

// struct for storing brem Recovered clustera and assocated track for target PFO
struct bremRecoveredClusterTrackMap {
    std::vector<edm4hep::Cluster> clusters;
    std::vector<edm4hep::Track> tracks;
};




struct ThreeVector {
    double x, y, z;
    
    ThreeVector() : x(0), y(0), z(0) {}
    ThreeVector(double x0, double y0, double z0) : x(x0), y(y0), z(z0) {}

    double mag() const { return std::sqrt(x*x + y*y + z*z); }
    double magXY() const { return std::sqrt(x*x + y*y); }
    double phi() const { return std::atan2(y,x); } // -pi to pi
    double theta() const { return std::acos(z / mag()); } // 0 to pi
    double eta() const { return -std::log(std::tan(theta()/2.)); }
  };

double deltaPhi(const double p1, const double p2) {
    double oo2pi = 1./(2.*M_PI);

    if (std::abs(p1-p2) <= M_PI)
      return p1 - p2;

    double n = std::round((p1-p2)*oo2pi);

    return p1 - p2 - n*2.*M_PI;
  }

/* decorating the track state at calorimeter level */
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

edm4hep::TrackState trackStateAtIP(const edm4hep::Track& track){
  edm4hep::TrackState trackState_AtIP;
  for (const auto& state : track.getTrackStates()) { 
     if (state.location == edm4hep::TrackState::AtIP) {
       trackState_AtIP = state;
       break;
     }
  }
  return trackState_AtIP;
}

/* track extrapolation to calorimeter cluster from the last track state */
edm4hep::TrackState helixAtLastTrack(const edm4hep::Track& track, double eCalBarrelInnerR, double eCalBarrelMaxZ,  double eCalEndCapInnerR, double eCalEndCapOuterR,
                                            double eCalEndCapInnerZ, double m_Bz){
                   
         std::cout<<"=======Inside helixAtLastTrack======== "<<std::endl;
         edm4hep::TrackState trackState_AtCalorimeter;
         for (const auto& state : track.getTrackStates()) {    //start of state loop
             
                 /* track needs to be extrapolated from its last hit */
                 if (state.location == edm4hep::TrackState::AtLastHit) {  // start of state.location if
                    
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
                    if (eCalBarrelInnerR > 0. || eCalEndCapInnerR > 0.) {   // start of eCalBarrel||eCalEndCap if
                       
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
                        if (eCalEndCapInnerR > 0) {                     // start of 1st eCalEndCapInnerR if
                            
                            float genericTime(std::numeric_limits<float>::max());
                            const pandora::StatusCode statusCode(helix.GetPointInZ(static_cast<float>(signPz) * eCalEndCapInnerZ,
                                                                   referencePoint, endCapProjection, genericTime));
                            float x = endCapProjection.GetX();
                            float y = endCapProjection.GetY();
                            float r = std::sqrt(x * x + y * y);
                            if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime) &&
                                   (r >= eCalEndCapInnerR) && (r <= eCalEndCapOuterR)) {  // start of eCalEndCapInnerR if
                                
                                minGenericTime = genericTime;
                                bestECalProjection = endCapProjection;
                                hasEndCapProjection = true;
                            }      // end of eCalEndCapInnerR if
                         }   // end of 1st eCalEndCapInnerR if
                         // Then project to barrel surface(s), and keep projection
                         // if extrapolation is within the z acceptance of the detector
                         pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
                         bool hasBarrelProjection = false;
                         if (eCalBarrelInnerR > 0) {              // start of eCalBarrelInnerR if
                            
                            float genericTime(std::numeric_limits<float>::max());
                            const pandora::StatusCode statusCode(
                            helix.GetPointOnCircle(eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));
                            if ((pandora::STATUS_CODE_SUCCESS == statusCode) &&
                                (std::fabs(barrelProjection.GetZ()) <= eCalBarrelMaxZ)) {  // start of barrel projection if
                                
                                hasBarrelProjection = true;
                                if (genericTime < minGenericTime) {   //start of generic time if
                                   
                                   minGenericTime = genericTime;
                                   secondBestECalProjection = bestECalProjection;
                                   bestECalProjection = barrelProjection;
                                }                                      //end of generic time if
                                else {   // start of else
                                   
                                   secondBestECalProjection = barrelProjection;
                                }        // end of else
                            }     // end of barrel projection if

                         }        // end of eCalBarrelInnerR if
                        
                        trackState_AtCalorimeter = getExtrapolationAtCalorimeter(bestECalProjection, helixAtLastHit, m_Bz);
                        auto ref_Calo = trackState_AtCalorimeter.referencePoint;
                         std::cout << "Track intersection1: ("
                                   << ref_Calo.x << ", "
                                   << ref_Calo.y << ", "
                                   << ref_Calo.z << ")\n";
                         
                    } // end of eCalBarrel||eCalEndCap if
                  break;
                 }  // end of state.location if
            }  //end of state loop

            return trackState_AtCalorimeter;
}


ThreeVector getMomentum(const edm4hep::TrackState& state, double Bz)
{
    double phi        = state.phi;
    double omega      = state.omega;
    double tanLambda  = state.tanLambda;

    double pt = (0.299792458 * Bz / std::abs(omega))*dd4hep::MeV/dd4hep::GeV;

    double px = pt * std::cos(phi);
    double py = pt * std::sin(phi);
    double pz = pt * tanLambda;

    return ThreeVector(px, py, pz);
}

double getFieldFromCompact() {
    dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
    const double position[3] = {0, 0, 0};      // position to calculate magnetic field at (the origin in this case)
    double magneticFieldVector[3] = {0, 0, 0}; // initialise object to hold magnetic field
    mainDetector.field().magneticField(position, magneticFieldVector); // get the magnetic field vector from DD4hep
    return magneticFieldVector[2] / dd4hep::tesla;                     // z component at (0,0,0)
  }


bool GetPointInZ(const double zPlane,        // z plane (endcap) coordinate
                   const ThreeVector& refPt,   // reference point
                   const ThreeVector& momVec,  // momentum
                   ThreeVector& intersectPt) { // output intersection point
    // protect from zero division
    if (std::abs(momVec.z) < std::numeric_limits<float>::epsilon()) {
      throw std::invalid_argument("BremRecoveryUtility: invalid momentum.z (zero division)");

      return false;
    }

    // check the assumption
    if (std::abs(refPt.z) > std::abs(zPlane)) {
      throw std::invalid_argument("BremRecoveryUtility: invalid refPoint.z (ref. point outside the endcap)");

      return false;
    }

    intersectPt.z = momVec.z > 0. ? std::abs(zPlane) : -std::abs(zPlane);
    //std::cout<<"intersectPt.z: "<<intersectPt.z<<std::endl;
    //std::cout<<"momVec.z: "<<momVec.z<<std::endl;
    const double genericTime = (intersectPt.z - refPt.z)/momVec.z;
    //std::cout<<"genericTime: "<<genericTime<<std::endl;
    //std::cout<<"momVec.x: "<<momVec.x<<std::endl;
    //std::cout<<"refPt.x: "<<refPt.x<<std::endl;
    intersectPt.x = refPt.x + momVec.x*genericTime;
    intersectPt.y = refPt.y + momVec.y*genericTime;

    return true;
  }


/* if track reherence point is not at origin */
bool GetPointOnCircle(const double radius,        // (barrel) radius
                        const ThreeVector& refPt,   // reference point
                        const ThreeVector& momVec,  // momentum
                        ThreeVector& intersectPt) { // output intersection point
    double refPtXY = refPt.magXY();
    double momVecXY = momVec.magXY();
    // protect from zero division
    if (std::abs(momVecXY) < std::numeric_limits<float>::epsilon()) {
      throw std::invalid_argument("BremRecoveryUtility: invalid pt (zero division)");

      return false;
    }

    // check the assumption
    if (std::abs(refPtXY) > radius) {
      throw std::invalid_argument("BremRecoveryUtility: invalid refPoint (ref. point outside the barrel)");

      return false;
    }
    // suppose a triangle OAP in the XY plane
    //      A        A is the reference point
    //      /\       O is the detector origin
    //     /  \      P is the intersection point
    //    ------     then the length OP = radius
    //   O      P

    // now do some middle school Euclidean geometry
    // cos(O+P) = refPt*momVec / |refPt|*|momVec|
    double normedInnerProd = (refPt.x*momVec.x + refPt.y*momVec.y) / (refPtXY*momVecXY);
    // sine rule : radius/sin(P+O) = |refPt|/sin(P)
    double sinP = refPtXY*std::sqrt(1. - normedInnerProd*normedInnerProd)/radius;
    double angOplusP = std::acos(normedInnerProd);
    double angP = std::asin(sinP); // always smaller than pi/2 by definition
    double angO = angOplusP - angP;
    // sine rule again : lenAP/sin(O) = |refPt|/sin(P)
    double lenAP = refPtXY*std::sin(angO)/sinP;
    double genericTime = lenAP/momVecXY; // this is what we wanted

    // now done with the geometry, come back to 3D
    intersectPt.x = refPt.x + momVec.x*genericTime;
    intersectPt.y = refPt.y + momVec.y*genericTime;
    intersectPt.z = refPt.z + momVec.z*genericTime;

    return true;
  }


/* if track reference at the origin */
bool GetPointOnCircleForRefAtOrigin(const double radius,        // (barrel) radius
                        const ThreeVector& refPt,   // reference point
                        const ThreeVector& momVec,  // momentum
                        ThreeVector& intersectPt) { // output intersection point
    double refPtXY = refPt.magXY();
    double momVecXY = momVec.magXY();
    // protect from zero division
    if (std::abs(momVecXY) < std::numeric_limits<float>::epsilon()) {
      throw std::invalid_argument("BremRecoveryUtility: invalid pt (zero division)");

      return false;
    }

    // check the assumption
    if (std::abs(refPtXY) > radius) {
      throw std::invalid_argument("BremRecoveryUtility: invalid refPoint (ref. point outside the barrel)");

      return false;
    }
    double genericTime = radius/momVecXY; // this is what we wanted

    // now done with the geometry, come back to 3D
    intersectPt.x = refPt.x + momVec.x*genericTime;
    intersectPt.y = refPt.y + momVec.y*genericTime;
    intersectPt.z = refPt.z + momVec.z*genericTime;

    return true;
  }


ElectronDresser::ElectronDresser( const std::string& name,ISvcLocator* svc)

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
ElectronDresser::operator()(const edm4hep::ReconstructedParticleCollection& input,
                           const edm4hep::MCParticleCollection& mcs,
                           const MCTruthRecoLinkCollection& links,
                           const edm4hep::ClusterCollection& clusters,
                           const edm4hep::TrackCollection& tracks) const {

    //==============DetectorRelatedOutput===============//

    auto detector = m_geoSvc->getDetector();
    /*std::cout << "Number of detectors = "
              << detector->detectors().size()
              << std::endl;*/

   // for (const auto& det : detector->detectors()) {
   //     std::cout << "Detector: " << det.first << std::endl;
   // }


    double eCalBarrelInnerR = detector->constantAsDouble("ECalBarrel_inner_radius")/dd4hep::mm;
     //std::cout << "eCalBarrelInnerR: " << eCalBarrelInnerR << std::endl;

    double eCalBarrelMaxZ = detector->constantAsDouble("ECalBarrel_half_length")/dd4hep::mm;
     //std::cout << "eCalBarrelMaxZ: " << eCalBarrelMaxZ << std::endl;

    double eCalEndCapInnerR = detector->constantAsDouble("ECalEndcap_inner_radius")/dd4hep::mm;
     //std::cout << "eCalEndCapInnerR: " << eCalEndCapInnerR << std::endl;

    double eCalEndCapOuterR = detector->constantAsDouble("ECalEndcap_outer_radius")/dd4hep::mm;
     //std::cout << "eCalEndCapOuterR: " << eCalEndCapOuterR << std::endl;

    double eCalEndCapInnerZ =detector->constantAsDouble("ECalEndcap_min_z")/dd4hep::mm;
     //std::cout << "eCalEndCapInnerZ: " <<eCalEndCapInnerZ << std::endl; 

    double m_Bz = getFieldFromCompact();
    //std::cout << "B field (T) is : " << m_Bz << std::endl;

    //double radius =
    //detector->constantAsDouble("ECalBarrel_inner_radius");

    //throw std::runtime_error("operator() reached");
    edm4hep::ReconstructedParticleCollection output;
   
    // Forming the cluster->track map
    std::vector<ClusterTrackMap> clusterTrackMaps;

    std::vector<edm4hep::Cluster> tempClusters;

    std::vector<std::vector<edm4hep::Cluster>> bremRecoverdClusters;

    std::cout << "==============start of new event ============ \n"<< std::endl;

    std::cout << "==============Number of PFOs ============ "<< input.size()<< std::endl;


    float Total_MC_Energy = 0.0;
    float Total_PFO_Energy = 0.0;
    float Total_MC_Truth_Energy = 0.0;
    float Total_Input_PFO_Cluster_Energy = 0.0;
    float Total_Input_PFO_Track_Energy = 0.0;

    for (const auto& mc : mcs) {

        //if(!mc.isCreatedInSimulation()){
           Total_MC_Truth_Energy= Total_MC_Truth_Energy+mc.getEnergy();
        //}
        std::cout << "=======MC particle Truth energy ==== "
                  << mc.getEnergy()
                  << " GeV" << std::endl;
        std::cout <<"mc.getPDG(): "<< mc.getPDG() <<std::endl;
        std::cout <<"mc.isCreatedInSimulation(): "<< mc.isCreatedInSimulation() <<std::endl;
    }
     std::cout << "=======Total MC particle Truth energy ==== "
                  << Total_MC_Truth_Energy
                  << " GeV" << std::endl;

    for (const auto& pfo : input) {  //PFO loop starts

        info()
            << "PFO Energy = "
            << pfo.getEnergy()
            << " Tracks = "
            << pfo.getTracks().size()
            << " Clusters = "
            << pfo.getClusters().size()
            << endmsg;

           Total_PFO_Energy = Total_PFO_Energy+ pfo.getEnergy();
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
              std::cout <<"mc.isCreatedInSimulation(): "<< mc.isCreatedInSimulation() << '\n';
              std::cout << "MC energy = "<< mc.getEnergy() << '\n';
              Total_MC_Energy = Total_MC_Energy+mc.getEnergy();
              break;
            }
          
        }
        
    
        // Forming the cluster->track map
        //std::vector<ClusterTrackMap> clusterTrackMaps;

        for (const auto& cluster : pfo.getClusters()) {  // PFO cluster loop starts

           if (!cluster.isAvailable())
              continue;

           Total_Input_PFO_Cluster_Energy = Total_Input_PFO_Cluster_Energy+cluster.getEnergy();
           tempClusters.push_back(cluster);           

           ClusterTrackMap entry;
           entry.cluster = cluster;

           for (const auto& track : pfo.getTracks()) { // PFO track loop starts
              if (track.isAvailable()) {
                  auto state = track.getTrackStates(edm4hep::TrackState::AtIP);
                  const auto momVec = getMomentum(state, m_Bz);

                  float p = std::sqrt(momVec.x*momVec.x + momVec.y*momVec.y + momVec.z*momVec.z);

                  // massless approximation

                  Total_Input_PFO_Track_Energy = Total_Input_PFO_Track_Energy+p;
                  entry.tracks.push_back(track);
              }
           }  //  PFO track loop ends

           clusterTrackMaps.push_back(std::move(entry));
        }  //  PFO cluster loop starts ends

     } // PFO  loop ends
      std::cout <<"Total_MC_Energy: "<<Total_MC_Energy <<std::endl;
      std::cout <<"Total_PFO_Energy: "<<Total_PFO_Energy <<std::endl;
      std::cout <<"diff_PFO_mc_Energy: "<<std::abs(Total_PFO_Energy-Total_MC_Energy) <<std::endl;

     // sorting the cluster track map
       std::sort(clusterTrackMaps.begin(), clusterTrackMaps.end(),
       [](const ClusterTrackMap& a, const ClusterTrackMap& b) {
              return a.cluster.getEnergy() > b.cluster.getEnergy(); 
       });


       std::sort(tempClusters.begin(),tempClusters.end(),
       [](const edm4hep::Cluster& a, const edm4hep::Cluster& b){
            return a.getEnergy() > b.getEnergy();
       });
      //sorting the tracks by distance from the cluster

       //std::cout << "==============Size of clusterTrackMaps ============ "<< clusterTrackMaps.size()<< std::endl;

       // std::cout << "==============Size of tempClusters ============ "<< tempClusters.size()<< std::endl;

       for (auto& entry : clusterTrackMaps) {  //  Looping over clusterTrackMaps starts

          if (entry.tracks.size() <= 1)
            continue;

          std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!there are more than one track==============="<<std::endl;

          const auto clusterPos = entry.cluster.getPosition();

          std::sort(entry.tracks.begin(), entry.tracks.end(),
            [&](const edm4hep::Track& a, const edm4hep::Track& b) {

               auto getDistance = [&](const edm4hep::Track& track) {

                 edm4hep::TrackState state=helixAtLastTrack(track, eCalBarrelInnerR, eCalBarrelMaxZ,  eCalEndCapInnerR, eCalEndCapOuterR,
                                             eCalEndCapInnerZ, m_Bz);

                 auto ref = state.referencePoint;
                 float dx = clusterPos.x - ref.x;
                 float dy = clusterPos.y - ref.y;
                 float dz = clusterPos.z - ref.z;

                 return std::sqrt(dx*dx + dy*dy + dz*dz);
               };

             return getDistance(a) < getDistance(b); // closest first
           });
        }  //Looping over clusterTrackMaps ends

        //std::vector<edm4hep::Cluster> bremCands;
        std::vector<int> clusterIndex;
        std::vector<bremRecoveredClusterTrackMap> bremRecoveredClusterTrackMaps;

        for (const auto& entry : clusterTrackMaps) {  //  Looping over clusterTrackMaps starts for brem recovery

            const auto& cluster = entry.cluster;

            //std::cout << "Cluster energy = " << cluster.getEnergy() << std::endl;

            if (entry.tracks.empty()) {
             //std::cout << "  No associated tracks." << std::endl;
             continue;
            }

            if (!entry.tracks[0].isAvailable())
                continue;

               // Use the track
               /*std::cout << "  Track has "
                  << entry.tracks[0].getTrackStates().size()
                  << " track states."
                  << std::endl;*/

            bool aStatus = false;
            const auto& track1st = entry.tracks[0];
            const auto& trackState1st = trackStateAtIP(track1st);
            const auto refPt = ThreeVector(trackState1st.referencePoint.x, trackState1st.referencePoint.y, trackState1st.referencePoint.z);
            const auto momVec = getMomentum( trackState1st, m_Bz);
            double zPlane = eCalEndCapInnerZ;
            double radius = eCalBarrelInnerR;
            ThreeVector intersectZ;
            aStatus = GetPointInZ(zPlane,refPt,momVec,intersectZ);
            //std::cout<<"refPt.x: "<<refPt.x<<", refPt.y: "<<refPt.y<<", refPt.z: "<<refPt.z<<std::endl;
            //std::cout<<"momVec.x: "<<momVec.x<<", momVec.y: "<<momVec.y<<", momVec.z: "<<momVec.z<<std::endl;

            if (!aStatus) {
              throw std::runtime_error("BremRecoveryAlgorithm: failed to run BremRecoveryUtility::GetPointInZ!");
            }
            //std::cout<<"intersectZ.x: "<<intersectZ.x<<", intersectZ.y: "<<intersectZ.y<<", intersectZ.z: "<<intersectZ.z<<std::endl;
            
            ThreeVector finalIntersectPt = intersectZ;

           if (radius <= 0.) { // this should never happen
              throw std::invalid_argument("BremRecoveryAlgorithm: ECAL barrel inner radius < 0!");
           } 

           //std::cout<<"intersectZ.magXY(): "<<intersectZ.magXY()<<" radius: "<<radius<<std::endl;

           if (intersectZ.magXY() > radius){ // the tangent intersects the barrel before reaching the endcap - recalculate it
              if(std::abs(refPt.x) < std::numeric_limits<float>::epsilon() &&
                 std::abs(refPt.y) < std::numeric_limits<float>::epsilon() &&
                 std::abs(refPt.z) < std::numeric_limits<float>::epsilon()){

                 aStatus = GetPointOnCircleForRefAtOrigin(radius,refPt,momVec,finalIntersectPt);

              }
              else{
                 aStatus = GetPointOnCircle(radius,refPt,momVec,finalIntersectPt);
              }
           }

          if (!aStatus) {
            throw std::runtime_error("BremRecoveryAlgorithm: failed to run BremRecoveryUtility::GetPointOnCircle!");
          }
          //std::cout<<"====aStatus======="<<aStatus<<std::endl;        
          //std::cout<<"finalIntersectPt.x: "<<finalIntersectPt.x<<", finalIntersectPt.y: "<<finalIntersectPt.y<<", finalIntersectPt.z: "<<finalIntersectPt.z<<std::endl;

          double dEtaMargin = m_dEtaMargin;
          double dPhiMargin = m_dPhiMargin;

          const auto& clusterPos1st = cluster.getPosition();
          const auto vecClus1 = ThreeVector(clusterPos1st.x,
                                            clusterPos1st.y,
                                            clusterPos1st.z);
          const double eta1st = vecClus1.eta();
          const double phi1st = vecClus1.phi();
          const double maxEta = eta1st + dEtaMargin;
          const double minEta = eta1st - dEtaMargin;
          const double bremBandPhi = deltaPhi(finalIntersectPt.phi(),phi1st);
          const double maxDPhi = bremBandPhi > 0. ? bremBandPhi + dPhiMargin : dPhiMargin;
          const double minDPhi = bremBandPhi > 0. ? -dPhiMargin : bremBandPhi - dPhiMargin;

           std::vector<edm4hep::Cluster> bremCands;
           //bremCands.push_back(cluster);
          // loop over clusters again
          for (const auto& cluster2nd:tempClusters) { //  Looping over temporary clusters
              
             if (cluster2nd.getObjectID().index == cluster.getObjectID().index){  // start if 
              // std::cout<<"cluster2nd.getObjectID().index: "<<cluster2nd.getObjectID().index<<std::endl;
              // std::cout<<"cluster.getObjectID().index: "<<cluster.getObjectID().index<<std::endl;
               continue;
             }      // end if 
             if (!cluster2nd.isAvailable())
                continue;

      // TODO apply a templated e/gamma ID
      // if (not pass ID...)
      //   continue;

      // check the strip
            //std::cout<<"cluster2nd.getObjectID().index: "<<cluster2nd.getObjectID().index<<std::endl;
            //std::cout<<"cluster.getObjectID().index: "<<cluster.getObjectID().index<<std::endl;
            const auto& clusterPos2nd = cluster2nd.getPosition();
            const auto vecClus2 = ThreeVector(clusterPos2nd.x,
                                              clusterPos2nd.y,
                                              clusterPos2nd.z);
            const double eta2nd = vecClus2.eta();
            const double phi2nd = vecClus2.phi();
            const double finalDeltaPhi = deltaPhi(phi2nd,phi1st);

      // 3. Incorporate the selected clusters and recluster the electron cluster
            if (minEta < eta2nd && eta2nd < maxEta && minDPhi < finalDeltaPhi && finalDeltaPhi < maxDPhi) {  //start if
                std::cout<<"===========found a brem candidate============"<<std::endl;
      
                if (std::find(clusterIndex.begin(), clusterIndex.end(), cluster2nd.getObjectID().index) == clusterIndex.end()){ //start if
                    bremCands.push_back(cluster2nd);
                    clusterIndex.push_back(cluster2nd.getObjectID().index);
                }//end if
                if (std::find(clusterIndex.begin(), clusterIndex.end(), cluster.getObjectID().index) == clusterIndex.end()){ //start if
                    bremCands.push_back(cluster);
                    clusterIndex.push_back(cluster.getObjectID().index);
                }//end if
             }//end if
          } //  end Looping over temporary clusters
          std::cout<<"===================size of bremCands==================:::  "<<bremCands.size()<<std::endl;
          if(bremCands.size()!=0)
             bremRecoverdClusters.push_back(bremCands);
          else{
             std::vector<edm4hep::Cluster> placeHolder;
             placeHolder.push_back(cluster);
             bremRecoverdClusters.push_back(placeHolder);
          }  
         //std::cout<<"===================size of bremCandsVector==================:::  "<<bremRecoverdClusters.size()<<std::endl;
        }  //  Looping over clusterTrackMaps ends for brem recovery
         //std::cout<<"===================size of clusterIndex==================:::  "<<clusterIndex.size()<<std::endl;
          for (int x : clusterIndex)
          {
            std::cout << x << " ";
          }
         std::cout << std::endl;
         //std::cout<<"===================size of bremCandsVector==================:::  "<<bremRecoverdClusters.size()<<std::endl;

         // for all selected and recovered clusters

         bremRecoveredClusterTrackMap recoveredMap;

         // Keep track IDs here so that tracks are not duplicated
         std::unordered_set<int> recoveredTrackIndices;

         for (const auto& ctm : clusterTrackMaps) {

           int clusterID = ctm.cluster.getObjectID().index;

           // Is this cluster one of the selected clusters?
           if (std::find(clusterIndex.begin(), clusterIndex.end(), clusterID)!= clusterIndex.end()){
           //if (selectedClusterIndices.count(clusterIndex)) {

             // Store the cluster
              recoveredMap.clusters.push_back(ctm.cluster);

             // Store associated tracks, but only once
              for (const auto& track : ctm.tracks) {
                //std::cout << "=========== Inside CTM tracks brem ============ " << std::endl;
                int trackIndex = track.getObjectID().index;
                //std::cout << "ctm.tracks = " << trackIndex << std::endl;

                if (recoveredTrackIndices.insert(trackIndex).second) {
                   recoveredMap.tracks.push_back(track);
                }
              }
           }
         }
         
       // --------------------------------------------------
       //  Add the recovered/selected group
       // --------------------------------------------------

       if (!recoveredMap.clusters.empty()) {
        bremRecoveredClusterTrackMaps.push_back(
        std::move(recoveredMap)
        );
      }

      // --------------------------------------------------
      //       Store every remaining cluster separately
      // --------------------------------------------------

      for (const auto& ctm : clusterTrackMaps) {

       int clusterID = ctm.cluster.getObjectID().index;

        // Skip clusters already included in recoveredMap
        if (std::find(clusterIndex.begin(), clusterIndex.end(), clusterID)!= clusterIndex.end()) {
          continue;
        }

        bremRecoveredClusterTrackMap remainingMap;

        // One cluster
        remainingMap.clusters.push_back(ctm.cluster);

        // Its associated tracks
        for (const auto& track : ctm.tracks) {
          //std::cout << "=========== Inside CTM tracks============ " << std::endl;
          //std::cout <<"track.getObjectID().index:::"<<track.getObjectID().index<<std::endl;
          remainingMap.tracks.push_back(track);
        }

        bremRecoveredClusterTrackMaps.push_back(
          std::move(remainingMap)
        );
      }

       for (const auto& entry : bremRecoveredClusterTrackMaps) {  //  Looping over clusterTrackMaps starts for brem recovery

            for(const auto& track: entry.tracks){

            //std::cout << "Retrieved Clusters ID = " << track.getObjectID().index << std::endl;
            }
       }
      
    //edm4hep::ReconstructedParticleCollection output;

    
    for (const auto& bremMap : bremRecoveredClusterTrackMaps) { // for bremMap starts

      // Create one ReconstructedParticle
      auto particle = output.create();
      
      float clusterEnergy = 0.0;
      float pfoEnergy = 0.0;

      // Add all associated tracks
      for (const auto& track : bremMap.tracks) { // for track starts 
        particle.addToTracks(track);
       }                                         // for track ends

      // Add all associated clusters
      for (const auto& cluster : bremMap.clusters) { // for cluster starts 
        particle.addToClusters(cluster);
        clusterEnergy += cluster.getEnergy();
        for (const auto& pfo : input) {              // for pfo starts

            for (const auto& pfoCluster : pfo.getClusters()) {  // for pfocluster starts

                if (pfoCluster.getObjectID().index ==
                    cluster.getObjectID().index) { // if starts 

                    pfoEnergy += pfo.getEnergy();

                    //std::cout << "Found PFO energy = "
                    //          << pfo.getEnergy() << std::endl;
                    
                    break;
                }                                  // if ends
            }                                       // for pfocluster ends
      }                                              // for pfo ends
    
    }                                                           // for cluster ends
     std::cout << "clusterEnergy = " << clusterEnergy << std::endl;                                     
     std::cout << "pfoEnergy = " << pfoEnergy << std::endl;
     particle.setEnergy(clusterEnergy);
    }// for bremMap ends

     

    std::cout << "====Number of PFO particles created:====="<< output.size() << std::endl;
    std::cout <<"\n\n"<< std::endl;
 
    float Total_PFO_Output_energy = 0.0;
    int outputSize = 0;
    for (const auto& pfoOut : output) {  //PFO loop starts
      if (outputSize==0){
       Total_PFO_Output_energy=Total_PFO_Output_energy+pfoOut.getEnergy();
      }
      outputSize=outputSize+1;
      std::cout << "PFO Output Energy = "<< pfoOut.getEnergy()<<" GeV"<<std::endl;
    }

    std::cout << "Total PFO Output Energy = "<< Total_PFO_Output_energy<<" GeV"<<std::endl;

    std::cout <<"Total_Input_PFO_Cluster_Energy = "<<Total_Input_PFO_Cluster_Energy<<" GeV"<<std::endl;

     std::cout <<"Total_Input_PFO_Track_Energy = "<<Total_Input_PFO_Track_Energy<<" GeV"<<std::endl;

    std::cout << "PFO resolution = "<< (Total_Input_PFO_Track_Energy-Total_MC_Truth_Energy)/Total_MC_Truth_Energy<<std::endl;

    std::cout << "::::::::::::::::::::::::end of an event ::::::::::::::::::::::::::::\n"<< std::endl;
    return output;
}
  
