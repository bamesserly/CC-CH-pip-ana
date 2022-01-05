#ifndef MICHEL_H
#define MICHEL_H

#include "CVUniverse.h"
#include "Cluster.h"

class Michel;

class Michel {
 public:
  // constructors
  Michel(const CVUniverse& univ, int ci);
  // Functions
  Michel(){};  // fill in most basic stuff in "generic info" . this is default
               // constructor you need it if you have public data members
  void DoMoreInitializationAndProcessing();  // fill in more complicated stuff
                                             // in "generic info"
  void DoMatching();                         // fill in "best matching stuff"
  void DoesMichelMatchVtx(const CVUniverse& univ);  // GEts info for Vtx Match
  void DoesMichelMatchClus(
      const CVUniverse& univ);  // Gets info for ClusterMatch
  void GetBestMatch();  // get type for best match out of all four saved matches
  void GetPionAngle(
      const CVUniverse&
          univ);  // Function to calculate pion angle with respect to beam.
                  // Simply gets angle between best michel endpoint and vertex.

  // std::vector<Michel*> CreateMichels(CVUniverse& univ);
  // generic info
  std::vector<double> up_location;    // upstream location 0 X 1 U 2 V 3 Z
  std::vector<double> down_location;  // downstream location
  double m_x1 = 9999.;                // Michel Endpoint 1 x
  double m_x2 = 9999.;                // Michel Endpoint 2 x
  double m_y1 = 9999.;                // Michel Endpoint 1 y
  double m_y2 = 9999.;                // Michel Endpoint 2 y
  double m_u1 = 9999.;                // Michel Endpoint 1 u
  double m_u2 = 9999.;                // Michel Endpoint 2 u
  double m_v1 = 9999.;                // Michel Endpoint 1 v
  double m_v2 = 9999.;                // Michel Endpoint 2 v
  double m_z1 = 9999.;                // Mihel Endpoint z 1
  double m_z2 = 9999.;                // Michel Endpoiint z2
  double energy = -999.;              // Michel energy
  double time = -999.;                // Michel Time
  int is_fitted = -1;                 // Is the Michel fitted? 0 no. 1 yes.
  // Following are 2D distances (were stored in vectors but now as explicit data
  // members)
  double up_to_vertex_XZ = 9999.;
  double up_to_vertex_UZ = 9999.;
  double up_to_vertex_VZ = 9999.;
  double down_to_vertex_XZ = 9999.;
  double down_to_vertex_UZ = 9999.;
  double down_to_vertex_VZ = 9999.;
  double down_to_clus_XZ = 9999.;
  double down_to_clus_UZ = 9999.;
  double down_to_clus_VZ = 9999.;
  double up_to_clus_XZ = 9999.;
  double up_to_clus_UZ = 9999.;
  double up_to_clus_VZ = 9999.;
  // Michel End point to Vertex distance  (up = endpoint 1 and down = endpoint
  // 2) TODO: actually find out which end point is upstream or downstream
  double up_to_vertex_dist3D = 9999.;
  double down_to_vertex_dist3D = 9999.;
  // Maybe keep a vector of clusters that matched to each endpoint?
  std::vector<Cluster> cluster_to_up_match;
  std::vector<Cluster> cluster_to_down_match;
  // 3D distances between cluster and Michel
  double up_to_cluster_dist3D = 9999.;  // Distance between vertex and cluster
                                        // that was matched to endpoint 1
  double down_to_cluster_dist3D =
      9999.;  // Distance between vertex and cluster matched to endpoint 2
  double up_clus_michel_dist3D =
      9999.;  // Distance between Michel endpoint 1 and clusters
  double down_clus_michel_dist3D =
      9999.;  // Distance between Michel endpoint 2 and clusters
  double up_clus_michvtx_dist3D =
      9999.;  // Distance between the Michel end point 1 that matched to
              // clusters and the vertex - this will be used as pion range
  double down_clus_michvtx_dist3D =
      9999.;  // Distance between the Michel endpoint 2 that matched to clusters
              // and the vertex - this will be used as pion range
  double vtx_michel_timediff = 9999.;

  double overlay_fraction =
      -1.0;  // Overlay fraction of the Michel, Default if Data. 0 if MC 1 if
             // Data... (maybe some events in between?)
  int nclusters = 0;      // number of (non-muon) clusters in the primary event
  int vtx_endpoint = 0;   // 1 or 2 for which Michel end point is closest
  int clus_endpoint = 0;  // 1 or 2 for which Michel endpoint is closest
  // best matching stuff
  // enum *best_cluster_match; // just tells you which of the four matches are
  // the best match enum BestClusterMatch {kUpVtxMatch, kDownVtxMatch,
  // kUpClusMatch, kDownClusMatch, kNClusterMatches}; the following is in place
  // until i figure out how to use the enum.
  int BestMatch = 0;  // 0 = null match, 1= kUpVtxMatch, 2 = kDownVtxMatch, 3 =
                      // kUpClusMatch, 4=  kDownClusMatch,
  int SecondBestMatch = 0;
  int tuple_idx;  // index of the Michel out of all the michels saved in the
                  // tuple
  double Best3Ddist =
      9999.;  // Best 3D distance out of eitehr a vertex or a cluster match
  // best 2D distance for the best type of match
  double best_XZ = 9999.;
  double best_UZ = 9999.;
  double best_VZ = 9999.;
  // Want to save the index of the clusters that the Michel best matched to.
  int xclus_idx;
  int uclus_idx;
  int vclus_idx;

  // True initial position of the michel  TODO: initial the other Michel truth
  // member data here  (energy, time, momentum etc)
  double true_angle = 9999.;
  double true_initialx = 9999.;
  double true_initialy = 9999.;
  double true_initialz = 9999.;
  double true_e = 9999.;
  double true_p = 9999.;
  double true_pdg = -1.0;
  int true_parentid = -1;
  int true_parentpdg = -1;
  double true_parent_energy = -9999.;
  double true_parent_p = -9999.;
  double true_parent_px = -9999.;
  double true_parent_py = -9999.;
  double true_parent_pz = -9999.;
  double true_parent_xi = -9999.;
  double true_parebt_yi = -9999.;
  double true_parent_zi = -9999.;
  double true_parent_xf = -9999.;
  double true_parent_yf = -9999.;
  double true_parent_zf = -9999.;
  // the following member data were created to investigate my weird convoluted
  // way of geting x and y values. Probably dont need them now. // TODO: check
  // and remove the following member data
  double best_angle = -9999.;
  double up_clus_x = 9999.;
  double up_clus_y = 9999.;
  double up_clus_z = 9999.;
  double down_clus_x = 9999.;
  double down_clus_y = 9999.;
  double down_clus_z = 9999.;
  double up_vtx_x = 9999.;
  double up_vtx_y = 9999.;
  double up_vtx_z = 9999.;
  double down_vtx_x = 9999.;
  double down_vtx_y = 9999.;
  double down_vtx_z = 9999.;
  double is_overlay = -1;
  double pionKE = -9999.;
  // Adding the following member data to determine the true endpoint position of
  // the Michel
  int trueEndpoint =
      -1;  // 0 = Overlay Michel, 1 = Endpoint 1 is correct intial position of
           // Michel, 2 = Endpoint 2 is correct Initial Position of Michel
  // ~DeleteMichel(){delete this;};  // This is going to be the main destructor.
  int recoEndpoint = -1;  // -1 is default/NULL like above. 1 = Endpoint 1 is
                          // better match, 2 = Endpoint 2 is better match

  std::vector<int> passable_matchtype{
      -1, -1, -1,
      -1};  // This vector will contain a value for each match type -1 or the
            // matchtype 1, 2, 3, 4 for UpVtx, DownVTx, Upclus, DownClus
            // depending of that match type passes our 2D distance cut. (if
            // distance is large, then it'll pass all of them).
};

// Setting the Michel data members directly from CV universe. No calculations
// made at this stage
Michel::Michel(const CVUniverse& univ, int ci) {
  energy = univ.GetVecElem("FittedMichel_michel_energy", ci);
  time = univ.GetVecElem("FittedMichel_michel_time", ci) / pow(10, 3);
  is_fitted = univ.GetVecElem("FittedMichel_michel_fitPass", ci);
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_x1", ci));
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_u1", ci));
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_v1", ci));
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_z1", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_x2", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_u2", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_v2", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_z2", ci));
  m_x1 = univ.GetVecElem("FittedMichel_michel_x1", ci);
  m_y1 = univ.GetVecElem("FittedMichel_michel_y1", ci);
  m_u1 = univ.GetVecElem("FittedMichel_michel_u1", ci);
  m_v1 = univ.GetVecElem("FittedMichel_michel_v1", ci);
  m_z1 = univ.GetVecElem("FittedMichel_michel_z1", ci);
  m_x2 = univ.GetVecElem("FittedMichel_michel_x2", ci);
  m_y2 = univ.GetVecElem("FittedMichel_michel_y2", ci);
  m_u2 = univ.GetVecElem("FittedMichel_michel_u2", ci);
  m_z2 = univ.GetVecElem("FittedMichel_michel_z2", ci);
  m_v2 = univ.GetVecElem("FittedMichel_michel_v2", ci);
  nclusters = univ.GetInt("FittedMichel_cluster_view_sz");
  overlay_fraction = univ.GetVecElem("FittedMichel_michel_datafraction", ci);

  true_initialx =
      univ.GetVecElem("FittedMichel_reco_micheltrajectory_initialx", ci);
  true_initialy =
      univ.GetVecElem("FittedMichel_reco_micheltrajectory_initialy", ci);
  true_initialz =
      univ.GetVecElem("FittedMichel_reco_micheltrajectory_initialz", ci);
  is_overlay = univ.GetVecElem("FittedMichel_michel_isoverlay", ci);
  true_e = univ.GetVecElem("FittedMichel_reco_micheltrajectory_energy", ci);
  true_pdg = univ.GetVecElem("FittedMichel_reco_micheltrajectory_pdg", ci);
  true_parentpdg = univ.GetVecElem("FittedMichel_true_primaryparent_pdg", ci);
  true_parentid =
      univ.GetVecElem("FittedMichel_true_primaryparent_trackID", ci);
  true_p = univ.GetVecElem("FittedMichel_reco_micheltrajectory_momentum", ci);

  double true_parentp =
      univ.GetVecElem("FittedMichel_true_primaryparent_momentum", ci);
  double true_parente =
      univ.GetVecElem("FittedMichel_true_primaryparent_energy", ci);
  double mass = mass = sqrt(pow(true_parente, 2) - pow(true_parentp, 2));
  pionKE = true_parente - mass;

  double true_parentpx =
      univ.GetVecElem("FittedMichel_true_primaryparent_momentumx", ci);
  double true_parentpy =
      univ.GetVecElem("FittedMichel_true_primaryparent_momentumy", ci);
  double true_parentpz =
      univ.GetVecElem("FittedMichel_true_primaryparent_momentumz", ci);

  TVector3 truep(true_parentpx, true_parentpy, true_parentpz);
  double true_theta = truep.Theta();
  true_angle = true_theta;  //*TMath::RadToDeg();

  // if (overlay_fraction < 0.5) std::cout << "True Parent of Michel is PDG: "
  // <<  true_parentpdg <<  " And Parent trackID: "  << true_parentid <<
  // std::endl;
  double end1diff =
      abs(true_initialz -
          m_z1);  // This gives a value for determining how close the
                  // reconstructed endpoint of the michel is to the true intial
                  // endpoint (the start point of where the michel decayed from)
  double end2diff =
      abs(true_initialz -
          m_z2);  // this is for endpoint 2. If you compare this to the endpoint
                  // that gets matched to a verted or cluster, you can determine
                  // which type of match ends up getting correcct matches or
                  // wrong matches.

  if (overlay_fraction > 0.5)
    trueEndpoint = 0;
  else if (true_parentpdg == 211 && end1diff < end2diff)
    trueEndpoint = 1;
  else if (true_parentpdg == 211 && end2diff < end1diff)
    trueEndpoint = 2;
  if (is_fitted == 1) {         // Do theMatching for Fitted Michels
    DoesMichelMatchVtx(univ);   // GEts info for Vtx Match
    DoesMichelMatchClus(univ);  // Gets info for ClusterMatch
    GetBestMatch();
    GetPionAngle(univ);
  }
}

// This function will get the angle between Michel endpoint that was matched and
// the vertex

void Michel::GetPionAngle(const CVUniverse& univ) {
  double vtx_x = univ.GetVertex().X();  // mm
  double vtx_y = univ.GetVertex().Y();  // mm
  double vtx_z = univ.GetVertex().Z();  // mm

  TVector3 vtx(vtx_x, vtx_y, vtx_z);
  TVector3 endpoint;

  if (this->BestMatch == 1 || this->BestMatch == 3)
    endpoint.SetXYZ(this->m_x1, this->m_y1, this->m_z1);
  else if (this->BestMatch == 2 || this->BestMatch == 4)
    endpoint.SetXYZ(this->m_x2, this->m_y2, this->m_z2);
  else
    endpoint.SetXYZ(9999., 9999., 9999.);

  TVector3 range = endpoint - vtx;
  double angle = univ.thetaWRTBeam(range.x(), range.y(), range.z());
  this->best_angle = angle;  //*TMath::RadToDeg(); // in Degrees
}

// This function will get an integer for the best match type of the Michel.
// It compares distance between MIchel and whatever it's best match is to find
// the Best type of Michel for a single Michel.

void Michel::GetBestMatch() {
  int upvtxmatch = 0;
  int downvtxmatch = 0;
  int upclusmatch = 0;
  int downclusmatch = 0;

  // if (this->up_to_vertex_dist3D != NULL && )

  // std::cout << "GET BEST MATCH FOR THIS MICHEL " << std::endl;
  // This is setting the values for which endpoint is a better match for each
  // type. TODO: Revise this function. There has to be a better way to compare
  // distances than what I wrote.

  std::vector<double> distances{
      this->up_to_vertex_dist3D, this->down_to_vertex_dist3D,
      this->up_clus_michel_dist3D, this->down_clus_michel_dist3D};
  // std::cout << "The distances for this michel are upvtx, downvtx, upclus,
  // downclus " << distances[0] << " - " << distances[1] << " - " << distances[2]
  // << " - " << distances[3] << " - "  << std::endl;
  std::sort(distances.begin(), distances.end());

  // std::cout << "The distances after sorting " << distances[0] << " - " <<
  // distances[1] << " - " << distances[2] << " - " << distances[3] << " - " <<
  // std::endl;
  // This bit of code will try to force the true Endpoint of hte Michel
  /*
   if (trueEndpoint == 1)
   {
    if (this->up_clus_michel_dist3D > this->up_to_vertex_dist3D){
      upvtxmatch = 1;
      this->best_XZ = this->up_to_vertex_XZ;
      this->best_UZ = this->up_to_vertex_UZ;
      this->best_VZ = this->up_to_vertex_VZ;
      this->BestMatch = 1;
      this->Best3Ddist = this->up_to_vertex_dist3D;

      }
    else if (this->up_clus_michel_dist3D < this->up_to_vertex_dist3D)
    {
      this->BestMatch = 3;
      this->best_XZ = this->up_to_clus_XZ;
      this->best_UZ = this->up_to_clus_VZ;
      this->best_VZ = this->up_to_clus_UZ;
      this->Best3Ddist = this->up_clus_michvtx_dist3D;
      upclusmatch = 1;
    }
   }
   else if (trueEndpoint == 2)
   {
    if (this->down_to_vertex_dist3D < this->down_clus_michel_dist3D)
    {
       this->BestMatch = 2;
       this->best_XZ = this->down_to_vertex_XZ;
       this->best_UZ = this->down_to_vertex_UZ;
       this->best_VZ = this->down_to_vertex_VZ;
       this->Best3Ddist = this->down_to_vertex_dist3D;
       downvtxmatch = 1;
    }
    else if (this->down_to_vertex_dist3D > this->down_clus_michel_dist3D)
    {
      this->BestMatch = 4;
      this->Best3Ddist = this->down_clus_michvtx_dist3D;
      this->best_XZ = this->down_to_clus_XZ;
      this->best_UZ = this->down_to_clus_UZ;
      this->best_VZ = this->down_to_clus_VZ;
      downclusmatch = 1;
    }
   }
   else{
      this->BestMatch = 0;
      this->Best3Ddist = 9999.;
      this->best_XZ = 9999.;
      this->best_UZ = 9999.;
      this->best_VZ = 9999.;
    }
*/
  // This bit of code will try to find the best 3D distance end point
  if (distances[0] == this->up_to_vertex_dist3D) {
    upvtxmatch = 1;
    this->best_XZ = this->up_to_vertex_XZ;
    this->best_UZ = this->up_to_vertex_UZ;
    this->best_VZ = this->up_to_vertex_VZ;
    this->BestMatch = 1;
    this->Best3Ddist = this->up_to_vertex_dist3D;
    // std::cout << "This  Michel is UPVTX and has true endpoint " <<
    // this->trueEndpoint << std::endl;
  } else if (distances[0] == this->down_to_vertex_dist3D) {
    this->BestMatch = 2;
    this->best_XZ = this->down_to_vertex_XZ;
    this->best_UZ = this->down_to_vertex_UZ;
    this->best_VZ = this->down_to_vertex_VZ;
    this->Best3Ddist = this->down_to_vertex_dist3D;
    downvtxmatch = 1;
    // std::cout << "This  Michel is DOWNVTX and has true endpoint " <<
    // this->trueEndpoint << std::endl;

  } else if (distances[0] == this->up_clus_michel_dist3D) {
    this->BestMatch = 3;
    this->best_XZ = this->up_to_clus_XZ;
    this->best_UZ = this->up_to_clus_VZ;
    this->best_VZ = this->up_to_clus_UZ;
    this->Best3Ddist = this->up_clus_michvtx_dist3D;
    upclusmatch = 1;
    //  std::cout << "This  Michel is UPCLUS and has true endpoint " <<
    //  this->trueEndpoint << std::endl;

  } else if (distances[0] == this->down_clus_michel_dist3D) {
    this->BestMatch = 4;
    this->Best3Ddist = this->down_clus_michvtx_dist3D;
    this->best_XZ = this->down_to_clus_XZ;
    this->best_UZ = this->down_to_clus_UZ;
    this->best_VZ = this->down_to_clus_VZ;
    downclusmatch = 1;
    //  std::cout << "This  Michel is DOWNCLUS and has true endpoint " <<
    //  this->trueEndpoint << std::endl;

  } else {
    this->BestMatch = 0;
    this->Best3Ddist = 9999.;
    this->best_XZ = 9999.;
    this->best_UZ = 9999.;
    this->best_VZ = 9999.;
  }

  if (distances[1] == this->up_to_vertex_dist3D)
    this->SecondBestMatch == 1;
  else if (distances[1] == this->down_to_vertex_dist3D)
    this->SecondBestMatch == 2;
  else if (distances[1] == this->up_clus_michel_dist3D)
    this->SecondBestMatch == 3;
  else if (distances[1] == this->down_clus_michel_dist3D)
    this->SecondBestMatch == 4;
  else {
    this->SecondBestMatch == 0;
  }
  /*
  if (upvtxmatch == 1 && (upclusmatch == 1 || downclusmatch == 1 )){
     if (this->up_to_vertex_dist3D < this->up_clus_michel_dist3D) {
        //std::cout << "UPVTX IS BEST " << std::endl;
        this->best_XZ = this->up_to_vertex_XZ;
        this->best_UZ = this->up_to_vertex_UZ;
        this->best_VZ = this->up_to_vertex_VZ;
        this->BestMatch = 1;
        this->Best3Ddist = this->up_to_vertex_dist3D;
     }
     else if (this->up_to_vertex_dist3D >= this->up_clus_michel_dist3D)
     {
        //std::cout << "UPCLUS IS BEST " << std::endl;
        this->BestMatch = 3;
        this->best_XZ = this->up_to_clus_XZ;
        this->best_UZ = this->up_to_clus_VZ;
        this->best_VZ = this->up_to_clus_UZ;
        this->Best3Ddist = this->up_clus_michvtx_dist3D;
     }
     else if (this->up_to_vertex_dist3D < this->down_clus_michel_dist3D) {
        //std::cout << "UPVTX IS BEST " << std::endl;
        this->BestMatch = 1;
        this->Best3Ddist = this->up_to_vertex_dist3D;
        this->best_XZ = this->up_to_vertex_XZ;
        this->best_UZ = this->up_to_vertex_UZ;
        this->best_VZ = this->up_to_vertex_VZ;
     }
     else if (this->up_to_vertex_dist3D >= this->down_clus_michel_dist3D)
     {
        //std::cout << "DOWNCLUSISBEST" << std::endl;
        this->BestMatch = 4;
        this->Best3Ddist = this->down_clus_michvtx_dist3D;
        this->best_XZ = this->down_to_clus_XZ;
        this->best_UZ = this->down_to_clus_UZ;
        this->best_VZ = this->down_to_clus_VZ;
     }
  }
  else if (downvtxmatch == 1 && (upclusmatch == 1 || downclusmatch == 1 )){
    if (this->down_to_vertex_dist3D < this->up_clus_michel_dist3D) {
        //std::cout << "DOWNVTX IS BEST " << std::endl;

        this->BestMatch = 2;
        this->best_XZ = this->down_to_vertex_XZ;
        this->best_UZ = this->down_to_vertex_UZ;
        this->best_VZ = this->down_to_vertex_VZ;
        this->Best3Ddist = this->down_to_vertex_dist3D;

    }
    else if (this->down_to_vertex_dist3D >= this->up_clus_michel_dist3D)
    {
        //std::cout << "UPCLUS IS BEST " << std::endl;
        this->BestMatch = 3;
        this->Best3Ddist = this->up_clus_michvtx_dist3D;
        this->best_XZ = this->up_to_clus_XZ;
        this->best_UZ = this->up_to_clus_UZ;
        this->best_VZ = this->up_to_clus_VZ;
    }
    else if (this->down_to_vertex_dist3D < this->down_clus_michel_dist3D)
    {
        //std::cout << "DOWNVTX IS BEST " << std::endl;

        this->BestMatch = 2;
        this->best_XZ = this->down_to_vertex_XZ;
        this->best_UZ = this->down_to_vertex_UZ;
        this->best_VZ = this->down_to_vertex_VZ;
        //std::cout << "Setting Best Dist" << std::endl;
        this->Best3Ddist = this->down_to_vertex_dist3D;

    }
    else if (this->down_to_vertex_dist3D >= this->down_clus_michel_dist3D)
    {

        //std::cout << "DOWNCLUSISBEST" << std::endl;
        this->BestMatch = 4;
        this->Best3Ddist = this->up_clus_michvtx_dist3D;
        this->best_XZ = this->down_to_clus_XZ;
        this->best_UZ = this->down_to_clus_UZ;
        this->best_VZ = this->down_to_clus_VZ;
    }

  }*/
  // else{
  // this->BestMatch = 0;
  // this->Best3Ddist = 9999.;
  // this->best_XZ = 9999.;
  // this->best_UZ = 9999.;
  // this->best_VZ = 9999.;
  // }
  //    std::cout << "The Best Match for this Michel is " << this->BestMatch <<
  //    std::endl;
  int matchtype = this->BestMatch;
  // Identifying the best reco endpoint based on the Best MAtch type.
  if (matchtype == 1 || matchtype == 3)
    this->recoEndpoint = 1;
  else if (matchtype == 2 || matchtype == 4)
    this->recoEndpoint = 2;
}

void Michel::DoesMichelMatchVtx(const CVUniverse& univ) {
  // std::cout << "GETTING VTX MATCH FOR MICHEL " << std::endl;

  // Getting Vertex Information
  double vtx_x = univ.GetVertex().X();               // mm
  double vtx_y = univ.GetVertex().Y();               // mm
  double vtx_z = univ.GetVertex().Z();               // mm
  double vtx_t = univ.GetVertex().T() / pow(10, 3);  // mus
  double vtx_u = (0.5 * (vtx_x - sqrt(3.) * vtx_y));
  double vtx_v = (0.5 * (vtx_x + sqrt(3.) * vtx_y));

  // std::cout << "VTX POSITION is (x, u , v, y, z) (" << vtx_x << " , " <<
  // vtx_u << " , " << vtx_v << " , " << vtx_y << " , " << vtx_z << std::endl;
  // Initializing all the distance comparisons I will need to make
  double zdiff1 = vtx_z - this->m_z1;
  double zdiff2 = vtx_z - this->m_z2;
  double xdiff = 9999.;
  double udiff = 9999.;
  double vdiff = 9999.;
  double XZdist = 9999.;
  double UZdist = 9999.;
  double VZdist = 9999.;

  double michely1 = this->m_y1;
  double michely2 = this->m_y2;
  double michelx1 = this->m_x1;
  double michelx2 = this->m_z2;
  double michelz1 = this->m_z1;
  double michelz2 = this->m_z2;
  double timediff = (this->time) - vtx_t;
  // std::cout << "Michel time " << this->time << " Vertex time " << vtx_t <<
  // "\n" << std::endl;
  this->vtx_michel_timediff = timediff;

  // 2D distance calculations for Endpoint 1
  xdiff = abs(vtx_x - this->m_x1);
  udiff = abs(vtx_u - this->m_u1);
  vdiff = abs(vtx_v - this->m_v1);

  XZdist = sqrt(xdiff * xdiff + zdiff1 * zdiff1);
  UZdist = sqrt(udiff * udiff + zdiff1 * zdiff1);
  VZdist = sqrt(vdiff * vdiff + zdiff1 * zdiff1);

  this->up_to_vertex_XZ = XZdist;
  this->up_to_vertex_UZ = UZdist;
  this->up_to_vertex_VZ = VZdist;

  // 2D Distance calculations for endpoint2
  xdiff = abs(vtx_x - this->m_x2);
  udiff = abs(vtx_u - this->m_u2);
  vdiff = abs(vtx_v - this->m_v2);
  XZdist = sqrt(xdiff * xdiff + zdiff2 * zdiff2);
  UZdist = sqrt(udiff * udiff + zdiff2 * zdiff2);
  VZdist = sqrt(vdiff * vdiff + zdiff2 * zdiff2);

  this->down_to_vertex_XZ = XZdist;
  this->down_to_vertex_UZ = UZdist;
  this->down_to_vertex_VZ = VZdist;

  // 3D distance calculations
  double xdiff1 = abs(vtx_x - this->m_x1);
  double xdiff2 = abs(vtx_x - this->m_x2);
  double ydiff1 = abs(vtx_y - this->m_y1);
  double ydiff2 = abs(vtx_y - this->m_y2);

  double dist1 = sqrt(zdiff1 * zdiff1 + xdiff1 * xdiff1 + ydiff1 * ydiff1);
  double dist2 = sqrt(zdiff2 * zdiff2 + xdiff2 * xdiff2 + ydiff2 * ydiff2);

  this->up_to_vertex_dist3D = dist1;
  this->down_to_vertex_dist3D = dist2;

  if (dist1 < dist2)
    this->vtx_endpoint = 1;
  else if (dist2 < dist1)
    this->vtx_endpoint = 2;

  ////std::cout << "END OF LOOP TO FIND VTX MATCH" << std::endl;
}

void Michel::DoesMichelMatchClus(const CVUniverse& univ) {
  // This is where the function for Cluster Matching goes

  ////std::cout << "STARTING SEARCH FOR CLUSTER MATCH " << std::endl;
  // Inititalizing vertex variables needed for cluster matching
  int nclusters = this->nclusters;
  // std::cout << "There are " <<  nclusters << "  available clusters for this
  // matching " << std::endl;
  double vtx_x = univ.GetVertex().X();               // mm
  double vtx_y = univ.GetVertex().Y();               // mm
  double vtx_z = univ.GetVertex().Z();               // mm
  double vtx_t = univ.GetVertex().T() / pow(10, 3);  // mus
  double vtx_u = (0.5 * (vtx_x - sqrt(3.) * vtx_y));
  double vtx_v = (0.5 * (vtx_x + sqrt(3.) * vtx_y));

  double closestdistance1x = 9999.;
  double closestdistance1u = 9999.;
  double closestdistance1v = 9999.;
  double closestdistance1z = 9999.;

  double closestdistance2x = 9999.;
  double closestdistance2u = 9999.;
  double closestdistance2v = 9999.;
  double closestdistance2z = 9999.;

  double michelx1 = this->m_x1;
  double michelx2 = this->m_x2;
  double michelu1 = this->m_u1;
  double michelu2 = this->m_u2;
  double michelv1 = this->m_v1;
  double michelv2 = this->m_v2;
  double michelz1 = this->m_z1;
  double michelz2 = this->m_z2;
  double michely1 = this->m_y1;
  double michely2 = this->m_y2;

  double micheltime = this->time;

  // std::cout << "Michel position 1 is (x, u, v, y, z) " << michelx1 << " , "
  // << michelu1 << " , " << michelv1 << " , "<< michely1 <<  " , " << michelz1
  // << std::endl;

  // std::cout << "Michel position 2 is (x, y, v, y, z) " << michelx2 << " , "
  // << michelu2 << " , " << michelv2 << " , " << michely2 << " , " << michelz2
  // << std::endl;

  std::vector<Cluster> endpoint1_clus;
  std::vector<Cluster> endpoint2_clus;

  // Get the closest distance for each view

  ////std::cout << "STARTING LOOP OVER CLUSTERS " << std::endl;
  int x1_idx = -1;  // want to save the index for each closest cluster
  int u1_idx = -1;
  int v1_idx = -1;
  int x2_idx = -1;
  int u2_idx = -1;
  int v2_idx = -1;

  for (int i = 0; i < nclusters; i++) {
    Cluster current_cluster = Cluster(univ, i);

    double energy = current_cluster.energy;
    double time = current_cluster.time;
    double pos = current_cluster.pos;
    double zpos = current_cluster.zpos;
    int view = current_cluster.view;
    double timediff = micheltime - time;

    if (energy < 2.) continue;  // only get clusters greater than 2 MeV

    // std::cout << "printing cluster info " << "energy " << energy << " time "
    // << time << " pos " << pos << " zpos " << zpos << std::endl;

    double zdiff1 = abs(zpos - michelz1);
    double zdiff2 = abs(zpos - michelz2);
    // Calculating 2D distance in X view
    if (view == 1) {
      // Endpoint 1 calculations
      double xdiff1 = abs(pos - michelx1);

      double x2Ddistance1 = sqrt(xdiff1 * xdiff1 + zdiff1 * zdiff1);

      // Endpoint 2 Calculations
      double xdiff2 = abs(pos - michelx2);

      double x2Ddistance2 = sqrt(xdiff2 * xdiff2 + zdiff2 * zdiff2);

      if (x2Ddistance1 <= closestdistance1x) {
        closestdistance1x =
            x2Ddistance1;  // this is redundant if I just use index instead
        x1_idx = i;
      }
      if (x2Ddistance2 <= closestdistance2x) {
        closestdistance2x = x2Ddistance2;
        x2_idx = i;
      }
    } else if (view == 2)  // Calculating 2D distance in U view
    {
      // Endpoint 2 Calculations
      double udiff1 = abs(pos - michelu1);

      double u2Ddistance1 = sqrt(udiff1 * udiff1 + zdiff1 * zdiff1);

      // Endpoint 1 Calculations
      double udiff2 = abs(pos - michelu2);

      double u2Ddistance2 = sqrt(udiff2 * udiff2 + zdiff2 * zdiff2);

      if (u2Ddistance1 < closestdistance1u) {
        closestdistance1u = u2Ddistance1;
        u1_idx = i;
      }
      if (u2Ddistance2 < closestdistance2u) {
        closestdistance2u = u2Ddistance2;
        u2_idx = i;
      }

    } else if (view == 3)  // Calculating 2D dsitance in V view
    {
      // Endpoint 1 Calculations
      double vdiff1 = abs(pos - michelv1);

      double v2Ddistance1 = sqrt(vdiff1 * vdiff1 + zdiff1 * zdiff1);
      // Endpoint 2 Calculations
      double vdiff2 = abs(pos - michelv2);

      double v2Ddistance2 = sqrt(vdiff2 * vdiff2 + zdiff2 * zdiff2);

      if (v2Ddistance1 < closestdistance1v) {
        closestdistance1v = v2Ddistance1;
        v1_idx = i;
      }
      if (v2Ddistance2 <= closestdistance2v) {
        closestdistance2v = v2Ddistance2;
        v2_idx = i;
      }
    }
  }

  // std::cout << "Printing closest clusters index to each end point: x1: " <<
  // x1_idx << " u1: " << u1_idx << " v1: " << v1_idx << " x2: " << x2_idx << "
  // u2: " << u2_idx << " v2: " << v2_idx << std::endl;

  // Now store the closest X, u, v clusters for each Michel Endpoint based on
  // the above closest distance

  // Closest cluster's index will be used to only

  std::vector<double> clusx1;
  std::vector<double> clusx2;

  std::vector<double> clusu1;
  std::vector<double> clusu2;

  std::vector<double> clusv1;
  std::vector<double> clusv2;
  for (int i = 0; i < nclusters; i++) {
    Cluster current_cluster = Cluster(univ, i);

    double energy = current_cluster.energy;
    double time = current_cluster.time;
    double pos = current_cluster.pos;
    double zpos = current_cluster.zpos;
    int view = current_cluster.view;
    double timediff = micheltime - time;
    if (energy < 2.) continue;
    //   std::cout << "Printing details about cluster "<< i << " : "  << energy
    //   << " : " << time << " : " << pos << " : " << zpos << " : " << view << "
    //   : " << timediff << std::endl;

    double zdiff1 = abs(zpos - michelz1);
    double zdiff2 = abs(zpos - michelz2);

    // Saving clusters with distances equal to the closest clusters. Again, this
    // is probably not the best way. I need to rewrite this section to just use
    // clusters from the index I saved in the previous loop. That way I reduce
    // the number of clusters that I have to loop over.

    if (view == 1) {
      // Endpoint 1

      if (i == x1_idx) {  // if current index is same as the closest endpoint 1
                          // cluster in X View
        endpoint1_clus.push_back(current_cluster);
        this->cluster_to_up_match.push_back(
            current_cluster);  // this one is redundant
                               // std::cout << "Printing details about Endpoint
                               // 1 X cluster index "<< i << " energy : "  <<
                               // energy << " time : " << time << " pos : " <<
                               // pos << " zpos : " << zpos << " view : " << view
                               // << " timedifference with Michel : " << timediff
                               // << std::endl;

        clusx1.push_back(pos);
        clusx1.push_back(zpos);
      }
      if (i == x2_idx) {  // if current index is same as the closest endpoint 2
                          // cluster in X View
        endpoint2_clus.push_back(current_cluster);
        this->cluster_to_down_match.push_back(current_cluster);  // TODO remove
        // std::cout << "Printing details about Endpoint 2 X cluster "<< i << "
        // energy : "  << energy << " time : " << time << " pos : " << pos << "
        // zpos: " << zpos << " view : " << view << " timediff : " << timediff <<
        // std::endl;

        clusx2.push_back(pos);
        clusx2.push_back(zpos);
      }

    } else if (view == 2) {
      if (i == u1_idx) {  // if current index is same as the closest endpoint 1
                          // cluster in U View
        endpoint1_clus.push_back(current_cluster);
        this->cluster_to_up_match.push_back(current_cluster);  // Remove?
        // std::cout << "Printing details about Endpoint 1 U  cluster "<< i << "
        // : "  << energy << " : " << time << " : " << pos << " : " << zpos << "
        // : " << view << " : " << timediff << std::endl;

        clusu1.push_back(pos);
        clusu1.push_back(zpos);
      }
      if (i == u2_idx)  // if current index is same as the closest endpoint 2
                        // cluster in U View
      {
        endpoint2_clus.push_back(current_cluster);
        this->cluster_to_down_match.push_back(current_cluster);  // remove?
        // std::cout << "Printing details about Endpoint 2 U cluster "<< i << "
        // : "  << energy << " : " << time << " : " << pos << " : " << zpos << "
        // : " << view << " : " << timediff << std::endl;

        clusu2.push_back(pos);
        clusu2.push_back(zpos);
      }
    } else if (view == 3) {
      if (i == v1_idx) {  // if current index is same as the closest endpoint 1
                          // cluster in V View
        endpoint1_clus.push_back(current_cluster);
        this->cluster_to_up_match.push_back(current_cluster);  // remove?
        // std::cout << "Printing details about Endpoint 1 V cluster "<< i << "
        // : "  << energy << " : " << time << " : " << pos << " : " << zpos << "
        // : " << view << " : " << timediff << std::endl;

        clusv1.push_back(pos);
        clusv1.push_back(zpos);
      }
      if (i == v2_idx) {  // if current index is same as the closest endpoint 2
                          // cluster in V View
        endpoint1_clus.push_back(current_cluster);
        this->cluster_to_down_match.push_back(current_cluster);  // remove?
        // std::cout << "Printing details about Endpoint 2 V cluster "<< i << "
        // : "  << energy << " : " << time << " : " << pos << " : " << zpos << "
        // : " << view << " : " << timediff << std::endl;

        clusv2.push_back(pos);
        clusv2.push_back(zpos);
      }
    }

  }  // End of loop over clusters

  // This is vector of positions for each endpoint cluster match
  std::vector<double> matchclus1;  // index [0] = x, [1] = y, [2] = z
  std::vector<double> matchclus2;  // index [0] = x, [1] = y, [2] = z

  // std::cout << "LOOPING OVER ENDPOINT1 CLUSTERS" << std::endl;

  double XZdist1 = 9999.;
  double UZdist1 = 9999.;
  double VZdist1 = 9999.;

  // Check if our cluster vectors are empty and calculate 2D distances for
  // endpoint 1
  if (!clusx1.empty()) {
    double xdif = abs(this->m_x1 - clusx1[0]);
    double zdif = abs(this->m_z1 - clusx1[1]);
    XZdist1 = sqrt(xdif * xdif + zdif * zdif);
  }
  if (!clusu1.empty()) {
    double udif = abs(this->m_u1 - clusu1[0]);
    double zdif = abs(this->m_z1 - clusu1[1]);
    UZdist1 = sqrt(udif * udif + zdif * zdif);
  }
  if (!clusv1.empty()) {
    double vdif = abs(this->m_v1 - clusv1[0]);
    double zdif = abs(this->m_z1 - clusv1[1]);
    VZdist1 = sqrt(vdif * vdif + zdif * zdif);
  }

  // std::cout << " XZ, UZ, VZ 1: " << XZdist1 << " , " << UZdist1 << " , " <<
  // VZdist1 << std::endl; Saving the 2D distances for endpoint 1
  this->up_to_clus_XZ = XZdist1;
  this->up_to_clus_UZ = UZdist1;
  this->up_to_clus_VZ = VZdist1;

  double XZdist2 = 9999.;
  double UZdist2 = 9999.;
  double VZdist2 = 9999.;

  if (!clusx2.empty()) {
    double xdif = abs(this->m_x2 - clusx2[0]);
    double zdif = abs(this->m_z2 - clusx2[1]);
    XZdist2 = sqrt(xdif * xdif + zdif * zdif);
  }
  if (!clusu2.empty()) {
    double udif = abs(this->m_u2 - clusu2[0]);
    double zdif = abs(this->m_z2 - clusu2[1]);
    UZdist2 = sqrt(udif * udif + zdif * zdif);
  }
  if (!clusv2.empty()) {
    double vdif = abs(this->m_v2 - clusv2[0]);
    double zdif = abs(this->m_z2 - clusv2[1]);
    VZdist2 = sqrt(vdif * vdif + zdif * zdif);
  }

  this->down_to_clus_XZ = XZdist2;
  this->down_to_clus_UZ = UZdist2;
  this->down_to_clus_VZ = VZdist2;
  this->down_clus_y = michely2;
  this->down_clus_x = michelx2;
  this->down_clus_z = michelz2;

  // std::cout << "GET 3D Information for Clusters - ENDPOINT1" << std::endl;
  // This is the convoluted system that calculates the Endpoint
  if (XZdist1 < UZdist1 && XZdist1 < VZdist1 &&
      UZdist1 < VZdist1) {  // XU views closest
    if (!clusu1.empty() && !clusx1.empty()) {
      double yclus = (1. / sqrt(3.)) * (clusx1[0] - 2 * clusu1[0]);
      matchclus1.push_back(clusx1[0]);
      matchclus1.push_back(yclus);      // y point of match 3D point
      matchclus1.push_back(clusx1[1]);  // setting the cluster 3D point z to be
                                        // of the closest view
      // std::cout << "The 2 closest clusters to EndPoint 1 are X and U with (x,
      // y, z) point ( " << clusx1[0] << " , " << yclus << " , " << clusx1[1] <<
      // " ) " << std::endl;
    }
  } else if (XZdist1 < UZdist1 && XZdist1 < VZdist1 &&
             UZdist1 > VZdist1) {  // XV closest
    if (!clusv1.empty() && !clusx1.empty()) {
      double yclus = (1. / sqrt(3.)) * (2 * clusv1[0] - clusx1[0]);
      // std::cout << "The 2 closest clusters to Endpoint 1 are X and V with (x,
      // y, z) point ( " << clusx1[0] << " , " << yclus << " , " << clusx1[1] <<
      // " ) " << std::endl;

      matchclus1.push_back(clusx1[0]);
      matchclus1.push_back(yclus);
      matchclus1.push_back(clusx1[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
    }
  } else if (UZdist1 < XZdist1 && UZdist1 < VZdist1 &&
             VZdist1 < XZdist1) {  // UV closest
    if (!clusv1.empty() && !clusu1.empty()) {
      double yclus = (1. / sqrt(3.)) * (clusv1[0] - clusu1[0]);
      double xclus = clusu1[0] + clusv1[0];
      // std::cout << "The 2 closest clusters to EndPoint 1 are U and V with (x,
      // y, z) point ( " << xclus << " , " << yclus << " , " << clusu1[1] << " )
      // " << std::endl;

      matchclus1.push_back(xclus);
      matchclus1.push_back(yclus);
      matchclus1.push_back(clusu1[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
    }
  } else if (UZdist1 < XZdist1 && UZdist1 < VZdist1 &&
             VZdist1 > XZdist1) {  // UX closest
    if (!clusu1.empty() && !clusx1.empty()) {
      double yclus = (1. / sqrt(3.)) * (clusx1[0] - 2 * clusu1[0]);
      // std::cout << "The 2 closest clusters to EndPoint 1 are U and X with (x,
      // y, z) point ( " << clusx1[0] << " , " << yclus << " , " << clusu1[1] <<
      // " ) " << std::endl;

      matchclus1.push_back(clusx1[0]);
      matchclus1.push_back(yclus);
      matchclus1.push_back(clusu1[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->up_clus_y = michely1;
    }
  } else if (VZdist1 < XZdist1 && VZdist1 < UZdist1 &&
             XZdist1 < UZdist1) {  // VX closest
    if (!clusv1.empty() && !clusx1.empty()) {
      double yclus = ((1. / sqrt(3.)) * (2 * clusv1[0] - clusx1[0]));
      // std::cout << "The 2 closest clusters to EndPoint 1 are V and X with (x,
      // y, z) point ( " << clusx1[0] << " , " << yclus << " , " << clusv1[1] <<
      // " ) " << std::endl;

      matchclus1.push_back(clusx1[0]);
      matchclus1.push_back(yclus);
      matchclus1.push_back(clusv1[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->up_clus_y = michely1;
    }
  } else if (VZdist1 < XZdist1 && VZdist1 < UZdist1 &&
             XZdist1 > UZdist1) {  // VU closest
    if (!clusu1.empty() && !clusv1.empty()) {
      double xclus = (1. / sqrt(3.)) * (clusv1[0] - clusu1[0]);
      double yclus = clusu1[0] + clusv1[0];
      // std::cout << "The 2 closest clusters to EndPoint 1 are V and U with (x,
      // y, z) point ( " << xclus << " , " << yclus << " , " << clusv1[1] << " )
      // " << std::endl;

      matchclus1.push_back(xclus);
      matchclus1.push_back(yclus);
      matchclus1.push_back(clusv1[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->up_clus_y = michely1;
      this->up_clus_x = michelx1;
    }
  }
  if (XZdist2 < UZdist2 && XZdist2 < VZdist2 && UZdist2 < VZdist2) {
    // std::cout << "XU closest to endpoint 2" << std::endl;
    if (!clusu2.empty() && !clusx2.empty()) {
      double yclus = (1. / sqrt(3.)) * (clusx2[0] - 2 * clusu2[0]);
      // std::cout << "The 2 closest clusters to EndPoint 2 are X and U with (x,
      // y, z) point ( " << clusx2[0] << " , " << yclus << " , " << clusx2[1] <<
      // " ) " << std::endl;

      matchclus2.push_back(clusx2[0]);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusx2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view

      this->down_clus_y = michely2;
    }
  } else if (XZdist2 < UZdist2 && XZdist2 < VZdist2 && UZdist2 > VZdist2) {
    // std::cout << "XV closest to endpoint 2" << std::endl;
    if (!clusv2.empty() && !clusx2.empty()) {
      double yclus = (1. / sqrt(3.)) * (2 * clusv2[0] - clusx2[0]);
      // std::cout << "The 2 closest clusters to EndPoint 2 are X and V with (x,
      // y, z) point ( " << clusx2[0] << " , " << yclus << " , " << clusx2[1] <<
      // " ) " << std::endl;

      matchclus2.push_back(clusx2[0]);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusx2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->down_clus_y = michely2;
    }
  } else if (UZdist2 < XZdist2 && UZdist2 < VZdist2 && VZdist2 < XZdist2) {
    // std::cout << "UV closest to endpoint 2" << std::endl;
    if (!clusu2.empty() && !clusv2.empty()) {
      double xclus = clusu2[0] + clusv2[0];
      double yclus = (1. / sqrt(3.)) * (clusv2[0] - clusu2[0]);
      // std::cout << "The 2 closest clusters to EndPoint 2 are U and V with (x,
      // y, z) point ( " << xclus << " , " << yclus << " , " << clusu2[1] << " )
      // " << std::endl;

      matchclus2.push_back(xclus);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusu2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->down_clus_y = michely2;
      this->down_clus_y = michelx2;
    }
  } else if (UZdist2 < XZdist2 && UZdist2 < VZdist2 && VZdist2 > XZdist2) {
    // std::cout << "UX closest to endpoint 2" << std::endl;
    if (!clusu2.empty() && !clusx2.empty()) {
      double yclus = (1. / sqrt(3.)) * (clusx2[0] - 2 * clusu2[0]);
      // std::cout << "The 2 closest clusters to EndPoint 2 are U and X with (x,
      // y, z) point ( " << clusx2[0] << " , " << yclus << " , " << clusu2[1] <<
      // " ) " << std::endl;

      matchclus2.push_back(clusx2[0]);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusu2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->down_clus_y = michely2;
    }
  } else if (VZdist2 < XZdist2 && VZdist2 < UZdist2 && XZdist2 < UZdist2) {
    // std::cout << "VX closest to endpoint 2" << std::endl;
    if (!clusv2.empty() && !clusx2.empty()) {
      double yclus = ((1. / sqrt(3.)) * (2 * clusv2[0] - clusx2[0]));
      // std::cout << "The 2 closest clusters to EndPoint 2 are V and X with (x,
      // y, z) point ( " << clusx2[0] << " , " << yclus << " , " << clusv2[1] <<
      // " ) " << std::endl;

      matchclus2.push_back(clusx2[0]);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusv2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->down_clus_y = michely2;
    }
  } else if (VZdist2 < XZdist2 && VZdist2 < UZdist2 && XZdist2 > UZdist2) {
    // std::cout << "VU closest to endpoint 2" << std::endl;
    if (!clusu2.empty() && !clusv2.empty()) {
      double xclus = (1. / sqrt(3.)) * (clusv2[0] - clusu2[0]);
      double yclus = clusu2[0] + clusv2[0];
      // std::cout << "The 2 closest clusters to EndPoint 2 are V and U with (x,
      // y, z) point ( " << xclus << " , " << yclus << " , " << clusv2[1] << " )
      // " << std::endl;

      matchclus2.push_back(xclus);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusv2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->down_clus_y = michely2;
      this->down_clus_x = michelx2;
    }
  }

  // std::cout << "GETTING 3D DISTANCES FOR THE CUSTERS " << std::endl;
  double clusx1diff = 9999.;
  double clusy1diff = 9999.;
  double clusz1diff = 9999.;

  double mclusx1diff = 9999.;
  double mclusy1diff = 9999.;
  double mclusz1diff = 9999.;

  double michvtx_x1diff = 9999.;
  double michvtx_x2diff = 9999.;
  double michvtx_y1diff = 9999.;
  double michvtx_y2diff = 9999.;
  double michvtx_z1diff = 9999.;
  double michvtx_z2diff = 9999.;

  if (!matchclus1.empty()) {
    clusx1diff = vtx_x - matchclus1[0];
    clusy1diff = vtx_y - matchclus1[1];
    clusz1diff = vtx_z - matchclus1[2];
    mclusx1diff = michelx1 - matchclus1[0];
    mclusy1diff = michely1 - matchclus1[1];
    mclusz1diff = michelz1 - matchclus1[2];
    michvtx_x1diff = michelx1 - vtx_x;
    michvtx_y1diff = michely1 - vtx_y;
    michvtx_z1diff = michelz1 - vtx_z;
    // std::cout << "3D point for Endpoint 1 Clusters is (x, y, z) " <<
    // matchclus1[0] << " , " << matchclus1[1] << " , " << matchclus1[2] <<
    // std::endl;
  }
  double clusx2diff = 9999.;
  double clusy2diff = 9999.;
  double clusz2diff = 9999.;

  double mclusx2diff = 9999.;
  double mclusy2diff = 9999.;
  double mclusz2diff = 9999.;

  if (!matchclus2.empty()) {
    clusx2diff = vtx_x - matchclus2[0];
    clusy2diff = vtx_y - matchclus2[1];
    clusz2diff = vtx_z - matchclus2[2];
    mclusx2diff = michelx2 - matchclus2[0];
    mclusy2diff = michely2 - matchclus2[1];
    mclusz2diff = michelz2 - matchclus2[2];
    michvtx_x2diff = michelx2 - vtx_x;
    michvtx_y2diff = michely2 - vtx_y;
    michvtx_z2diff = michelz2 - vtx_z;
    // std::cout << "3D point for Endpoint 2 Clusters is (x, y, z) " <<
    // matchclus2[0] << " , " << matchclus2[1] << " , " << matchclus2[2] <<
    // std::endl;
  }
  /// 2 types of 3D distance calculations for Cluster matching:
  // 1. Cluster to Vertex
  double dist1 =
      sqrt(pow(clusx1diff, 2) + pow(clusy1diff, 2) + pow(clusz1diff, 2));
  double dist2 =
      sqrt(pow(clusx2diff, 2) + pow(clusy2diff, 2) + pow(clusz2diff, 2));
  // 2. Cluster to Michel
  double mdist1 =
      sqrt(pow(mclusx1diff, 2) + pow(mclusy1diff, 2) + pow(mclusz1diff, 2));
  double mdist2 =
      sqrt(pow(mclusx2diff, 2) + pow(mclusy2diff, 2) + pow(mclusz2diff, 2));
  // std::cout << " The michel endpoint 1 - cluster 3D point distance is " <<
  // mdist1 << std::endl; std::cout << " The michel endpoint 2 - cluster 3D point
  // distance is " << mdist2 << std::endl; Saving all the distances to the Michel
  // member data
  this->down_clus_michel_dist3D = mdist2;
  this->up_clus_michel_dist3D = mdist1;
  this->up_to_cluster_dist3D = dist1;
  this->down_to_cluster_dist3D = dist2;
  // std::cout << "Printing 3D distances to vertex for cluster matches " <<
  // dist1 << " and " << dist2 << std::endl; std::cout << "Printing 3D distances
  // to michel for cluster matches " << mdist1 << " and " << mdist2 << std::endl;
  double michdist1 = sqrt(pow(michvtx_x1diff, 2) + pow(michvtx_y1diff, 2) +
                          pow(michvtx_z1diff, 2));
  double michdist2 = sqrt(pow(michvtx_x2diff, 2) + pow(michvtx_y2diff, 2) +
                          pow(michvtx_z2diff, 2));
  this->up_clus_michvtx_dist3D = michdist1;
  this->down_clus_michvtx_dist3D = michdist2;
  // Marks which endpoint is the closest match for this match type
  if (dist1 < dist2)
    this->clus_endpoint = 1;
  else if (dist1 > dist2)
    this->clus_endpoint = 2;
}

#endif  // MICHEL_H
