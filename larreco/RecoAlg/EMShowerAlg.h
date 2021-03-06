////////////////////////////////////////////////////////////////////
// Implementation of the EMShower algorithm
//
// Forms EM showers from clusters and associated tracks.
// Also provides methods for finding the vertex and further
// properties of the shower.
//
// Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
////////////////////////////////////////////////////////////////////

#ifndef EMShowerAlg_hxx
#define EMShowerAlg_hxx

// framework
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
namespace fhicl { class ParameterSet; }

// larsoft
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/ShowerEnergyAlg.h"

// C++
#include <map>

// ROOT
#include "RtypesCore.h"
#include "TVector2.h"
#include "TVector3.h"

//temp
#include "larsim/MCCheater/BackTrackerService.h"

class TH1I;
class TProfile;

namespace shower {
  class EMShowerAlg;
  class HitPosition;
}

class shower::EMShowerAlg {
public:

  EMShowerAlg(fhicl::ParameterSet const& pset);

  /// Map associated tracks and clusters together given their associated hits
  void AssociateClustersAndTracks(std::vector<art::Ptr<recob::Cluster> > const& clusters,
				  art::FindManyP<recob::Hit> const& fmh,
				  art::FindManyP<recob::Track> const& fmt,
				  std::map<int,std::vector<int> >& clusterToTracks,
                                  std::map<int,std::vector<int> >& trackToClusters) const;

  /// Map associated tracks and clusters together given their associated hits, whilst ignoring certain clusters
  void AssociateClustersAndTracks(std::vector<art::Ptr<recob::Cluster> > const& clusters,
				  art::FindManyP<recob::Hit> const& fmh,
				  art::FindManyP<recob::Track> const& fmt,
				  std::vector<int> const& clustersToIgnore,
				  std::map<int,std::vector<int> >& clusterToTracks,
                                  std::map<int,std::vector<int> >& trackToClusters) const;

  /// Takes the initial showers found and tries to resolve issues where one bad view ruins the event
  std::vector<int> CheckShowerPlanes(std::vector<std::vector<int> > const& initialShowers,
  				     std::vector<art::Ptr<recob::Cluster> > const& clusters,
                                     art::FindManyP<recob::Hit> const& fmh) const;

  /// Constructs a recob::Track from sets of hits in two views. Intended to be used to construct the initial first part of a shower.
  /// All PMA methods taken from the pma tracking algorithm (R. Sulej and D. Stefan).
  /// This implementation also orients the track in the correct direction if a map of shower centres (units [cm]) in each view is provided.
  std::unique_ptr<recob::Track> ConstructTrack(std::vector<art::Ptr<recob::Hit> > const& track1,
					       std::vector<art::Ptr<recob::Hit> > const& track2,
                                               std::map<int,TVector2> const& showerCentreMap) const;

  /// Constructs a recob::Track from sets of hits in two views. Intended to be used to construct the initial first part of a shower.
  /// All methods taken from the pma tracking algorithm (R. Sulej and D. Stefan).
  std::unique_ptr<recob::Track> ConstructTrack(std::vector<art::Ptr<recob::Hit> > const& track1,
                                               std::vector<art::Ptr<recob::Hit> > const& track2) const;

  /// Finds the initial track-like part of the shower and the hits in all views associated with it
  void FindInitialTrack(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& hits,//art::PtrVector<recob::Hit> const& hits,
			std::unique_ptr<recob::Track>& initialTrack,
                        std::map<int,std::vector<art::Ptr<recob::Hit> > >& initialTrackHits, int plane) const;

  /// Makes showers given a map between tracks and all clusters associated with them
  std::vector<std::vector<int> > FindShowers(std::map<int,std::vector<int> > const& trackToClusters) const;

  /// Makes a recob::Shower object given the hits in the shower and the initial track-like part
  recob::Shower MakeShower(art::PtrVector<recob::Hit> const& hits,
			   std::unique_ptr<recob::Track> const& initialTrack,
                           std::map<int,std::vector<art::Ptr<recob::Hit> > > const& initialTrackHits) const;

  /// <Tingjun to document>
  recob::Shower MakeShower(art::PtrVector<recob::Hit> const& hits,
			   art::Ptr<recob::Vertex> const& vertex,
                           int & iok) const;

  /// Makes space points from the shower hits in each plane
  std::vector<recob::SpacePoint> MakeSpacePoints(std::map<int,std::vector<art::Ptr<recob::Hit> > > hits, std::vector<std::vector<art::Ptr<recob::Hit> > >& hitAssns) const;

  /// Takes the hits associated with a shower and orders them so they follow the direction of the shower
  //std::map<int,std::vector<art::Ptr<recob::Hit> > > OrderShowerHits(art::PtrVector<recob::Hit> const& shower, int plane);
  std::map<int,std::vector<art::Ptr<recob::Hit> > > OrderShowerHits(const art::PtrVector<recob::Hit>& shower, int plane) const;

  /// <Tingjun to document>
  void FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >const& showerHits,
			    art::Ptr<recob::Vertex> const& vertex,
                            std::vector<art::Ptr<recob::Hit> >& trackHits) const;

  /// <Tingjun to document>
  Int_t WeightedFit(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,  Double_t *parm) const;

  /// <Tingjun to document>
  bool isCleanShower(std::vector<art::Ptr<recob::Hit> > const& hits) const;

  int fDebug;

private:

  /// Checks the hits across the views in a given shower to determine if there is one in the incorrect TPC
  void CheckIsolatedHits(std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap) const;

  /// Takes the shower hits in all views and ensure the ordering is consistent
  /// Returns bool, indicating whether or not everything makes sense!
  bool CheckShowerHits(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& showerHitsMap) const;

  /// Constructs a 3D point (in [cm]) to represent the hits given in two views
  TVector3 Construct3DPoint(art::Ptr<recob::Hit> const& hit1, art::Ptr<recob::Hit> const& hit2) const;

  /// Finds dE/dx for the track given a set of hits
  double FinddEdx(std::vector<art::Ptr<recob::Hit> > const& trackHits, std::unique_ptr<recob::Track> const& track) const;

  /// Orders hits along the best fit line through the charge-weighted centre of the hits.
  /// Orders along the line perpendicular to the least squares line if perpendicular is set to true.
  std::vector<art::Ptr<recob::Hit> > FindOrderOfHits(std::vector<art::Ptr<recob::Hit> > const& hits, bool perpendicular = false) const;

  /// Takes a map of the shower hits on each plane (ordered from what has been decided to be the start)
  /// Returns a map of the initial track-like part of the shower on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > FindShowerStart(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& orderedShowerMap) const;

  /// Takes all the shower hits, ready ordered, and returns information to help with the orientation of the shower in each view
  /// Returns map of most likely permutations of reorientation
  /// Starts at 0,0,0 (i.e. don't need to reorient any plane) and ends with 1,1,1 (i.e. every plane needs reorienting)
  /// Every permutation inbetween represent increasing less likely changes to satisfy the correct orientation criteria
  std::map<int,std::map<int,bool> > GetPlanePermutations(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap) const;

  /// Find the global wire position
  double GlobalWire(const geo::WireID& wireID) const;

  /// Return the coordinates of this hit in global wire/tick space
  TVector2 HitCoordinates(art::Ptr<recob::Hit> const& hit) const;

  /// Return the coordinates of this hit in units of cm
  TVector2 HitPosition(art::Ptr<recob::Hit> const& hit) const;

  /// Return the coordinates of this hit in units of cm
  TVector2 HitPosition(TVector2 const& pos, geo::PlaneID planeID) const;

  /// Takes initial track hits from multiple views and forms a track object which best represents the start of the shower
  std::unique_ptr<recob::Track> MakeInitialTrack(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& initialHitsMap,
                                                 std::map<int,std::vector<art::Ptr<recob::Hit> > > const& showerHitsMap) const;

  /// Takes the hits associated with a shower and orders then so they follow the direction of the shower
  void OrderShowerHits(std::vector<art::Ptr<recob::Hit> > const& shower,
		       std::vector<art::Ptr<recob::Hit> >& orderedShower,
                       art::Ptr<recob::Vertex> const& vertex) const;

  /// Projects a 3D point (units [cm]) onto a 2D plane
  /// Returns 2D point (units [cm])
  TVector2 Project3DPointOntoPlane(TVector3 const& point, int plane, int cryostat = 0) const;

  /// Determines the 'relative wire width', i.e. how spread a shower is across wires of each plane relative to the others
  /// If a shower travels along the wire directions in a particular view, it will have a smaller wire width in that view
  /// Returns a map relating these widths to each plane
  std::map<double,int> RelativeWireWidth(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap) const;

  /// Returns the charge-weighted shower centre
  TVector2 ShowerCentre(const std::vector<art::Ptr<recob::Hit> >& showerHits) const;

  /// Returns a rough charge-weighted shower 'direction' given the hits in the shower
  TVector2 ShowerDirection(const std::vector<art::Ptr<recob::Hit> >& showerHits) const;

  /// Returns the RMS of the hits from the central shower 'axis' along the length of the shower
  double ShowerHitRMS(const std::vector<art::Ptr<recob::Hit> >& showerHits) const;

  /// Returns the gradient of the RMS vs shower segment graph
  double ShowerHitRMSGradient(const std::vector<art::Ptr<recob::Hit> >& showerHits, TVector2 trueStart = TVector2(0,0)) const;

  /// Returns the plane which is determined to be the least likely to be correct
  int WorstPlane(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap) const;


  // Parameters
  double fMinTrackLength;
  double fdEdxTrackLength;
  double fSpacePointSize;
  // Parameters to fit wire vs time
  unsigned int         fNfitpass;
  std::vector<unsigned int>     fNfithits;
  std::vector<double>  fToler;

  // Services used by this class
  art::ServiceHandle<geo::Geometry const> fGeom;
  detinfo::DetectorProperties const* fDetProp;
  art::ServiceHandle<art::TFileService const> tfs;

  // Algs used by this class
  shower::ShowerEnergyAlg fShowerEnergyAlg;
  calo::CalorimetryAlg fCalorimetryAlg;
  pma::ProjectionMatchingAlg fProjectionMatchingAlg;

  std::string fDetector;


  // tmp
  int FindTrueParticle(const std::vector<art::Ptr<recob::Hit> >& showerHits) const;
  int FindParticleID(const art::Ptr<recob::Hit>& hit) const;
  art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
  TH1I* hTrueDirection;
  TProfile* hNumHitsInSegment, *hNumSegments;
  //  void MakePicture();
  bool fMakeGradientPlot, fMakeRMSGradientPlot;
  int fNumShowerSegments;

};

class shower::HitPosition {
 public:

  TVector2 WireTick;
  TVector2 Cm;

  // default constructor
  HitPosition();
  // contructors
  HitPosition(TVector2 wiretick, TVector2 cm): HitPosition()
  {
    WireTick = wiretick;
    Cm = cm;
  }
  HitPosition(TVector2 wiretick, geo::PlaneID planeID): HitPosition() {
    WireTick = wiretick;
    Cm = ConvertWireTickToCm(wiretick, planeID);
  }

  TVector2 ConvertWireTickToCm(TVector2 wiretick, geo::PlaneID planeID) {
    return TVector2(wiretick.X() * fGeom->WirePitch(planeID),
		    fDetProp->ConvertTicksToX(wiretick.Y(), planeID));
  }

 private:

  geo::GeometryCore const* fGeom = nullptr;
  detinfo::DetectorProperties const* fDetProp = nullptr;

};

#endif
