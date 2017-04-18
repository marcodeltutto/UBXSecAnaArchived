#ifndef SelectionTools_h
#define SelectionTools_h

#include <iostream>
#include <iomanip>
#include <string>

#include "AnaTree.h"

using namespace std;

class SelectionTools {

  public :
  int test;
  AnaTree * fanatree;

  void Test();
  SelectionTools(AnaTree*);
  virtual ~SelectionTools();
  bool FlashTag(int &);
  void CreateVertexTrackAssociation(std::map< int,std::vector<int> > &VertexTrackCollection);
  void SelectVertexTrackForward(std::map< int,std::vector<int> > &VertexTrackCollection, int & vertexCandidate);
  bool InFV(std::string type, int vertexCandidate);
  static bool InFV(double x, double y, double z);
  static double GetTrackLength(AnaTree *anatree, int trackIndex);
  int GetBestTrack(int vertexCandidate, std::map< int,std::vector<int> > &VertexTrackCollection); 
  bool IsFlashMatched(int trackCandidate, int theFlash);
  bool IsLongTrack(int trackCandidate);
  int GetEquivalentTrackWithPandoraCosmic(int TrackID);

  // Cut variables
  double flashwidth = 80; //cm. Distance flash-track
  double distcut = 5; //cm. Distance track start/end to vertex
  double lengthcut = 75; //cm. Length of longest track
  double beammin = 3.55/*-0.36*/; //us. Beam window start
  double beammax = 5.15/*-0.36*/; //us. Beam window end
  double PEthresh = 50; //PE
  double MCTrackToMCVtxDist = 0.5; //cm. distance between mc track start and mc vertex
  double TrackToMCDist = 5; //cm. Distance track start/end to mcvertex

  /*This defines our current settings for the fiducial volume
  double FVx = 256.35;
  double FVy = 233;
  double FVz = 1036.8;
  double borderx = 10.;
  double bordery = 20.;
  double borderz = 10.;
  double cryoradius = 191.61;
  double cryoz = 1086.49 + 2*67.63;
  */
};

#endif
#ifdef SelectionTools_cxx

SelectionTools::SelectionTools(AnaTree * anatree) {

  fanatree = anatree;

}

SelectionTools::~SelectionTools() {
  
}
#endif // #ifdef SelectionTools_cxx
