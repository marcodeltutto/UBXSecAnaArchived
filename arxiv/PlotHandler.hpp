#ifndef PlotHandler_h
#define PlotHandler_h

#include <iostream>
#include <iomanip>
#include <string>

#include "Spectrum.hpp"
#include "AnaTree.h"

#include "THStack.h"

using namespace std;

struct CutPlots {

  Spectrum* Snuenergy_numu;
  Spectrum* Snuenergy_anumu;
  Spectrum* Snuenergy_nue;
  Spectrum* Snuenergy_nc;
  Spectrum* Snuenergy_cosmics;
  Spectrum* Snuenergy_outfv;

  Spectrum* Strklen_numu;
  Spectrum* Strklen_anumu;
  Spectrum* Strklen_nue;
  Spectrum* Strklen_nc;
  Spectrum* Strklen_cosmics;
  Spectrum* Strklen_outfv;

  Spectrum* Scostheta_numu;
  Spectrum* Scostheta_anumu;
  Spectrum* Scostheta_nue;
  Spectrum* Scostheta_nc;
  Spectrum* Scostheta_cosmics;
  Spectrum* Scostheta_outfv;
};

class PlotHandler {

  public :
  void Test();
  PlotHandler();
  virtual ~PlotHandler();
  void InstantiateIntermidiatePlots(double totalPOT);
  std::string MakeIntermidiatePlots(std::string cutString, AnaTree * anatree, int bestTrackID);
  THStack* MakeStackHisto(std::string name, std::string label, TLegend* leg, Spectrum *s1, Spectrum *s2, Spectrum *s3, Spectrum *s4, Spectrum *s5, Spectrum *s6);
  void SaveIntermidiatePlots();
  void MakeRatioPlots(std::string,std::string);

  private:
  std::map<std::string,CutPlots> cutToPlotsMap; // This is a map: cutname <-> spectrum array
};

#endif
#ifdef PlotHandler_cxx

PlotHandler::PlotHandler() {


}

PlotHandler::~PlotHandler() {
  
}
#endif // #ifdef PlotHandler_cxx
