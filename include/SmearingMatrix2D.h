/**
 * \file SmearingMatrix2D.h
 *
 * 
 * \brief Class def header for a class SmearingMatrix2D
 *
 * @author Marco Del Tutto
 */

/** \addtogroup UBXSec

    @{*/
#ifndef SMEARINGMATRIX2D_H
#define SMEARINGMATRIX2D_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <map>
#include <time.h>

#include <TSystem.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TRandom1.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TApplication.h>
#include <TChain.h>
#include "TThread.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraph.h"
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TMath.h>
#include "Math/SMatrix.h"
#include "TMatrix.h"


//#include "PlottingTools.h"

namespace ubana{
  
  /**
   \class SmearingMatrix2D
   User defined class SmearingMatrix2D ... these comments are used to generate
   doxygen documentation!
 */

  class SmearingMatrix2D {
    
  public:
    
    /// Default constructor
    SmearingMatrix2D();

    /// Default destructor
    ~SmearingMatrix2D(){}

    /// Configure function parameters
    void SetScaleFactors(double bnbcosmic, double bnbon, double extbnb, double intimecosmic = 0);

    ///
    void CalculateSmearingMatrix(); 

    ///
    void SetTTree(TTree*);

    ///
    void SetBins(double *, int, double *, int);

    ///
    void SetOutputFileName(std::string name);

    ///
    void PrintSmearingMatrixLatex();

    ///
    void PrintSmearingMatrixLatex(int true_m, int true_n);

  
  protected:

    bool _configured = false;

    TTree *_tree;

    std::vector<std::pair<double, double>> _var1_bins;
    std::vector<std::pair<double, double>> _var2_bins;

    TH2D *_reco_per_true;

    std::vector<std::vector<std::vector<std::vector<double>>>> _S; ///< The smearing matrix

    std::ofstream _f_out; ///< The output file

    
  };
}

#endif
/** @} */ // end of doxygen group 

