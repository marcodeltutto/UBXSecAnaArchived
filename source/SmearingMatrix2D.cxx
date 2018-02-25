#ifndef SMEARINGMATRIX2D_CXX
#define SMEARINGMATRIX2D_CXX

#include "SmearingMatrix2D.h"

namespace ubana {

  SmearingMatrix2D::SmearingMatrix2D()
  {
  }


  void SmearingMatrix2D::SetScaleFactors(double bnbcosmic, double bnbon, double extbnb, double intimecosmic)
  {
  
    _configured = true;
  }

  void SmearingMatrix2D::SetOutputFileName(std::string name) 
  {
    _f_out.open(name, std::ios::out | std::ios::trunc);
  }

  void SmearingMatrix2D::SetTTree(TTree *t)
  {
    _tree = t;
  }

  void SmearingMatrix2D::SetBins(double *var1_b, int n_var1_bins, double *var2_b, int n_var2_bins)
  {

    for (int i = 0; i < n_var1_bins; i++) 
    {
      _var1_bins.push_back(std::make_pair(var1_b[i], var1_b[i+1]));
    }

    for (int i = 0; i < n_var2_bins; i++) 
    {
      _var2_bins.push_back(std::make_pair(var2_b[i], var2_b[i+1]));
    }

    std::cout << "Number of var1 bins: " << _var1_bins.size() << std::endl;
    std::cout << "Number of var2 bins: " << _var2_bins.size() << std::endl;

    _reco_per_true = new TH2D("reco_per_true", "reco_per_true", n_var1_bins, var1_b, n_var2_bins, var2_b);
  }


  void SmearingMatrix2D::CalculateSmearingMatrix() 
  {
    Double_t        mom_tree_true;
    Double_t        mom_tree_mcs;
    //Bool_t          mom_tree_contained;
    //Bool_t          mom_tree_selected;
    Double_t        mom_tree_angle_true;
    Double_t        mom_tree_angle_reco;

    TBranch        *b_mom_tree_true;   //!
    TBranch        *b_mom_tree_mcs;   //!
    //TBranch        *b_mom_tree_contained;   //!
    //TBranch        *b_mom_tree_selected;   //!
    TBranch        *b_mom_tree_angle_true;   //!
    TBranch        *b_mom_tree_angle_reco;   //!

    _tree->SetMakeClass(1);

    _tree->SetBranchAddress("mom_tree_true", &mom_tree_true, &b_mom_tree_true);
    _tree->SetBranchAddress("mom_tree_mcs", &mom_tree_mcs, &b_mom_tree_mcs);
    //_tree->SetBranchAddress("mom_tree_contained", &mom_tree_contained, &b_mom_tree_contained);
    //_tree->SetBranchAddress("mom_tree_selected", &mom_tree_selected, &b_mom_tree_selected);
    _tree->SetBranchAddress("mom_tree_angle_true", &mom_tree_angle_true, &b_mom_tree_angle_true);
    _tree->SetBranchAddress("mom_tree_angle_reco", &mom_tree_angle_reco, &b_mom_tree_angle_reco);

    Long64_t nentries = _tree->GetEntriesFast();

    std::cout << "MMMMMMMMMMMMMMMMMMMMMM Number of entries: " << nentries << std::endl;

    // Resize the smearing matrix
    _S.resize( _var1_bins.size(), 
              std::vector<std::vector<std::vector<double>>> (_var2_bins.size(),
                                                             std::vector<std::vector<double>>(_var1_bins.size(),
                                                                                             std::vector<double> (_var2_bins.size(), 0.
                                                                                                                 )
                                                                                            )
                                                            )
              );

    int counter = 0;

    for (int i = 0; i < _var1_bins.size(); i++) {
      for (int j = 0; j < _var2_bins.size(); j++) {
        for (int m = 0; m < _var1_bins.size(); m++) {
          for (int n = 0; n < _var2_bins.size(); n++) { 
            std::cout << "(" << i << ", " << j << ", " << m << ", " << n << ") => " << _S[i][j][m][n] << std::endl;
            counter++;
          }
        }
      }
    }
    std::cout << "Total entries: " << counter << std::endl;


    // True bin m, n
    int m = 0, n = 0;

    for (int m = 0; m < _var1_bins.size(); m++) {
      for (int n = 0; n < _var2_bins.size(); n++) {

        auto v1_bin = _var1_bins.at(m);
        auto v2_bin = _var2_bins.at(n);

        std::cout << "b1: " << v1_bin.first << " - " << v1_bin.second << std::endl;
        std::cout << "b2: " << v2_bin.first << " - " << v2_bin.second << std::endl;

        _reco_per_true->Reset();

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
          _tree->GetEntry(jentry);

         
        if (  mom_tree_angle_true > v1_bin.first && mom_tree_angle_true < v1_bin.second
           && mom_tree_true > v2_bin.first       && mom_tree_true < v2_bin.second) {

            // Filling reco bin i, j
            _reco_per_true->Fill(mom_tree_angle_reco, mom_tree_mcs);
          }
        }

        // Normalize to get a probability
        _reco_per_true->Scale(1./_reco_per_true->Integral());

        // Set values to matrix
        TCanvas *c = new TCanvas();
        _reco_per_true->Draw("colz text");
        for (int i = 0; i < _var1_bins.size(); i++) {
          for (int j = 0; j < _var2_bins.size(); j++) {
            std::cout << "(" << i << ", " << j << ")" << _reco_per_true->GetBinContent(i+1, j+1) << std::endl;

            double value = _reco_per_true->GetBinContent(i, j);
            if (std::isnan(value))
              value = 0.;
  
            _S[i][j][m][n] = value;
          }
        }

        // Saving the plot
        std::stringstream sstm;
        sstm << "True Bin (" << m << ", " << n << ")";
        std::string str = sstm.str();
        _reco_per_true->SetTitle(str.c_str());
        _reco_per_true->GetXaxis()->SetTitle("Reco Bin i");
        _reco_per_true->GetYaxis()->SetTitle("Reco Bin j");

        sstm.str("");
        sstm << "smearing_matrix_true_" << m << "_" << n << "_";


        TString name = sstm.str();
        c->SaveAs(name + ".pdf");
        c->SaveAs(name + ".C","C");
      }
    }


    for (int i = 0; i < _var1_bins.size(); i++) {
      for (int j = 0; j < _var2_bins.size(); j++) {
        for (int m = 0; m < _var1_bins.size(); m++) {
          for (int n = 0; n < _var2_bins.size(); n++) { 
            std::cout << "(" << i << ", " << j << ", " << m << ", " << n << ") => " << _S[i][j][m][n] << std::endl;
          }
        }
      }
    }

  }


  void SmearingMatrix2D::PrintSmearingMatrixLatex()
  {

    for (int m = 0; m < _var1_bins.size(); m++) {
      for (int n = 0; n < _var2_bins.size(); n++) {
        this->PrintSmearingMatrixLatex(m, n);
      }
    }

  }



  void SmearingMatrix2D::PrintSmearingMatrixLatex(int true_m, int true_n)
  {


    if (!_f_out.is_open()) {
      std::cout << "File not opened." << std::endl;
      return;
    }

    _f_out << "\\begin{equation}" << std::endl;
    _f_out << "S_{ij" << true_m << true_n << "} =" << std::endl;
    _f_out << "\\begin{bmatrix}" << std::endl;

    for (int i = 0; i < _var1_bins.size(); i++) {
      for (int j = 0; j < _var2_bins.size(); j++) {

        _f_out << _S[i][j][true_m][true_n] << "  &  ";

      }

      _f_out << " \\\\" << std::endl;
    }

    _f_out << "\\end{bmatrix}" << std::endl;
    _f_out << "\\end{equation}" << std::endl << std::endl;

  }



}


#endif
