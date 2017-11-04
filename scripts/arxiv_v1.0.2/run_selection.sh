#! /bin/bash

./EventSelection --filename files/mcc8.3_v1.0.2/ubxsec_output_mc_bnbcosmic_mcc8.3_v1.0.2.root --maxentries -1 --flashstart 3.2 --flashend 4.8 -p
mv ubxsecana_output.root ubxsecana_output_bnbcosmic_mcc8.3_v1.0.2.root

./EventSelection --filename files/mcc8.3_v1.0.2/ubxsec_output_data_bnbon_mcc8.3_v1.0.2.root --maxentries -1 --flashstart 3.3 --flashend 4.9
mv ubxsecana_output.root ubxsecana_output_bnbon_mcc8.3_v1.0.2.root

./EventSelection --filename files/mcc8.3_v1.0.2/ubxsec_output_data_extbnb_mcc8.3_v1.0.2.root --maxentries -1 --flashstart 3.65 --flashend 5.25 -s 0.406
mv ubxsecana_output.root ubxsecana_output_extbnb_mcc8.3_v1.0.2.root


