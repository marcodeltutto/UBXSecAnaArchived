#! /bin/bash

echo
echo "Running selection on mc_intimecosmic"
echo

./EventSelection --filename files/mcc8.3_v2.0.1/ubxsec_output_mc_intimecosmic_mcc8.3_v2.0.1.root --maxentries -1 --flashstart 3.65 --flashend 5.25 -s 0.406
mv ubxsecana_output.root ubxsecana_output_intimecosmic_mcc8.3_v2.0.1.root


echo 
echo "Running selection on data_extbnb"
echo

./EventSelection --filename files/mcc8.3_v2.0.1/ubxsec_output_data_extbnb_mcc8.3_v2.0.1.root --maxentries -1 --flashstart 3.65 --flashend 5.25 -s 0.406
mv ubxsecana_output.root ubxsecana_output_extbnb_mcc8.3_v2.0.1.root


