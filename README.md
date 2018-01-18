# UBXSecAna

## Build

```
make (to build all)
make evtsel (to build event selection)
make draw (to build plotting tools)
```

## Run
Update the lines below changing the input file names appropriately. Then run with:

```
MC: ./EventSelection --filename ubxsec_output_mc_bnbcosmic.root --maxentries -1 --flashstart 3.2 --flashend 4.8 -g 198 -p
DATA ON: ./EventSelection --filename ubxsec_output_data_bnbon.root --maxentries -1 --flashstart 3.3 --flashend 4.9 -g 243
DATA OFF: ./EventSelection --filename ubxsec_output_data_extbnb.root --maxentries -1 --flashstart 3.65 --flashend 5.25 -s 0.406 -g 243

Comments: 
-p in MC allows to count POTs
-g value is the dQ/dx calib constant
-s is the shift to apply to flashes
```
