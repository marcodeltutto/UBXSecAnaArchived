#! /bin/bash

echo
echo "***************************"
echo "* Running Event Selection"
echo "***************************"
echo

source ./scripts/run_selection.sh

echo
echo "***************************"
echo "* Making plots"
echo "***************************"
echo

source ./scripts/draw_data_mc.sh
