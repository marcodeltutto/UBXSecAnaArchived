#! /bin/bash

echo
echo "***************************"
echo "* Running Event Selection"
echo "***************************"
echo

source run_selection.sh

echo
echo "***************************"
echo "* Making plots"
echo "***************************"
echo

source draw_data_mc.sh
