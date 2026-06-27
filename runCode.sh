#!/bin/bash

HOMEStella=`pwd`
export HOMEStella

cd $HOMEStella/run/strad || exit 1

# run the code
nohup ./xstella6y12m.exe >& st.log &
echo Stella is started, watch
echo tail -f $HOMEStella/run/strad/st.log

