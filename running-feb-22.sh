#!/bin/bash
if [ -f ./csm ]; then
    echo "Clearing CSM executable"
    rm csm
fi
g++ -std=c++11 -o csm -O3 algorithm-feb-22.cpp DFScode-feb-22.cpp replica-feb-22.cpp #replica-dec-2612-NEW.cpp 
#g++ -std=c++11 -o csm -O3 algorithm.cpp DFScode.cpp replica-sep.cpp
#time ./csm mico.txt 12000 output_mico.txt 1 1
#time ./csm mico.txt $1 output_mico_10400.txt 1 1
#time ./csm DBLPorigparsed.txt $1 output_dblp.txt 1 1 
time ./csm lastfm_parsed.txt $1 output_lastfm.txt 1 $2 2>/dev/null 
#time ./csm chemical.txt $1 output_chemical.txt 1 $2
#time ./csm DBLP_coauthor.txt 12500 output_dblp.txt 1 1
#time ./csm YeastUniform.txt 40 output1.txt 1 1
#time ./csm Yeast_Function_Graph.txt 200 output1.txt 1 1
notify-send 'CSM status' 'Program complete'
