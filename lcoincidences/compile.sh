#!/bin/bash
cd $(dirname $0)

g++ -O3 -march=native -o PE.exe pseudo_exp.C `root-config --cflags --libs`