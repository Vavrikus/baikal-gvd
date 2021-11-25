#!/bin/bash
cd $(dirname $0)

g++ -O3 -o PE.exe pseudo_exp.C `root-config --cflags --libs`