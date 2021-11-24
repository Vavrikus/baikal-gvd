#!/bin/bash

g++ -O3 -std=c++2a -o PE.exe pseudo_exp.C `root-config --cflags --libs`