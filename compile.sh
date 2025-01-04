#!/bin/bash

g++ -c source/routines.cpp -Iinclude
g++ -o main.exe main.cpp routines.o -Iinclude

