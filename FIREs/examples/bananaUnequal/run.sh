#!/usr/bin/env bash

# due to more variables, we can increase -tp with small -lmt or -len

tp=8 # thread pool size

echo M/FIRE -c bananaUnequal -tp $tp
time ../../M/FIRE -c bananaUnequal -tp $tp
echo
