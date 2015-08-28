#!/bin/bash

source config.ini
cfile="config.h"

if [ $use64MT19937 -eq 1 ]; then
  sed -i 's/\/\/\ \(#define\ USE_64BIT_MT19937.*$\)/\1/' $cfile
else
  sed -i 's/^.*\(#define\ USE_64BIT_MT19937.*$\)/\/\/\ \1/' $cfile
fi

if [ $debugMode -eq 1 ]; then
  sed -i 's/\/\/\ \(#define\ DEBUG_MODE.*$\)/\1/' $cfile
else
  sed -i 's/^.*\(#define\ DEBUG_MODE.*$\)/\/\/\ \1/' $cfile
fi

exit 0
