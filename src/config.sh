#!/bin/bash

source config.ini
cfile="config.h"

if [ $data_visGL_output -eq 1 ]; then
  sed -i 's/\/\/\ \(#define\ DATA_VISGL_OUTPUT.*$\)/\1/' $cfile
else
  sed -i 's/^.*\(#define\ DATA_VISGL_OUTPUT.*$\)/\/\/\ \1/' $cfile
fi

if [ $debugMode -eq 1 ]; then
  sed -i 's/\/\/\ \(#define\ DEBUG_MODE.*$\)/\1/' $cfile
else
  sed -i 's/^.*\(#define\ DEBUG_MODE.*$\)/\/\/\ \1/' $cfile
fi

if [ $verboseMode -eq 1 ]; then
  sed -i 's/\/\/\ \(#define\ VERBOSE_MDOE.*$\)/\1/' $cfile
else
  sed -i 's/^.*\(#define\ VERBOSE_MDOE.*$\)/\/\/\ \1/' $cfile
fi

if [ $use64MT19937 -eq 1 ]; then
  sed -i 's/\/\/\ \(#define\ PRNG_64BIT_MT19937.*$\)/\1/' $cfile
  sed -i 's/^.*\(#define\ PRNG_32BIT_MT19937.*$\)/\/\/\ \1/' $cfile
  sed -i 's/^.*\(#define\ PRNG_DRAND48.*$\)/\/\/\ \1/' $cfile
else
  sed -i 's/^.*\(#define\ PRNG_64BIT_MT19937.*$\)/\/\/\ \1/' $cfile
fi

if [ $use32MT19937 -eq 1 ]; then
  sed -i 's/\/\/\ \(#define\ PRNG_32BIT_MT19937.*$\)/\1/' $cfile
  sed -i 's/^.*\(#define\ PRNG_64BIT_MT19937.*$\)/\/\/\ \1/' $cfile
  sed -i 's/^.*\(#define\ PRNG_DRAND48.*$\)/\/\/\ \1/' $cfile
else
  sed -i 's/^.*\(#define\ PRNG_32BIT_MT19937.*$\)/\/\/\ \1/' $cfile
fi

if [ $useDRAND48 -eq 1 ] || [ $use32MT19937 -eq 0 -a $use64MT19937 -eq 0 ]; then
  sed -i 's/\/\/\ \(#define\ PRNG_DRAND48.*$\)/\1/' $cfile
  sed -i 's/^.*\(#define\ PRNG_32BIT_MT19937.*$\)/\/\/\ \1/' $cfile
  sed -i 's/^.*\(#define\ PRNG_64BIT_MT19937.*$\)/\/\/\ \1/' $cfile
else
  sed -i 's/^.*\(#define\ PRNG_DRAND48.*$\)/\/\/\ \1/' $cfile
fi

exit 0
