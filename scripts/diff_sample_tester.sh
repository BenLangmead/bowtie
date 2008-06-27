#!/bin/sh

#
# Ask the diff_sample.cpp code to produce and sanity check a large
# range of difference samples.
#

make diff_sample-with-asserts
./diff_sample-with-asserts -v -s 3 4096
