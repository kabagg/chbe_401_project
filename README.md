README
================
Keith Baggerly
12/5/2019

This repository houses files I used in preparing my summary of
geothermal loops for fluid mechanics.

There are two R script files to be run in sequence:

01\_load\_parameter\_values.R - this file collects various properties of
ethane and ammonia, such as critical temperatures and pressures, and
values which were read by eye from Eavor’s PH diagram.

02\_complete\_vertex\_values.R - this uses isobaric data at key pressure
levels from NIST to infer other properties (e.g., temperatures) of the
vertices on Eavor’s PH diagram using linear interpolation. This file
actually sources the first, so running this shouldn’t break things.

data/

data files from NIST:

  - c2h6\_0p32.tsv
  - c2h6\_0p78.tsv
  - c2h6\_0p94.tsv

isobaric properties of ethane (C\_2H\_6) evaluated at reduced pressures
of 0.32, 0.78, and 0.94, in unit Kelvin increments across the range NIST
provides

  - nh3\_0p32.tsv
  - nh3\_0p78.tsv
  - nh3\_0p94.tsv

isobaric properties of ammonia (NH\_3) evaluated at reduced pressures of
0.32, 0.78, and 0.94, in unit Kelvin increments across the range NIST
provides

other files:

  - vertex\_values\_used.csv

complete table of values both read off and inferred for the vertices in
Eavor’s PH diagram
