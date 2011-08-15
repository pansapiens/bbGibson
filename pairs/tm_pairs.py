#!/usr/bin/env python
#
# Usage: python tm_pairs.py >pairs.dat
#    or: python tm_pairs.py myprimers.fasta >pairs.dat
#
# Tests all-against-all melting temperature for a set of oligonucleotide pairs.
# Takes sequences from a FASTA file, or generates a random set of sequences
# if no file is provided.
#
# Outputs a table suitable for input into R.
# eg:
# # can be plotted in R like:
# > tms<-read.table("pairs.dat", sep="\t", header=TRUE)
# > hist(tms$Tm, xlab="Tm", ylab="number of pairs", main="All-against-all melting temperatures of 1000 primer pairs")

# NOTE: needed to increase the gap opening penalty in SynBio.MeltingTemp.Calculate
#       line 294 shouls be:
#         pairwise2.align.localms(seq_one, seq_two, 5, -1, -100,
#                  -50, one_alignment_only = 1)

import sys, random

from Bio.Seq import Seq
from Bio import SeqIO

# TODO:
#       For pair melting temp determination, consider trying UNAFold ( http://mfold.rna.albany.edu/?q=DINAMelt/software ) and MELTING ( http://www.ebi.ac.uk/compneur-srv/melting/ )

from Bio.Alphabet import IUPAC

from SynBio.MeltingTemp.Calculate import MeltingTemperature
from SynBio.Primers.Analysis import OligoSetAnalysis
from SynBio.Primers import *

dna_conc = 0.01
salt_conc = 0.05
target_tm = 70.0
number_of_connectors = 1000
connector_length = 30
# if true, don't output Tms for primer pairs that have no complementarity
# (Tm = -1.0)
ignore_non_complementary = True

if len(sys.argv) == 1:
  RANDOM = True
else:
  RANDOM = False

def complement(seq):
  return Seq(seq, IUPAC.unambiguous_dna).complement().tostring()

def generate_random_seq(length):
  seq = ""
  for i in range(length):
    n = random.choice(['A', 'T', 'G', 'C'])
    seq = seq + n
  return seq


melter = MeltingTemperature(dna_conc, salt_conc)
#tm = melter.calculate_melting_plus_align(seqX, seqY)
#print tm

connectors = []
if RANDOM:
  # generate random connector sequences
  for i in range(number_of_connectors):
    connectors.append(generate_random_seq(connector_length))
else:
  # read connectors from fasta file, commandline arg
  fh = open(sys.argv[1], 'r')
  for s in SeqIO.parse(fh, "fasta"):
    connectors.append(s.seq.tostring().upper())
  fh.close()

#print connectors

# all-against-all find Tm of connector pairs
tm_pairs = {}
for a in range(len(connectors)):
  for b in range(len(connectors)):
    try:
      tm = melter.calculate_melting_plus_align(connectors[a], complement(connectors[b]))[0]
    except KeyError:
      # occasionally the synbio melting temp code fails on
      # something to do with loop entropy calculation.
      # we just skip those sequences, since it's infrequent &
      # should have any significant effect on overall stats for
      # large sets of primer pairs
      sys.stderr.write( "FAIL :"+connectors[a]+","+complement(connectors[b])+"\n" )
      continue
    #print "%i-%i" % (a, b), tm
    if ignore_non_complementary and tm < 0.0:
      pass
    else:
      tm_pairs["%i-%i" % (a, b)] = tm

# output connector pair Tms for input into R,
# ignore self priming pairs
# can be plotted in R like:
# > tms<-read.table("pairs.dat", sep="\t", header=TRUE)
# > hist(tms$Tm, xlab="Tm", ylab="number of pairs", main="All-against-all melting temperatures of 1000 primer pairs")

print "\tTm"
for pair in tm_pairs:
  pa, pb = pair.split("-")
  if pa == pb:
    continue
  else:
    print "%s\t%.1f" % (pair, tm_pairs[pair])

