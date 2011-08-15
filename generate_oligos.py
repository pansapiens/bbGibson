import sys, os, random
from optparse import OptionParser
import re

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# NOTE: The SynBio package requires Biopython 1.54
#       I install Biopython from source in a virtualenv
#       The SynBio source is at: https://bitbucket.org/chapmanb/synbio
#       (commit 7b1b3a972b7e is what we have used)
from SynBio.MeltingTemp.Calculate import MeltingTemperature
from SynBio.Primers.Analysis import OligoSetAnalysis
from SynBio.Primers import *

about="""A tool for generating random connector-pair primers,
for use in Biobrick-directed Gibson assembly (bbGibson).
"""
parser = OptionParser(description=about)
parser.add_option("-n", "--number", type="int", dest="NUM_TO_GEN",
                  default="10",
                  help="number of connectors to generate. Actual number of oligos output may be less, since not all can pass requirements")
parser.add_option("-v", "--verbose",
                  action="store_true", dest="VERBOSE", default=False,
                  help="print status messages to stdout")
parser.add_option("-l", "--length", type="int", dest="CONNECTOR_LENGTH", 
                  default="40",
                  help="length of the connector region")
parser.add_option("--min-tm", type="float", dest="MIN_CONNECTOR_TM",
                  default="60.0",
                  help="minimum accepted melting temperature for connector region")
parser.add_option("--max-tm", type="float", dest="MAX_CONNECTOR_TM",
                  default="200.0",
                  help="maximum allowed melting temperature for connector region")
parser.add_option("--hairpin-energy", type="float", dest="LOWEST_HAIRPIN_ENERGY",
                  default="-6.0",
                  help="free energy (delta-G) of most stable predicted hairpin in oligo (UNAFold)")
parser.add_option("-u", "--run-unafold",
                  action="store_true", dest="RUN_UNAFOLD", default=False,
                  help="run UNAFold to predict hairpin free-energy")
parser.add_option("--unafold-bin",
                  action="store", type="string", dest="UNAFOLD_BIN",
                  default="/usr/bin/UNAFold.pl",
                  help="full path to UNAFold.pl executable")
parser.add_option("--max-repeats", type="int", dest="MAX_N_REPEATS", 
                  default="5",
                  help="maximum number of nucleotide repeats allowed (eg AAAAA)")
parser.add_option("-f", "--fasta", dest="OUTPUT_FASTA", action="store_true", 
                  default=False,
                  help="output oligo sequences in FASTA format")

(opts, args) = parser.parse_args()


NUM_TO_GEN = opts.NUM_TO_GEN
#CONNECTOR_LENGTH = 30
CONNECTOR_LENGTH = opts.CONNECTOR_LENGTH
MIN_CONNECTOR_TM = opts.MIN_CONNECTOR_TM
#MAX_CONNECTOR_TM = 90.0 # used if set, ignored if None
MAX_CONNECTOR_TM = opts.MAX_CONNECTOR_TM
LOWEST_HAIRPIN_ENERGY = opts.LOWEST_HAIRPIN_ENERGY
MAX_N_REPEATS = opts.MAX_N_REPEATS # longest allowed single nucleotide repeat (eg AAAAA)

OUTPUT_FASTA = opts.OUTPUT_FASTA
#OUTPUT_FASTA = False


#VERBOSE = True
VERBOSE = opts.VERBOSE

RUN_UNAFOLD = opts.RUN_UNAFOLD
UNAFOLD_BIN = opts.UNAFOLD_BIN
#UNAFOLD_BIN = "/usr/local/bin/UNAFold.pl" #OSX location

# Mutated biobrick prefix:: cAATTCGCGGCCGCTgCTAG - Right hand part - 5' -> 3'
# Mutated biobrick suffix:: TgCTAGTAGCGGCCGCgGCAG - Left hand part - 5' -> 3'
mbbPrefix = "cAATTCGCGGCCGCTgCTAG".lower()
mbbSuffix = "TgCTAGTAGCGGCCGCgGCAG".lower()

# recognition sequences for various restriction sites
restriction_sites = {"EcoRI":"GAATTC", \
           "XbaI":"TCTAGA", \
           "SpeI":"ACTAGT",\
           "PstI":"CTGCAG", \
           "NheI":"GCTAGC", \
           "NotI":"GCGGCCGC", \
           "PvuII":"CAGCTG", \
           "XhoI":"CTCGAG", \
           "AvrII":"CCTAGG", \
           "SapI.1":"GCTCTTC", \
           "SapI.2":"GAAGAGC", \
           "BglII":"AGATCT", \
           "BamHI":"GGATCC", \
           "NgoMIV":"GCCGGC", \
           "AgeI":"ACCGGT" \
          }

# disallowed restriction sites for various assembly standards 
# (these are disallowed in the whole oligo, not just connector)
rfc_sites = {"RFC10":["EcoRI", "XbaI", "SpeI", "PstI", "NotI"], \
       "RFC12":["EcoRI", "NheI", "SpeI", "PstI", "NotI", "PvuII", \
       "XhoI", "AvrII", "XbaI", "SapI.1", "SapI.2"], \
       "RFC21":["EcoRI", "BglII", "BamHI", "XhoI"], \
       "RFC23":["EcoRI", "XbaI", "SpeI", "PstI"], \
       "RFC25":["EcoRI", "XbaI", "NgoMIV", "AgeI", "SpeI", "PstI"]\
      }

allowed_forward_exceptions = {CONNECTOR_LENGTH+6:"NotI"}
allowed_reverse_exceptions = {6:"NotI"}

#Begin helper functions

def err(msg):
  """
  Prints to stderr, only if global VERBOSE is True.
  """
  if VERBOSE: sys.stderr.write(msg)

def GenerateRandomSeq(length):
  seq = ""
  for i in range(length):
    n = random.choice(['A', 'T', 'G', 'C'])
    seq = seq + n
  return seq

def Complement(seq):
  return Seq(seq, IUPAC.unambiguous_dna).complement().tostring()

def ReverseComplement(seq):
  return Seq(seq, IUPAC.unambiguous_dna).complement().tostring()[::-1]


#Checks seq and the reverse complement of seq for illegal restriction sites from all RFC's
def all_restriction_check(seq):
  
  discard = False
  seq = seq.upper()
  seqRC = ReverseComplement(seq)
  
  for rfc in rfc_sites:
    for site in rfc_sites[rfc]:
      
      matches = list(re.finditer(restriction_sites[site], seq))
      
      # (if there are no matches, this loop will never turn)
      for match in matches:
        i = match.start()
        #print "i: {} site: {}" .format(i, site)
        for e in allowed_forward_exceptions:
          if i == e and allowed_forward_exceptions[e] == site:
            err("Found allowed expection site %s::%s (%i) site\n" % (rfc, site, e))
            pass
          else:
            err("Found %s site: %s - discarding primer\n" % (rfc, site))
            discard = True
      
      matches = list(re.finditer(restriction_sites[site], seqRC))
      
      for match in matches:
        i = match.start()
        #print "i: {} site: {}" .format(i, site)
        for e in allowed_reverse_exceptions:
          if i == e and allowed_reverse_exceptions[e] == site:
            err("Found allowed expection site %s::%s (%i) site\n" % (rfc, site, e))
            pass
          else:
            err("Found %s site: %s - discarding primer\n" % (rfc, site))
            discard = True
  
  return discard
  
def unafold_hairpin_check(seq, unafold_bin=UNAFOLD_BIN):
  """
  Run UNAfold to determine delta-G of RNA hairpin.
  Unix only, requires UNAfold installed somewhere 
  (defined with the keyword argument unafold_bin).
  """

  rnd = random.randint(0, 10**8)
  os.system('echo "%s" >/tmp/gibsonator_seq-%i; cd /tmp; %s gibsonator_seq-%i >/dev/null' % (seq, rnd, unafold_bin, rnd))
  f = open("/tmp/gibsonator_seq-%i.dG" % (rnd), 'r')
  f.readline()
  dG_hairpin = float(f.readline().split()[1])
  f.close()

  # cleanup
  os.system("rm /tmp/gibsonator_seq-%i*" % (rnd))

  return dG_hairpin
  
  
#End helpers



# MAIN
#Begin sequence generation

print "\n"

z = 0

for i in range(NUM_TO_GEN):
  
  discard = False
  Sequence = GenerateRandomSeq(CONNECTOR_LENGTH)
  
  ##### Lets setup the sequences here
  
  # Sequence + mutated prefix :: 5' -> 3'
  SeqPrefix = Sequence + mbbPrefix
  # Sequence + mutated Suffix :: 3' -> 5'
  SeqSuffix = Complement(mbbSuffix) + Complement(Sequence)
  
  # Sequence + mutated suffix :: 5' -> 3' ##### Use this as the working sequence!
  SeqSuffixReverse = SeqSuffix[::-1]
  
  seqs = [SeqPrefix, SeqSuffixReverse]
  
  #Look for start codons
  for s in seqs:
    if s.find("ATG") != -1:
      err("Found start codon - discarding primer\n")
      discard = True
  if discard:
    continue
  
  #Restriction enzyme check
  discard = all_restriction_check(SeqPrefix)
  
  if discard:
    continue
  
  #Calculate melting temps for paired connector regions
  melter = MeltingTemperature(0.01, 0.05)
  tm = melter.calculate_melting_plus_align(Sequence, Complement(Sequence))
  #print tm[0]
  
  if tm[0] < MIN_CONNECTOR_TM:
    continue
  
  if MAX_CONNECTOR_TM and tm[0] > MAX_CONNECTOR_TM:
    continue

  #Check connector for repetitive sequences
  repeatPat = re.compile(r"(.)\1{%i,}" % (MAX_N_REPEATS))
  m = re.search(repeatPat, SeqPrefix)
  if m:
    err("Matched sequence of repeats >= 5")
    continue

  if RUN_UNAFOLD:
    hairpin_dG = unafold_hairpin_check(SeqPrefix)
    if hairpin_dG <= LOWEST_HAIRPIN_ENERGY:
      err("RNA Hairpin detected below energy threshold (%.1f): %.1f" % \
                                   (LOWEST_HAIRPIN_ENERGY, hairpin_dG))
      continue

  if not OUTPUT_FASTA:
    spacer = "                     "
    print "Tm of connector region: ", tm[0]
    print "A.r: 3' ", SeqSuffix, " 5'"
    print "B.f: 5' ",spacer + SeqPrefix, " 3'"

    if RUN_UNAFOLD:
      print "dG hairpin (RNA): ", hairpin_dG
    print "\n"

  else:
    if not RUN_UNAFOLD: hairpin_dG=999.0 # dummy number if hairpin_dG not calculated
    print ">%i-A.r [5'-3'] (Tm=%.1f; dG_hairpin=%.1f)\n%s" % \
                  (z, tm[0], hairpin_dG, SeqSuffix[::-1])
    print ">%i-B.f [5'-3'] (Tm=%.1f; dG_hairpin=%.1f)\n%s" % \
                            (z, tm[0], hairpin_dG, SeqPrefix)
  
  z += 1

sys.stderr.write("Sequences passed checks: %s/%s\n" % (z, NUM_TO_GEN))
