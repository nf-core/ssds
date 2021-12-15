#!/usr/bin/env python
from collections import defaultdict
from pybedtools import BedTool
import argparse
import os
import sys
import math

## START FUNCTIONS
def log_fr(feature):
    count_fwd = float(feature[-2]) + 1
    count_rev = float(feature[-1]) + 1
    if count_fwd == 1 and count_rev == 1:
        return False

    val = math.log(count_fwd / count_rev, 2)
    feature = feature[:3]
    feature.append(str(val))
    return feature


def sum_fr(feature):
    count_fwd = float(feature[-2])
    count_rev = float(feature[-1])
    if count_fwd == 0 and count_rev == 0:
        return False

    val = count_fwd + count_rev
    feature = feature[:3]
    feature.append(str(val / total_fragments * 1000000))
    return feature


def no_zero(feature):
    count = float(feature[-1])
    if count == 0:
        return False
    else:
        return feature


def rpm(feature):
    count = float(feature[-1])
    if count == 0:
        return False
    else:
        feature[3] = str(count / total_fragments * 1000000)
        return feature


## END FUNCTIONS

parser = argparse.ArgumentParser(description="Make FWD/REV bedgraphs.")

parser.add_argument("--fwd", help="FWD fragments BED")
parser.add_argument("--rev", help="REV fragments BED")
parser.add_argument("--win", help="Windows BED")
parser.add_argument("--name", help="Output file name/suffix")
parser.add_argument("--g", help="Genome index file")
parser.add_argument("--vers", action="store_true", help="Echo version")
parser.add_argument("--v", action="store_true", help="Verbose mode")

args = parser.parse_args()

## "Just show version" mode
if args.vers == True:
    print("2.0.0")
    quit()

fwd_has_reads = True
rev_has_reads = True

total_fragments = 0

if os.path.getsize(args.fwd) == 0:
    fwd_has_reads = False
    print(args.fwd + " is empty ")
else:
    bed_fwd = BedTool(args.fwd).sort(g=args.g)
    total_fragments += bed_fwd.count()

if os.path.getsize(args.rev) == 0:
    rev_has_reads = False
    print(args.rev + " is empty ")
else:
    bed_rev = BedTool(args.rev).sort(g=args.g)
    total_fragments += bed_rev.count()

if args.win:
    wins = BedTool(args.win).sort(g=args.g)

if fwd_has_reads:
    map_fwd = wins.map(b=bed_fwd, c=4, o="count", g=args.g)
    fin_fwd = map_fwd.each(rpm).sort()
    fin_fwd.saveas(args.name + ".FWD.bedgraph")

if rev_has_reads:
    map_rev = wins.map(b=bed_rev, c=4, o="count", g=args.g)
    fin_rev = map_rev.each(rpm).sort()
    fin_rev.saveas(args.name + ".REV.bedgraph")

if fwd_has_reads and rev_has_reads:
    map_both = map_fwd.map(b=bed_rev, c=4, o="count", g=args.g)
    map_log2 = map_both.sort().each(log_fr)
    map_tot = map_both.sort().each(sum_fr)

    map_log2.saveas(args.name + ".FR.bedgraph")
    map_tot.saveas(args.name + ".TOT.bedgraph")
