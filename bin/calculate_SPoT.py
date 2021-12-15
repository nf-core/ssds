#!/usr/bin/env python
from collections import defaultdict
from pybedtools import BedTool
import argparse
import os
import re
import sys
import csv
import math

parser = argparse.ArgumentParser(
    description="Calculate the signal portion of tags / fragments (SPoT)."
)

## Allow "Just show version" mode
parser.add_argument("--reads_bed", help="Reads/Fragments BED")
parser.add_argument("--reads_bam", help="Reads BAM")
parser.add_argument("--intervals_bed", help="Genomic intervals BED")
parser.add_argument("--name", help="Reads name")
parser.add_argument("--iname", help="Intervals name")
parser.add_argument("--g", help="Genome index file")
parser.add_argument("--o", default="no_name", help="Output prefix")
parser.add_argument(
    "--rand", action="store_true", help="Also show SPoT for a randomized set"
)
parser.add_argument("--vers", action="store_true", help="Echo version")
parser.add_argument("--v", action="store_true", help="Verbose mode")

args = parser.parse_args()

if args.vers == True:
    print("1.0.0")
    quit()

out_report = open(args.o + ".SSDS_SPoT_report.txt", "w")


def calculate_spot(bed, bam, reads_id):
    counts = {}

    ## Allow checking of lots of interval BEDs at once
    if args.intervals_bed == "all":
        intervals_bigbed = open("allintervals.tab", "w")

        directory = r"."
        for entry in os.listdir("."):
            if entry.endswith(".bed"):
                #            if entry.path.endswith(".bed") and entry.is_file():
                id = re.sub("\.bed$", "", entry)
                counts[id] = 0
                if args.rand:
                    counts[id + "(R)"] = 0
                with open(entry) as tsv:
                    for l in csv.reader(tsv, delimiter="\t"):
                        intervals_bigbed.write(
                            "\t".join([l[0], l[1], l[2], id, id, "+\n"])
                        )

        intervals_bigbed.close()
        intervals_bed = BedTool("allintervals.tab")
        intervals_name = "all"
    else:
        counts[re.sub("\.bed$", "", args.intervals_bed)] = 0
        if args.rand:
            counts[re.sub("\.bed$", "", args.intervals_bed) + "(R)"] = 0
        intervals_bed = BedTool(args.intervals_bed)
        intervals_name = args.iname

    rand_intervals = intervals_bed.shuffle(g=args.g, chrom=True, seed=42)

    if bed is not None:
        reads = BedTool(bed)
        intervals_bed = intervals_bed.sort()
        rand_intervals_bed = rand_intervals.sort()
    elif bam is not None:
        reads = BedTool(bam)
        intervals_bed = intervals_bed.sort(g=args.g)
        rand_intervals_bed = rand_intervals.sort(g=args.g)
    else:
        sys.exit("Reads BED / BAM not provided")

    at_intervals = intervals_bed.intersect(
        b=reads, c=True, sorted=True, stream=True
    ).sort()

    for i in at_intervals:
        counts[i.name] = counts[i.name] + i.count

    at_random = rand_intervals_bed.intersect(b=reads, c=True, sorted=True, stream=True)

    for i in at_random:
        counts[i.name + "(R)"] = counts[i.name + "(R)"] + i.count

    num_reads = float(reads.count())

    for c in counts:
        if counts[c] > 0:
            SPoT = str(float(counts[c]) / num_reads)
        else:
            SPoT = 0

        report_line = "\t".join([reads_id + "_SPoT", c, str(SPoT)])

        print(report_line)

        out_report.write(report_line + "\n")


bams = re.sub("[\[\]\s]", "", args.reads_bam).split(",")
names = re.sub("[\[\]\s]", "", args.name).split(",")

for i in range(0, len(bams)):
    calculate_spot(None, bams[i], names[i])
