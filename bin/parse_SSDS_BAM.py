#!/usr/bin/env python

from collections import defaultdict
from Bio.Seq import Seq
import pysam
import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Parse ssDNA from SSDS BAM file.")

## Allow "Just show version" mode
parser.add_argument("--bam", help="BAM file to process")
parser.add_argument("--name", help="Output file name/suffix")
parser.add_argument("--mode", help="ss (ss only) or all (5 sub-types)", default="ss")
parser.add_argument("--vers", action="store_true", help="Echo version")
parser.add_argument("--v", action="store_true", help="Verbose mode")

args = parser.parse_args()

if args.vers == True:
    print("2.0.0")
    quit()

## functions ---------------------
def READ_REMOVECLIPPING(read):
    if read.is_unmapped:
        return 1

    newcigar = []
    clip_5 = 0
    clip_3 = 0

    changed = False
    inseq = False

    for op, length in read.cigar:
        if op == 5:  # H
            changed = True
        elif op == 4:  # S
            changed = True
            if not inseq:
                clip_5 = length
            else:
                clip_3 = length
        else:
            inseq = True
            newcigar.append((op, length))

    if not changed:
        return 0

    read.cigar = newcigar
    orig_length = len(read.seq)

    s = read.seq
    q = read.qual

    if clip_3:
        read.seq = s[clip_5:-clip_3]
        if q:
            read.qual = q[clip_5:-clip_3]
    else:
        read.seq = s[clip_5:]
        if q:
            read.qual = q[clip_5:]

    newtags = []
    if clip_5:
        newtags.append(("C5", clip_5))
    if clip_3:
        newtags.append(("C3", clip_3))

    newtags.append(("OL", float(clip_5 + clip_3) / orig_length))

    read.tags = read.tags + newtags

    return 2


def SORT_AND_INDEX(prefix):
    bamid = prefix + ".bam.tmp"
    out = prefix + ".bam"
    pysam.sort(bamid, "-o", prefix + ".bam")
    pysam.index("%s.bam" % prefix)


def PRINT_SSDS_REPORT(prefix, report_dict):
    out_report = open(prefix + ".SSDS_parse_report.txt", "w")

    for field in [
        "ssDNA_fragments",
        "ssLow_fragments",
        "type2_fragments",
        "dsDNA_fragments",
        "unclassified_fragments",
        "total_fragments",
        "valid_fragments",
        "unmapped",
        "filtered_fragments",
        "fragment_too_long",
        "fragment_too_short",
        "alignment_too_short",
        "supplementary_alignment",
    ]:
        out_report.write(
            "totinfo" + "\t" + field + "\t" + str(report_dict["totinfo"][field]) + "\n"
        )

    for type in [
        "ssDNA",
        "ssLow",
        "type2",
        "dsDNA",
        "unclassified",
    ]:
        for rep_value in ["ITR", "uH", "FillIn", "Fragment"]:
            for key in sorted(report_dict[type + "_" + rep_value].keys()):
                out_report.write(
                    "\t".join(
                        [
                            type + "_" + rep_value,
                            str(key),
                            str(report_dict[type + "_" + rep_value][key]) + "\n",
                        ]
                    )
                )

    out_report.close()


## END functions ---------------------

samfile = pysam.AlignmentFile(args.bam, "rb")  # BAM file reader.

names = {
    "ss": args.name + "_ssDNA",
    "sl": args.name + "_ssLow",
    "t2": args.name + "_type2",
    "ds": args.name + "_dsDNA",
    "uc": args.name + "_unclassified",
}

counts = {"ss": 0, "sl": 0, "t2": 0, "ds": 0, "uc": 0}

report_details = {
    "totinfo": {
        "ssDNA_fragments": 0,
        "ssLow_fragments": 0,
        "type2_fragments": 0,
        "dsDNA_fragments": 0,
        "unclassified_fragments": 0,
        "total_fragments": 0,
        "valid_fragments": 0,
        "unmapped": 0,
        "filtered_fragments": 0,
        "fragment_too_long": 0,
        "fragment_too_short": 0,
        "alignment_too_short": 0,
        "supplementary_alignment": 0,
    }
}

for frag_type in [
    "ssDNA",
    "ssLow",
    "type2",
    "dsDNA",
    "unclassified",
]:
    for detail in ["ITR", "uH", "FillIn", "Fragment"]:
        report_details[frag_type + "_" + detail] = {}

out_bam_ssDNA = pysam.AlignmentFile(
    names["ss"] + ".queryname_sorted.bam", "wb", template=samfile
)
out_bam_ssLow = pysam.AlignmentFile(
    names["sl"] + ".queryname_sorted.bam", "wb", template=samfile
)
out_bam_type2 = pysam.AlignmentFile(
    names["t2"] + ".queryname_sorted.bam", "wb", template=samfile
)
out_bam_ds = pysam.AlignmentFile(
    names["ds"] + ".queryname_sorted.bam", "wb", template=samfile
)
out_bam_unc = pysam.AlignmentFile(
    names["uc"] + ".queryname_sorted.bam", "wb", template=samfile
)

out_bed_ssDNA = open(names["ss"] + ".queryname_sorted.bed", "w")
out_bed_ssLow = open(names["sl"] + ".queryname_sorted.bed", "w")
out_bed_type2 = open(names["t2"] + ".queryname_sorted.bed", "w")
out_bed_ds = open(names["ds"] + ".queryname_sorted.bed", "w")
out_bed_unc = open(names["uc"] + ".queryname_sorted.bed", "w")

read1 = None
read2 = None

tot_count = 0

# Iterate through reads.
for read in samfile:

    if read.mate_is_unmapped or not read.is_paired:
        continue

    if read.is_read2:
        read2 = read
    else:
        read1 = read
        read2 = None
        continue

    report_details["totinfo"]["total_fragments"] += 1

    if not read1 is None and not read2 is None and read1.query_name == read2.query_name:

        if read1.is_unmapped or read2.is_unmapped:
            report_details["totinfo"]["unmapped"] += 1
            report_details["totinfo"]["filtered_fragments"] += 1
            continue

        if (abs(read1.reference_start - read2.reference_start)) > 10000:
            report_details["totinfo"]["fragment_too_long"] += 1
            report_details["totinfo"]["filtered_fragments"] += 1
            continue

        if abs(read1.template_length) < 35 or abs(read2.template_length) < 35:
            report_details["totinfo"]["fragment_too_short"] += 1
            report_details["totinfo"]["filtered_fragments"] += 1
            continue

        if abs(read1.reference_length) < 35 or abs(read2.reference_length) < 35:
            report_details["totinfo"]["alignment_too_short"] += 1
            report_details["totinfo"]["filtered_fragments"] += 1
            continue

        if read1.is_supplementary or read2.is_supplementary:
            report_details["totinfo"]["supplementary_alignment"] += 1
            report_details["totinfo"]["filtered_fragments"] += 1
            continue

        if read1.is_reverse:
            if args.v:
                print("R1 rev")
            seq1 = Seq(read1.query_sequence).reverse_complement()
            seq2 = read2.query_sequence
            map = read2.get_reference_positions(full_length=True)
            fragment_strand = "-"
            fragment_start = read2.reference_start
            fragment_end = read1.reference_start + read1.reference_length
        else:
            if args.v:
                print("R1 fwd")
            seq1 = read1.query_sequence
            seq2 = Seq(read2.query_sequence).reverse_complement()
            map = read2.get_reference_positions(full_length=True)[::-1]
            fragment_strand = "+"
            fragment_start = read1.reference_start
            fragment_end = read2.reference_start + read2.reference_length

        fragment_q = str(read1.mapping_quality) + "_" + str(read2.mapping_quality)
        fragment_length = fragment_end - fragment_start

        report_details["totinfo"]["valid_fragments"] += 1
        if args.v:
            print("Valid pair ! ")
            print(seq1)
            print(seq2)

        itr = 0
        fillin = 0
        uhomology = 0
        itr_mismatches = 0
        aln_string = ""

        for i in range(0, len(seq1)):
            if i >= len(seq2):
                break

            if not seq1[i] == seq2[i]:
                itr_mismatches = itr_mismatches + 1

            if itr_mismatches > 1:
                break

            itr = itr + 1

            if map[i] == None:
                aln_string = aln_string + "0"
                fillin = fillin + 1
            else:
                aln_string = aln_string + "1"
                uhomology = uhomology + 1

        if args.v:
            print(aln_string)
            print(str(seq1[0:i]))
            print(
                "ITR = "
                + str(itr)
                + " ; uH = "
                + str(uhomology)
                + "; Fill-in = "
                + str(fillin)
            )
            print()

        code1 = READ_REMOVECLIPPING(read1)
        code2 = READ_REMOVECLIPPING(read2)

        read1.set_tag("it", itr)
        read1.set_tag("mm", itr_mismatches)
        read1.set_tag("uh", uhomology)
        read1.set_tag("fi", fillin)
        read1.set_tag("os", fillin)

        read2.set_tag("it", itr)
        read2.set_tag("mm", itr_mismatches)
        read2.set_tag("uh", uhomology)
        read2.set_tag("fi", fillin)
        read2.set_tag("os", fillin)

        read2.set_tags(
            [
                ("it", itr),
                ("mm", itr_mismatches),
                ("uh", uhomology),
                ("fi", fillin),
                ("os", fillin),
            ]
        )

        if fragment_start >= fragment_end:
            continue

        bed_line = (
            "\t".join(
                [
                    read1.reference_name,
                    str(fragment_start),
                    str(fragment_end),
                    fragment_q,
                    str(itr) + "_" + str(uhomology) + "_" + str(fillin),
                    fragment_strand,
                ]
            )
            + "\n"
        )

        if args.v:
            print(bed_line)

        tot_count += 1
        frag_len = fragment_end - fragment_start

        if itr > 5 and fillin > 2 and uhomology >= 0:
            bam_out = out_bam_ssDNA
            bed_out = out_bed_ssDNA
            frag_type = "ssDNA"
            counts["ss"] = counts["ss"] + 1

        elif fillin > 0 and args.mode == "all":
            bam_out = out_bam_ssLow
            bed_out = out_bed_ssLow
            frag_type = "ssLow"
            counts["sl"] = counts["sl"] + 1

        elif itr >= 4 and fillin == 0 and args.mode == "all":
            bam_out = out_bam_type2
            bed_out = out_bed_type2
            frag_type = "type2"
            counts["t2"] = counts["t2"] + 1

        elif itr < 3 and fillin == 0 and args.mode == "all":
            bam_out = out_bam_ds
            bed_out = out_bed_ds
            frag_type = "dsDNA"
            counts["ds"] = counts["ds"] + 1

        else:
            bam_out = out_bam_unc
            bed_out = out_bed_unc
            frag_type = "unclassified"
            counts["uc"] = counts["uc"] + 1

        bam_out.write(read1)
        bam_out.write(read2)
        bed_out.write(bed_line)

        stat = {"ITR": itr, "FillIn": fillin, "uH": uhomology, "Fragment": frag_len}

        report_details["totinfo"][frag_type + "_fragments"] += 1

        for rep_val in ["ITR", "uH", "FillIn", "Fragment"]:
            if stat[rep_val] in report_details[frag_type + "_" + rep_val]:
                report_details[frag_type + "_" + rep_val][stat[rep_val]] += 1
            else:
                report_details[frag_type + "_" + rep_val][stat[rep_val]] = 1


out_bam_ssDNA.close()
out_bam_ssLow.close()
out_bam_type2.close()
out_bam_ds.close()
out_bam_unc.close()

out_bed_ssDNA.close()
out_bed_ssLow.close()
out_bed_type2.close()
out_bed_ds.close()
out_bed_unc.close()

for key in counts:
    sys.stderr.write(
        "## ITR_id.py: " + key + " : " + str(counts[key]) + " reads found !!\n"
    )
    if counts[key] == 0:
        sys.stderr.write("## ITR_id.py: WARNING: NO " + key + " reads found !!\n")
        os.remove(names[key] + ".queryname_sorted.bam")
        os.remove(names[key] + ".queryname_sorted.bed")

## PRINT FINAL REPORT
PRINT_SSDS_REPORT(args.name, report_details)
