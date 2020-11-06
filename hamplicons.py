#!/usr/bin/env python
# FIXME:
# pylint: disable=line-too-long
import collections
import json
import logging
import multiprocessing
import os
import shutil
import sys
import datetime

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from gzip import open as gzopen
from itertools import combinations
from subprocess import Popen

import coloredlogs
import xlsxwriter


__version_tuple__ = (4, 0, 0)
__version__ = ".".join(map(str, __version_tuple__))


## Input constants
# Expected format/extension of input files
READ_FORMAT = ".fastq.gz"
# Minimum read length, used to exclude primer dimers, etc.
MIN_READ_LENGTH = 50

# minimum number of reads a target must have to be reported
TARGET_TOTAL_THRESHOLD = 100
# minimum number of reads a specific target variation must have
TARGET_LOCAL_THRESHOLD_READS = 10
# minimum percentage of reads a specific indel must contribute to the total indels to be reported
TARGET_LOCAL_THRESHOLD_PERCENTAGE = 0.05

## Unmatched reporting constants
# report significant products with >=this reads
SIGNIFICANT_UNMATCHED_THRESHOLD_READS = 500
# require these products to be at least this long (otherwise its probably primer-dimers)
SIGNIFICANT_UNMATCHED_THRESHOLD_LENGTH = 50

## Barcode matching constants
# bp in start/end that must match to categorize fastq reads as specific product from fasta
BARSIZE = 30
# when performing additional matching
REDUCED_BAR = 15

## Other constants
# Size of regions compared between PCR products to ensure we dont have way too similar sequences
WORRY_ZONE = 40


TargetSequence = collections.namedtuple(
    "TargetSequence",
    ["name", "seq", "len", "barcode_5p", "barcode_3p", "ham_5p", "ham_3p"],
)


def run_command(command, **kwargs):
    proc = Popen(command, **kwargs)
    return proc.wait()


def setup_logging(filename, level=logging.INFO):
    fmt = "%(asctime)s - %(message)s"
    if not sys.stdout.isatty():
        fmt = "%(asctime)s - %(levelname)s - %(message)s"

    coloredlogs.install(
        level=level,
        fmt=fmt,
        datefmt="%H:%M:%S",
        level_styles={
            "info": {"color": "white"},
            "critical": {"color": "red", "bold": True},
            "verbose": {"color": "blue"},
            "error": {"color": "red"},
            "debug": {"color": "green"},
            "warning": {"color": "yellow"},
        },
    )

    # setup a global warning logger
    formatter = logging.Formatter("%(asctime)s : %(levelname)s : %(message)s")
    handler = logging.FileHandler(filename, mode="w")
    handler.setFormatter(formatter)
    handler.setLevel(level)
    logger = logging.getLogger("")
    logger.addHandler(handler)


def convert_index_to_well(index, rows=8, columns=12):  # 1 based index
    thisrow = (index - 1) / columns
    thisplate = thisrow / rows
    thiscolumn = index - thisrow * columns
    thisrow -= rows * thisplate
    return "{}{}_p{}".format(chr(65 + thisrow), thiscolumn, thisplate)


def hamming(seq1, seq2):
    return sum(nuc1 != nuc2 for nuc1, nuc2 in zip(seq1, seq2))


def run_flash(args):
    timestamp = datetime.datetime.now()
    destination, sample = args

    input_1, input_2 = sample["R1"], sample["R2"]
    output = os.path.join(destination, "index%03i" % (sample["#"],))
    logfile = output + ".stdout"

    command = ["flash", "-t", "1", "-M", "151", "-z", "-o", output, input_1, input_2]

    log = logging.getLogger("flash")
    log.debug("Sample %03i: Merging reads into %s.*", sample["#"], output)
    with open(logfile, "w") as loghandle:
        if run_command(command, stdout=loghandle):
            log.error("error while running FLASH; see log at %r", logfile)
            return False

    delta = (datetime.datetime.now() - timestamp).total_seconds()
    log.info("Sample %03i merged in %.1f seconds", sample["#"], delta)

    return True


def merge_fastq(destination, samples, threads):
    log = logging.getLogger("flash")
    log.info("Merging FASTQ files")

    args = [(destination, sample) for sample in samples]
    pool = multiprocessing.Pool(processes=threads)
    for result in pool.imap(run_flash, args):
        if not result:
            pool.terminate()
            return False

    return True


def count_fastq_sequences(filename, min_length=MIN_READ_LENGTH):
    counts = collections.defaultdict(int)

    with gzopen(filename, mode="rt") as handle:
        for idx, line in enumerate(handle):
            if idx % 4 == 1:
                read = line.strip().upper()
                if len(read) > min_length:
                    counts[read] += 1

    return dict(counts)


def collect_fastq_files(log, root, recursive):
    log.info("Locating FASTQ files in %s", root)

    fastq_files = collections.defaultdict(dict)
    for (dirpath, _, filenames) in os.walk(root):
        for filename in filenames:
            if filename.endswith(READ_FORMAT):
                # someName_Si_L001_Rx_001.fastq.gz
                fields = filename.split("_")
                name = fields[0]
                sample = fields[-4]  # x = sequential sample numbering prefixed with 'S'
                pair = fields[-2]  # i = R1 or R2 (paired end)

                pairs = fastq_files[(int(sample[1:]), name)]
                if pair in pairs:
                    log.error("Multiple %s files for %s (%s)", pair, name, sample)
                    return None
                elif pair not in ("R1", "R2"):
                    log.error("Unexpected %s file for %s (%s)", pair, name, sample)
                    return None

                pairs[pair] = os.path.join(dirpath, filename)

        if not recursive:
            break

    samples = []
    for idx, ((sample, name), pairs) in enumerate(sorted(fastq_files.items()), start=1):
        if len(pairs) != 2:
            log.error("Missing files for sample %r (%i)", sample, name)
            return None

        samples.append({"#": idx, "name": name, "R1": pairs["R1"], "R2": pairs["R2"]})

    log.info("Identified %g FASTQ file pairs", len(samples))

    return samples


##############################################################################
# Reading, checking, and post-processing of FASTA target sequences


def read_fasta(fasta_file):
    name = None
    sequence = []
    records = []

    with open(fasta_file) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence)))

                # Leave out header meta-data, if present
                name = line[1:].split(None, 1)[0]
                sequence = []
            else:
                sequence.append(line.strip())

        if name is not None:
            records.append((name, "".join(sequence)))

    return records


def merge_identical_targets(targets):
    log = logging.getLogger("fasta")

    targets_by_barcodes = collections.defaultdict(list)
    for target in targets:
        targets_by_barcodes[(target.barcode_5p, target.barcode_3p)].append(target)

    filtered_targets = []
    for group in targets_by_barcodes.values():
        if len(group) > 1:
            if any(group[0].seq != target.seq for target in group[1:]):
                log.error(
                    "Found targets with identical barcodes but different sequences:"
                )
                for target in group:
                    log.error("    - %s", target.name)
                return None

            log.critical("Merging targets with identical sequence:")
            for target in group:
                log.critical("    - %s", target.name)

            name = "__OR__".join(sorted(set(target.name for target in group)))
            group = [group[0]._replace(name=name)]

        filtered_targets.extend(group)

    return filtered_targets


def compare_targets_5p(target1, target2, pos, max_ham, log):
    ham = hamming(target1.barcode_5p, target2.seq[pos : pos + BARSIZE])
    if ham < max_ham:
        log.warning(
            "Barcode from %r appears in %r, ham = %s", target1.name, target2.name, ham
        )
        log.warning(
            "    %s%s%s...",
            "." * pos,
            target1.barcode_5p,
            "." * (WORRY_ZONE - BARSIZE - pos),
        )
        log.warning("    %s...", target2.seq[:WORRY_ZONE])


def compare_targets_3p(target1, target2, pos, max_ham, log):
    offset = len(target2.seq) - WORRY_ZONE + pos
    ham = hamming(target1.barcode_3p, target2.seq[offset : offset + BARSIZE])
    if ham < max_ham:
        log.warning(
            "Barcode from %r appears in %r, ham = %s", target1.name, target2.name, ham
        )
        log.warning(
            "    ...%s%s%s",
            "." * pos,
            target1.barcode_3p,
            "." * (WORRY_ZONE - BARSIZE - pos),
        )
        log.warning("    ...%s", target2.seq[-WORRY_ZONE:])


def check_targets(targets, max_ham):
    # 1. Check for identical sequences and sequences with identical barcodes
    targets = merge_identical_targets(targets)
    if targets is None:
        return

    # 2. Check for FASTA records with the same name but different sequences
    log = logging.getLogger("fasta")
    targets_by_name = {}
    for target in targets:
        if target.name in targets_by_name:
            log.error(
                "Multiple, different target sequences with the name %r", target.name
            )
            return

        targets_by_name[target.name] = target

    # 2: itereate each barcode through worry zone
    for key1, key2 in combinations(targets_by_name, 2):
        target1 = targets_by_name[key1]
        target2 = targets_by_name[key2]

        for pos in range(WORRY_ZONE - BARSIZE + 1):
            compare_targets_5p(target1, target2, pos, max_ham, log)
            compare_targets_5p(target2, target1, pos, max_ham, log)
            compare_targets_3p(target1, target2, pos, max_ham, log)
            compare_targets_3p(target2, target1, pos, max_ham, log)

        # get nearest ham
        ham_5p = hamming(target1.barcode_5p, target2.barcode_5p)
        ham_3p = hamming(target1.barcode_3p, target2.barcode_3p)

        # set nearest ham
        targets_by_name[target1.name] = target1._replace(
            ham_5p=min(target1.ham_5p, ham_5p), ham_3p=min(target1.ham_3p, ham_3p)
        )
        targets_by_name[target2.name] = target2._replace(
            ham_5p=min(target2.ham_5p, ham_5p), ham_3p=min(target2.ham_3p, ham_3p)
        )

    targets = sorted(targets_by_name.values(), key=lambda item: item.name.lower())

    name_len = max(len(x.name) for x in targets)
    log.info(
        "{:<{}} {:>8} {:>8} {:>8}".format(
            "Name", name_len, "Length", "5p Ham", "3p Ham"
        )
    )
    for target in targets:
        log.info(
            "{:<{}} {:>8} {:>8} {:>8}".format(
                target.name, name_len, target.len, target.ham_5p, target.ham_3p
            )
        )

    return targets


def build_targets(fasta, max_ham):
    targets = []
    for (name, seq) in fasta:
        seq = seq.upper()

        targets.append(
            TargetSequence(
                name=name,
                seq=seq,
                len=len(seq),
                barcode_5p=seq[:BARSIZE],
                barcode_3p=seq[-BARSIZE:],
                ham_5p=max_ham,
                ham_3p=max_ham,
            )
        )

    return check_targets(targets, max_ham)


##############################################################################


def strict_matching(targets, read_counts, results, forlog):
    for target in targets:
        for read, count in tuple(read_counts.items()):
            if read.startswith(target.barcode_5p) and read.endswith(target.barcode_3p):
                results[target.name][len(read) - target.len] += count
                read_counts.pop(read)  # remove analyzed reads

                forlog[target.name][0] += count


def relaxed_matching(args, targets, read_counts, results, forlog):
    for target in targets:
        # long = l = BARSIZE, short = s = REDUCED_BAR
        # relaxed = r = ham < nearest ham, ultra relaxed = u = ham < max_ham
        # never short short (off-targets)
        # ultru ultra (uncertain origin) will be done in 'ultra relaxed matching' which is very slow

        # lbar, lham, rbar, rham
        parameters = (
            # 1. lrsu: 5' long relaxed, 3' short ultra
            (BARSIZE, target.ham_5p, REDUCED_BAR, args.hd),
            # 2. sulr: 5' short ultra, 3' long relaxed
            (REDUCED_BAR, args.hd, BARSIZE, target.ham_3p),
            # 3. srlu: 5' short relaxed, 5' long ultra
            (REDUCED_BAR, target.ham_5p, BARSIZE, args.hd),
            # 4. lusr: 5' long ultra, 5' short relaxed
            (BARSIZE, args.hd, REDUCED_BAR, target.ham_3p),
        )

        for idx, (lbar, lham, rbar, rham) in enumerate(parameters, start=1):
            for read, count in tuple(read_counts.items()):
                if (
                    hamming(read[:lbar], target.barcode_5p[:lbar]) < lham
                    and hamming(read[-rbar:], target.barcode_3p[-rbar:]) < rham
                ):
                    results[target.name][len(read) - target.len] += count
                    forlog[target.name][idx] += count
                    read_counts.pop(read)


def ultra_relaxed_matching(args, targets, read_counts, results, forlog):
    for read, count in tuple(read_counts.items()):
        # find out where a read MIGHT originate from
        matching_pcrs = []
        min_score = args.hd
        for target in targets:
            score = hamming(
                target.barcode_5p[:REDUCED_BAR], read[:REDUCED_BAR]
            ) + hamming(target.barcode_3p[-REDUCED_BAR:], read[-REDUCED_BAR:])
            if score < min_score:
                matching_pcrs = [target]
            elif score == min_score:
                matching_pcrs.append(target)

        if matching_pcrs:
            # ok this read might be something.
            # check lusu or sulu
            maybe_product = set()
            for target in matching_pcrs:
                # lusu
                if (
                    hamming(target.barcode_5p, read[:BARSIZE]) < args.hd
                    and hamming(target.barcode_3p[-REDUCED_BAR:], read[-REDUCED_BAR:])
                    < args.hd
                ):
                    maybe_product.add(target)
                # sulu
                if (
                    hamming(target.barcode_5p[:REDUCED_BAR], read[:REDUCED_BAR])
                    < args.hd
                    and hamming(target.barcode_3p, read[-BARSIZE:]) < args.hd
                ):
                    maybe_product.add(target)

            if len(maybe_product) == 1:
                target = maybe_product.pop()
                results[target.name][target.len - len(read)] += count
                read_counts.pop(read)  # remove analyzed reads
                forlog[target.name][5] += count
            elif len(maybe_product) > 1:
                min_score = args.hd * 2
                winner = []
                for target in maybe_product:
                    score = hamming(target.barcode_5p, read[:BARSIZE]) + hamming(
                        target.barcode_3p, read[-BARSIZE:]
                    )
                    if score < min_score:
                        min_score = score
                        winner = [target]
                    elif score == min_score:
                        winner.append(target)

                if len(winner) == 1:
                    target = winner[0]
                    results[target.name][target.len - len(read)] += count
                    read_counts.pop(read)  # remove analyzed reads
                    forlog[target.name][5] += count
                elif len(winner) > 1:
                    log = logging.getLogger("matching")
                    log.warning("ambigious reads:\n%s", read)
                    for target in winner:
                        log.warning(target.name)


def analyze_sample(parameters):
    args, sample, targets = parameters
    log = logging.getLogger("matching")

    results = {}
    for target in targets:
        results[target.name] = collections.defaultdict(int)

    filename = os.path.join(
        args.output_merged, "index%03i.extendedFrags.fastq.gz" % (sample["#"],)
    )
    log.debug("Sample %03i: Reading index from %s", sample["#"], filename)
    read_counts = count_fastq_sequences(filename)
    total_reads = sum(read_counts.values())

    logstats = {}
    for target in targets:
        # strict, lrsu, sulr, srlu, lusr, ur total
        logstats[target.name] = [0, 0, 0, 0, 0, 0]

    # STRICT matching (very fast, takes out anything we can easily assign to a fasta target)
    log.debug(
        "Sample %03i: Strict matching %i products with %i reads",
        sample["#"],
        len(read_counts),
        sum(read_counts.values()),
    )
    strict_matching(targets, read_counts, results, logstats)

    if not args.dr:
        # RELAXED matching (slower, but not very)
        log.debug(
            "Sample %03i: Relaxed matching %i products with %i reads",
            sample["#"],
            len(read_counts),
            sum(read_counts.values()),
        )
        relaxed_matching(args, targets, read_counts, results, logstats)

        # ULTRA RELAXED MATCHING (probably slow)
        # now anything left will be ultra relaxed matched. This opens up for a read to match multiple targets  so we need to check for this and
        # also any offtargets these primers may be producing will also show up here...maybe that should be a number to report...
        log.debug(
            "Sample %03i: Ultra relaxed matching %i products with %i reads",
            sample["#"],
            len(read_counts),
            sum(read_counts.values()),
        )
        ultra_relaxed_matching(args, targets, read_counts, results, logstats)

    return {
        "#": sample["#"],
        "name": sample["name"],
        "#reads": total_reads,
        "#unmatched": sum(read_counts.values()),
        "targets": results,
        "logstats": logstats,
        "unmatched": read_counts,
    }


def log_sample_stats(sample):
    log = logging.getLogger("matching")

    name_length = max(len(x) for x in sample["logstats"])
    log.info(
        "{:<{}} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}".format(
            "Target",
            name_length,
            "Total",
            "Strict",
            "lrsu",
            "sulr",
            "srlu",
            "lusr",
            "ur",
        )
    )
    for header, numbers in sample["logstats"].items():
        log.info(
            "{:<{}} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}".format(
                header, name_length, sum(numbers), *numbers
            )
        )

    unmatched_reads = []
    for read, count in sample["unmatched"].items():
        if (
            count >= SIGNIFICANT_UNMATCHED_THRESHOLD_READS
            and len(read) >= SIGNIFICANT_UNMATCHED_THRESHOLD_LENGTH
        ):
            unmatched_reads.append((count, read))

    if unmatched_reads:
        log.warning("Unmatched reads, index: %s", sample["#"])
        for count, read in sorted(unmatched_reads):
            log.warning(" % 6i: %s", count, read)

    log.info("Unmatched unique reads: %s", len(sample["unmatched"]))
    log.info("Unmatched total reads: %s", sum(sample["unmatched"].values()))


def analyze(args, samples, fasta):
    log = logging.getLogger("matching")

    results = []  # results are [{name: defaultdict(int), ...}, ...]
    pool = multiprocessing.Pool(processes=args.threads)
    pool_args = [(args, sample, fasta) for sample in samples]
    for output in pool.imap(analyze_sample, pool_args):
        log.info("Completed analysis of sample %i ..", output["#"])
        log_sample_stats(sample=output)

        output.pop("unmatched")
        results.append(output)

    if results:
        output_summary = {
            "#": "*",
            "name": "summary",
            "logstats": results[0].pop("logstats"),
            "unmatched": {},
        }

        for result in results[1:]:
            for name, numbers in result.pop("logstats").items():
                for idx, count in enumerate(numbers):
                    output_summary["logstats"][name][idx] += count

        log.info("Completed analysis of %i samples!", len(results))
        log_sample_stats(sample=output_summary)

    return results


def write_report(args, results, output_filename):
    workbook = xlsxwriter.Workbook(output_filename)
    thousandsepstyle = workbook.add_format({"num_format": "#,##0"})
    percentstyle = workbook.add_format({"num_format": "0.00%"})

    worksheet = workbook.add_worksheet("results")
    worksheet.set_column(2, 2, 20)
    worksheet.write(0, 0, "Index/Well")
    worksheet.write(0, 1, "Sample")
    worksheet.write(0, 2, "Target")
    worksheet.write(0, 3, "Reads")
    worksheet.write(0, 4, "WT")
    worksheet.write(0, 5, "Indel")
    worksheet.write(0, 6, "%Indel")
    # peak counting
    max_peaks = 0
    for sample in results:
        for target, peaks in sample["targets"].items():
            sigpeaks = len(
                [
                    p
                    for s, p in peaks.items()
                    if s
                    and p >= TARGET_LOCAL_THRESHOLD_READS
                    and (float(p) / sum(peaks.values()))
                    >= TARGET_LOCAL_THRESHOLD_PERCENTAGE
                ]
            )
            if sigpeaks > max_peaks:
                max_peaks = sigpeaks
    # add a header per peak
    peak_col_start = 7
    if max_peaks:
        worksheet.set_column(peak_col_start, max_peaks + peak_col_start - 1, 16)
    for col in range(peak_col_start, max_peaks + peak_col_start):
        worksheet.write(0, col, "peak{}".format(col - (peak_col_start - 1)))
        worksheet.write(
            0, col + max_peaks, "peak{}%".format(col - (peak_col_start - 1))
        )
    row = 1
    for sample in results:
        for target, peaks in sample["targets"].items():
            if sum(peaks.values()) >= TARGET_TOTAL_THRESHOLD:
                if args.cw:
                    worksheet.write(row, 0, convert_index_to_well(sample["#"]))
                else:
                    worksheet.write(row, 0, sample["#"])
                worksheet.write(row, 1, sample["name"])
                worksheet.write(row, 2, target)
                worksheet.write(row, 3, sum(peaks.values()), thousandsepstyle)
                worksheet.write(row, 4, peaks[0], thousandsepstyle)
                worksheet.write(
                    row, 5, sum(peaks.values()) - peaks[0], thousandsepstyle
                )
                worksheet.write(row, 6, "=F{0}/D{0}".format(row + 1), percentstyle)
                column = 7
                for indel_size, num in sorted(peaks.items()):
                    if indel_size:
                        if (
                            num > TARGET_LOCAL_THRESHOLD_READS
                            and (float(num) / sum(peaks.values()))
                            >= TARGET_LOCAL_THRESHOLD_PERCENTAGE
                        ):
                            if indel_size % 3 == 0:
                                indel_size = "{} (inframe)".format(indel_size)
                            worksheet.write(row, column, "{}".format(indel_size))
                            worksheet.write_formula(
                                row,
                                column + max_peaks,
                                "={0}/D{1}".format(num, row + 1),
                                percentstyle,
                            )
                            column += 1
                row += 1
    worksheet.autofilter(
        0, 0, row, (max_peaks * 2) + (peak_col_start - 1)
    )  # Same as above.

    workbook.close()


def write_json(args, targets, results, output_filename):
    data = {
        "samples": {},
        "targets": {target.name: target.seq for target in targets},
        "version": __version__,
        "settings": {
            "-cw": args.cw,
            "-dr": args.dr,
            "-hd": args.hd,
            "-pl": args.pl,
        },
    }

    for sample in results:
        result = {}
        for target_name, target_data in sample["targets"].items():
            peaks = dict(target_data)
            total_reads = sum(peaks.values())
            total_wildtype = peaks.pop(0, 0)
            total_indels = total_reads - total_wildtype

            if total_reads:
                result[target_name] = {
                    "reads": total_reads,
                    "reads_wildtype": total_wildtype,
                    "reads_mutant": total_indels,
                    "peaks": peaks,
                }

        data["samples"][sample["#"]] = {
            "name": sample["name"],
            "reads": sample["#reads"],
            "unidentified": sample["#unmatched"],
            "targets": result,
        }

    with open(output_filename, "wt") as handle:
        handle.write(json.dumps(data))


def parse_args(argv):
    parser = ArgumentParser(
        description="hamplicons - estimating indel sizes in amplicons using Hamming distances",
        epilog="Needs flash available in PATH " "(http://ccb.jhu.edu/software/FLASH/)",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "output_prefix", nargs="?", default="output", help="Prefix for output files"
    )

    parser.add_argument(
        "-sm", action="store_true", help="skip merging of fastq paired end reads"
    )
    parser.add_argument("-dr", action="store_true", help="Deactivate relaxed matching")
    parser.add_argument("-cw", action="store_true", help="Convert Index to Well")
    parser.add_argument(
        "-tf",
        type=str,
        default="targets.fa",
        help="Fasta file with WT PCR",
        metavar="targets",
    )
    parser.add_argument(
        "-dd",
        type=str,
        default="fastq",
        help="data directory. Directory containing fastq.gz files",
        metavar="directory",
    )
    parser.add_argument(
        "-dd-recursive",
        action="store_true",
        help="Search data directoryrecursively",
    )
    parser.add_argument(
        "-pl",
        type=list,
        nargs=2,
        help="Define layout for correct index translation, e.g. -pl 8 12",
        default=[8, 12],
        metavar="int",
    )
    parser.add_argument(
        "-hd",
        type=int,
        help="Allow this many mismatches in relaxed matching",
        default=4,
        metavar="int",
    )

    parser.add_argument(
        "--threads",
        type=int,
        help="Number of threads for read merging/analysis",
        default=multiprocessing.cpu_count(),
        metavar="N",
    )

    parser.add_argument(
        "--version", action="version", version="%(prog)s v" + __version__
    )

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    # if we are merging, check that flash exists
    if not (args.sm or shutil.which("flash")):
        print("\nFatal error:\n")
        print(
            "This program requires flash to be available on path.\nI.e. 'which flash' needs to return something.\nFlash is available here: http://ccb.jhu.edu/software/FLASH/"
        )
        return 1

    args.output_merged = args.output_prefix + ".merged"
    os.makedirs(args.output_merged, exist_ok=True)

    setup_logging(args.output_prefix + ".log")
    log = logging.getLogger("main")

    # read targets
    log.info("Reading FASTA file %r", args.tf)
    targets = build_targets(read_fasta(args.tf), max_ham=args.hd)
    if targets is None:
        return 1

    # Collect sample files (FASTQ) and names
    samples = collect_fastq_files(log, args.dd, recursive=args.dd_recursive)
    if samples is None:
        return 1
    # merge reads
    if not args.sm:
        if not merge_fastq(args.output_merged, samples, threads=args.threads):
            return 1
    else:
        log.info("Skipping FASTQ merging")

    # analyze reads
    log.info("Analyzing FASTQ files")
    results = analyze(args, samples, targets)

    log.info("Writing report to %r", args.output_prefix + ".xlsx")
    write_report(
        args=args,
        results=results,
        output_filename=args.output_prefix + ".xlsx",
    )

    log.info("Writing JSON to %r", args.output_prefix + ".json")
    write_json(
        args=args,
        targets=targets,
        results=results,
        output_filename=args.output_prefix + ".json",
    )

    log.info("Done")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
