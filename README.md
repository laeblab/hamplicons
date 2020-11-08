# hamplicons - estimating indel sizes in amplicons using Hamming distances

## Requirements

`hamplicons` requires Python3 and [FLASH](https://ccb.jhu.edu/software/FLASH/) (Fast Length Adjustment of SHort reads). The `flash` executable must be located in the current user's `PATH`.

## Installation

``` bash
$ git clone https://github.com/laeblab/hamplicons.git
$ cd hamplicons
$ python3 setup.py install
$ hamplicons -h
```

## Usage

By default `hamplicons` expects a FASTA file named `targets.fa` (see below) and a directory containing FASTQ files named `fastq` to be located in the current working directory. Additionally, the FASTQ files are expected to follow the standard Illumina naming scheme, namely `.*_S[0-9]+_L[0-9]+_R[12]_[0-9]+.fastq.gz`. Output is written using an `output` prefix by default

``` bash
$ ls -p
fastq/  targets.fa
$ hamplicons
$ ls -p
fastq/  output.xlsx  output.log  output.merged/  targets.fa
```

The location of both the FASTA file and the FASTQ directory, as well as the output prefix can be changed using command-line options:

``` bash
$ hamplicons --data-directory /path/to/fastqs/ --targets-fasta /path/to/targets.fa my_ouput_prefix
$ ls -p
my_output_prefix.xlsx  my_output_prefix.log  my_output_prefix.merged/
```

## Target Amplicons

The `targets.fa` file is expected to contain one or more wild-type amplicons in FASTA format. The names of amplicon are used in the output report and must be unique. In addition, as `hamplicons` assigns sequences to wild-type amplicons by comparing them to the first and last bp of the wild-type amplicon, `hamplicon` cannot differentiate wild-type amplicons that only differ outside of the first/last 30 bp.
