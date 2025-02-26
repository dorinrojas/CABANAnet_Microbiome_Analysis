# Class 5: Ecological characterization of the microbiome

- - - -

## Taxon abundance and depth estimation (Bowtie2, InStrain)

Once the taxonomica annotation and (optional) dereplication is performed it is require to conduct the estimation of abundance and depth of each MAG. This bioinformatic process calculates the number of reads mapped to the bins divided by the total number of reads in the sample. In easy words, this represent the amount of times the bin was sequenced.

Combined with the taxonomical classification of each bins, this value (known as relative abundance) or percentage represents the abundace of identified species in the gut. This can be interpreted as how predominant or rare different bacteria are in the gastrointestinal tract.

For research purposes, the is one of the main aspects to evaluate in human microbiome studies. The diversity and the abundance of bacterial species can drastically change according to the diet, environmental factors, and disease. For instance, an abundance above 0.2% of *Akkermansia* sp. seems to decrease the risk of metabolic syndrome [(Zhou et al. 2021)](https://www.tandfonline.com/doi/full/10.2147/DMSO.S311388). Additionally, the abundance or presence of certain taxa can influence the abundance of other species. In this regard, the presence or Ruminococcaceae and Lachnospiraceae were associated to *Akkermansia* abundance [(Zhou et al. 2021)](https://www.tandfonline.com/doi/full/10.2147/DMSO.S311388).

In general, the analysis of abundance is relevant to access the relation between environmental factors and human health and the microbiome composition. Hence, this comes of one of the final aims of metagenomic pipelines.

As mentioned above, the main method to conduct this estimation is by mapping the reads to the bins resulting from the dereplication process or the bins bins after quality filtering (the ones used for taxonomic classification). This generates a `.bam` or `.cram` file which contains the information of those mapped reads. These files are used to extract quantitative information as a table were the abundance (amount of reads mapped) per samples is presented.

There are several methods to perform the aligment and the abundance estimation. During this workshop, the tools Bowtie2 and inStrain will be employed.

[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is an ultrafast and memory-efficient tools for alignments. It is particularly good at relatively long genomes and metagenomes. This memory-efficiency is due to the creation of its own index through the option `bowtie2-build` that favors the aligments. The Bowtie2 official page holds record of representative indexes for long genomes (e.g. human genomes). However, the tools provides the options of building an index based on specific data. In this case, an index of the analyzed pangenome is going to be created.

[inStrain](https://github.com/MrOlm/inStrain) is a python program that performs several analysis within metagenomic datasets. It can analyze co-occurence of genome populations, coverage, microdiversity, single nucleotide polymorphisms (SNPs), among others. This workshop is particularly interested in the analysis of microdiversity, which is perform with the option `profile` of the software. For this, inStrain take the representative genomes and use the previosuly generated aligment files (`.bam`) to calculate metrics for each genomes.

The output from this analysis are abundance per bin that are later combined to the taxonomy classification of GTDB-tk to allows the interpretation and visualization in different graph. Additional information, such as the functional annotations can also be associated to this interpretation. This will be discussed as the last section of this workshop.

> Notice that taxonomic classification is not required for the abundance estimation as it is later appended to the abundace table. Therefore, the abundance estimation can be performed before or after the this analysis, depending on the user-preferences.

### Running bowtie2

The first step towards running Bowtie2 for mapping is creating the pangenome index. This is an easy command to be used within the same container.

```bash
(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/bowtie2-2.5.4.sif bowtie2-build -h
Bowtie 2 version 2.5.4 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>
    reference_in            comma-separated list of files with ref sequences
    bt2_index_base          write bt2 data to files with this dir/basename
*** Bowtie 2 indexes will work with Bowtie v1.2.3 and later. ***
Options:
    -f                      reference files are Fasta (default)
    -c                      reference sequences given on cmd line (as
                            <reference_in>)
    --large-index           force generated index to be 'large', even if ref
                            has fewer than 4 billion nucleotides
    --debug                 use the debug binary; slower, assertions enabled
    --sanitized             use sanitized binary; slower, uses ASan and/or UBSan
    --verbose               log the issued command
    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting
    -p/--packed             use packed strings internally; slower, less memory
    --bmax <int>            max bucket sz for blockwise suffix-array builder
    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)
    --dcv <int>             diff-cover period for blockwise (default: 1024)
    --nodc                  disable diff-cover (algorithm becomes quadratic)
    -r/--noref              don't build .3/.4 index files
    -3/--justref            just build .3/.4 index files
    -o/--offrate <int>      SA is sampled every 2^<int> BWT chars (default: 5)
    -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)
    --threads <int>         # of threads
    --seed <int>            seed for random number generator
    -q/--quiet              verbose output (for debugging)
    --h/--help              print this message and quit
    --version               print version information and quit
```

This option created several files within a specific folder with the indicated prefix. It is recommended to create a singular folder with in the `14-bowtie` path for this files alone. The prefix is the `<bt2_index_name>` indicated in the usage example. Code for the index to be created and specify the use of all available threads.

> If the job takes a while running, the results in the common repository `/home/public/met-workshop`

The output from the index creation pipeline should look something similar to this, depending on the prefix designated in the code.

```bash
(base) [dorian.rojas@accessnode test]$ ls 14-bowtie/1-index/
pangenome.1.bt2  pangenome.2.bt2  pangenome.3.bt2  pangenome.4.bt2  pangenome.rev.1.bt2  pangenome.rev.2.bt2
```

After this, the code for mapping can be designed. The command for mapping is just `bowtie2`.

```bash
/usr/local/env-execute: line 3: exec: bowtie: not found
(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/bowtie2-2.5.4.sif bowtie2 -h
Bowtie 2 version 2.5.4 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
Usage:
  bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]

  <bt2-idx>  Index filename prefix (minus trailing .X.bt2).
             NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.
  <m1>       Files with #1 mates, paired with files in <m2>.
             Could be gziped (extension: .gz) or bzip2ed (extension: .bz2).
  <m2>       Files with #2 mates, paired with files in <m1>.
             Could be gziped (extension: .gz) or bzip2ed (extension: .bz2).
  <r>        Files with unpaired reads.
             Could be gziped (extension: .gz) or bzip2ed (extension: .bz2).
  <i>        Files with interleaved paired-end FASTQ/FASTA reads
             Could be gziped (extension: .gz) or bzip2ed (extension: .bz2).
  <bam>      Files are unaligned BAM sorted by read name.
  <sam>      File for SAM output (default: stdout)

  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.

Options (defaults in parentheses):

 Input:
  -q                 query input files are FASTQ .fq/.fastq (default)
  --tab5             query input files are TAB5 .tab5
  --tab6             query input files are TAB6 .tab6
  --qseq             query input files are in Illuminas qseq format
  -f                 query input files are (multi-)FASTA .fa/.mfa
  -r                 query input files are raw one-sequence-per-line
  -F k:<int>,i:<int> query input files are continuous FASTA where reads
                     are substrings (k-mers) extracted from the FASTA file
                     and aligned at offsets 1, 1+i, 1+2i ... end of reference
  -c                 <m1>, <m2>, <r> are sequences themselves, not files
  -s/--skip <int>    skip the first <int> reads/pairs in the input (none)
  -u/--upto <int>    stop after first <int> reads/pairs (no limit)
  -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
  -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)
  --trim-to [3:|5:]<int> trim reads exceeding <int> bases from either 3 or 5 end
                     If the read end is not specified then it defaults to 3 (0)
  --phred33          qualities are Phred+33 (default)
  --phred64          qualities are Phred+64
  --int-quals        qualities encoded as space-delimited integers

 Presets:                 Same as:
  For --end-to-end:
   --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
   --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50
   --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)
   --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

  For --local:
   --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
   --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
   --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
   --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

 Alignment:
  -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
  -L <int>           length of seed substrings; must be >3, <32 (22)
  -i <func>          interval between seed substrings w/r/t read len (S,1,1.15)
  --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)
  --dpad <int>       include <int> extra ref chars on sides of DP table (15)
  --gbar <int>       disallow gaps within <int> nucs of read extremes (4)
  --ignore-quals     treat all quality values as 30 on Phred scale (off)
  --nofw             do not align forward (original) version of read (off)
  --norc             do not align reverse-complement version of read (off)
  --no-1mm-upfront   do not allow 1 mismatch alignments before attempting to
                     scan for the optimal seeded alignments
  --end-to-end       entire read must align; no clipping (on)
   OR
  --local            local alignment; ends might be soft clipped (off)

 Scoring:
  --ma <int>         match bonus (0 for --end-to-end, 2 for --local)
  --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)
  --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
  --rdg <int>,<int>  read gap open, extend penalties (5,3)
  --rfg <int>,<int>  reference gap open, extend penalties (5,3)
  --score-min <func> min acceptable alignment score w/r/t read length
                     (G,20,8 for local, L,-0.6,-0.6 for end-to-end)

 Reporting:
  (default)          look for multiple alignments, report best, with MAPQ
   OR
  -k <int>           report up to <int> alns per read; MAPQ not meaningful
   OR
  -a/--all           report all alignments; very slow, MAPQ not meaningful

 Effort:
  -D <int>           give up extending after <int> failed extends in a row (15)
  -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)

 Paired-end:
  -I/--minins <int>  minimum fragment length (0)
  -X/--maxins <int>  maximum fragment length (500)
  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
  --no-mixed         suppress unpaired alignments for paired reads
  --no-discordant    suppress discordant alignments for paired reads
  --dovetail         concordant when mates extend past each other
  --no-contain       not concordant when one mate alignment contains other
  --no-overlap       not concordant when mates overlap at all

 BAM:
  --align-paired-reads
                     Bowtie2 will, by default, attempt to align unpaired BAM reads.
                     Use this option to align paired-end reads instead.
  --preserve-tags    Preserve tags from the original BAM record by
                     appending them to the end of the corresponding SAM output.

 Output:
  -t/--time          print wall-clock time taken by search phases
  --un <path>        write unpaired reads that didnt align to <path>
  --al <path>        write unpaired reads that aligned at least once to <path>
  --un-conc <path>   write pairs that didnt align concordantly to <path>
  --al-conc <path>   write pairs that aligned concordantly at least once to <path>
    (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g.
    --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)
  --quiet            print nothing to stderr except serious errors
  --met-file <path>  send metrics to file at <path> (off)
  --met-stderr       send metrics to stderr (off)
  --met <int>        report internal counters & metrics every <int> secs (1)
  --no-unal          suppress SAM records for unaligned reads
  --no-head          suppress header lines, i.e. lines starting with @
  --no-sq            suppress @SQ header lines
  --rg-id <text>     set read group id, reflected in @RG line and RG:Z: opt field
  --rg <text>        add <text> ("lab:value") to @RG line of SAM header.
                     Note: @RG line only printed when --rg-id is set.
  --omit-sec-seq     put '*' in SEQ and QUAL fields for secondary alignments.
  --sam-no-qname-trunc
                     Suppress standard behavior of truncating readname at first whitespace
                     at the expense of generating non-standard SAM.
  --xeq              Use '='/'X', instead of 'M,' to specify matches/mismatches in SAM record.
  --soft-clipped-unmapped-tlen
                     Exclude soft-clipped bases when reporting TLEN.
  --sam-append-comment
                     Append FASTA/FASTQ comment to SAM record.
  --sam-opt-config <config>
                     Use <config>, example '-MD,YP,-AS', to toggle SAM Optional fields.

 Performance:
  -p/--threads <int> number of alignment threads to launch (1)
  --reorder          force SAM output order to match order of input reads
  --mm               use memory-mapped I/O for index; many 'bowtie's can share

 Other:
  --qc-filter        filter out reads that are bad according to QSEQ filter
  --seed <int>       seed for random number generator (0)
  --non-deterministic
                     seed rand. gen. arbitrarily instead of using read attributes
  --version          print version information and quit
  -h/--help          print this usage message
```

The code is simple and requires the `-x` function with the path towards the index with the prefix, `-1` and `-2` with the path to the raw reads, and the output files that are defined as a `.sam` file. Additionally, it is only required to add the command for threads to use.

Notice the output from Bowtie2 is in the `.sam` format instead of the required `.bam`. Hence three additional steps need to be performed after running the mapping. The transformation from `.sam` to `.sam`, the sorting, and the indexing of the `.bam` file. In this regard, several options from the toolkit `samtools` are required.

```bash
(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/samtools-1.21.sif samtools --help

Program: samtools (Tools for alignments in the SAM format)
Version: 1.21 (using htslib 1.21)

Usage:   samtools <command> [options]

Commands:
  -- Indexing
     dict           create a sequence dictionary file
     faidx          index/extract FASTA
     fqidx          index/extract FASTQ
     index          index alignment

  -- Editing
     calmd          recalculate MD/NM tags and '=' bases
     fixmate        fix mate information
     reheader       replace BAM header
     targetcut      cut fosmid regions (for fosmid pool only)
     addreplacerg   adds or replaces RG tags
     markdup        mark duplicates
     ampliconclip   clip oligos from the end of reads

  -- File operations
     collate        shuffle and group alignments by name
     cat            concatenate BAMs
     consensus      produce a consensus Pileup/FASTA/FASTQ
     merge          merge sorted alignments
     mpileup        multi-way pileup
     sort           sort alignment file
     split          splits a file by read group
     quickcheck     quickly check if SAM/BAM/CRAM file appears intact
     fastq          converts a BAM to a FASTQ
     fasta          converts a BAM to a FASTA
     import         Converts FASTA or FASTQ files to SAM/BAM/CRAM
     reference      Generates a reference from aligned data
     reset          Reverts aligner changes in reads

  -- Statistics
     bedcov         read depth per BED region
     coverage       alignment depth and percent coverage
     depth          compute the depth
     flagstat       simple stats
     idxstats       BAM index stats
     cram-size      list CRAM Content-ID and Data-Series sizes
     phase          phase heterozygotes
     stats          generate stats (former bamcheck)
     ampliconstats  generate amplicon specific stats

  -- Viewing
     flags          explain BAM flags
     head           header viewer
     tview          text alignment viewer
     view           SAM<->BAM<->CRAM conversion
     depad          convert padded BAM to unpadded BAM
     samples        list the samples in a set of SAM/BAM/CRAM files

  -- Misc
     help [cmd]     display this help message or help for [cmd]
     version        detailed version information
```

For transforming the output the `view` command is required. While for sorting and indexing, `sort` and `index` are used.

Particularly, the transformation from `.sam` to `.bam` allows the extraction and indication of the information of interest. These are defined using the `-F, -f` commands to exclude or include (respectively) depending on the code. This flags work under numbers codes that state for different alignment information. The complete list of combination is defined on the official webpage of [Samtools](https://www.htslib.org/doc/samtools-flags.html). This can be complex to understand; however, 'calculators' have been developed to facilitate this process (check this ['Decoding SAM flags page'](https://broadinstitute.github.io/picard/explain-flags.html)).

The commands for all these options are simple. However, the toolkit holds a great variety of arguments. To avoid misunderstanding, here are the basic commands to use.

```bash
#Tranforming from SAM to BAM
$CTN_PATH/samtools-1.21.sif samtools view \
    -F 4 -bS $output/${sample}.sam \
        > $output/${sample}-raw.bam

# `view` transforms. 
# `-F 4` indicates extraction of aligned sequences
# `-S` indicated to autodetec the input file extension and `-b` to generate a .bam output

#Sort
$CTN_PATH/samtools-1.21.sif samtools sort -@ 64 \
    $output/${sample}-raw.bam \
    -o $output/${sample}.sorted.bam

# `-@` indictes threads
# `-o` for output

#Indexing
$CTN_PATH/samtools-1.21.sif samtools index $output/${sample}.sorted.bam

#This index the sorted file, the output if default as the file name with `.bai` at the end

rm $output/${sample}-raw.bam 
#Removes intermediate files
```

Write the code for mapping, transforming, sorting, and indexing the reads to the pangenome. Follow the template of code provided above. All these functions can be included into the same file to favor reproducibiltiy. In this case, the raw reads should be mapped individually. Hence, the `accessions.txt`, `batch.sh`, and `for` loop are required.

> If the job takes a while running, the results in the common repository `/home/public/met-workshop`

The output from this mapping should look something similar to the following. There should a `.sam`, `.bam`, and `.bam.bai` per sample. The `.bai` is the index of the mapping and the `.bam` is the aligment that will be used as input in inStrain.

```bash
(base) [dorian.rojas@accessnode test]$ ll 14-bowtie/
total 790844
drwxrwxr-x 2 dorian.rojas dorian.rojas       180 feb 24 13:08 1-index
-rw-rw-r-- 1 dorian.rojas dorian.rojas 360140221 feb 25 12:08 SRR8555091.sam
-rw-rw-r-- 1 dorian.rojas dorian.rojas  18426729 feb 25 12:08 SRR8555091.sorted.bam
-rw-rw-r-- 1 dorian.rojas dorian.rojas    127128 feb 25 12:08 SRR8555091.sorted.bam.bai
-rw-rw-r-- 1 dorian.rojas dorian.rojas 390287793 feb 25 12:07 SRR9988205.sam
-rw-rw-r-- 1 dorian.rojas dorian.rojas  40698950 feb 25 12:08 SRR9988205.sorted.bam
-rw-rw-r-- 1 dorian.rojas dorian.rojas    117496 feb 25 12:08 SRR9988205.sorted.bam.bai
```

### Running instrain

Prior to running inStrain, the creation of a `scaffolds to bin` file is necessary. This table indicate the contigs and the bins in which they are located. This is a required command in inStrain as the mapped file is a pangenome. If this option and file are not provided, the tools considers that all contigs are part of the sample genome, resulting in innacurate estimation.

To generate this table, a command prepared by the coding team of dRep is recommended in the inStrain documentation. This correspond to a python code (presented below) that runs by the command `python parse_stb.py --reverse -f <path/to/genomes> -o <path/to/output>`.

```python
#!/usr/bin/env python

import os
import sys
import gzip
import textwrap
import argparse

def extract_bins(fasta, stb_file, out_base):
    if (out_base != '') & (out_base[-1] != '/'):
        out_base += '_'

    # Make scaffold to bin dictionary
    stb = {}
    with open(stb_file) as stb_reader:
        for line in stb_reader:
            line = line.strip()

            if line.startswith('#') or line.startswith('scaffold_name'):
                continue

            scaffold, bin = line.split('\t')[:2]
            stb[scaffold.strip()] = bin.strip()

    # Parse the FASTA file and write to bins
    if fasta.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'

    with open_func(fasta, mode) as handle:
        record_id = ''
        record_seq = ''
        for line in handle:
            if line.startswith('>'):
                if record_id and record_id in stb:
                    bin = stb[record_id]
                    with open(f"{out_base}{bin}.fa", 'a') as nw:
                        nw.write(f">{record_id}\n{record_seq}\n")
                record_id = line[1:].strip().split(' ')[0]
                record_seq = ''
            else:
                record_seq += line.strip()

        # Write the last record
        if record_id and record_id in stb:
            bin = stb[record_id]
            with open(f"{out_base}{bin}.fa", 'a') as nw:
                nw.write(f">{record_id}\n{record_seq}\n")

import os
import gzip

def gen_stb(fastas):
    if ((len(fastas) == 1) & (not fastas[0].endswith('.gz'))):
        # See if this is a text file, not a fasta file
        text_list = True
        genomes = []
        with open(fastas[0], 'r') as o:
            for line in o.readlines():
                if line.startswith('>'):
                    text_list = False
                    break
                else:
                    genomes.append(line.strip())
        if text_list:
            print("Treating .fasta input as list")
            fastas = genomes

    stb = {}
    for fasta in fastas:
        bin = os.path.basename(fasta)
        if fasta.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        with open_func(fasta, mode) as handle:
            for line in handle:
                if line.startswith('>'):
                    id = line[1:].strip().split(' ')[0]
                    stb[id] = bin

    return stb


def print_stb(stb, output):
    file = open(output,'w')
    for key in sorted(stb, key=stb.get):
        file.write("{0}\t{1}\n".format(key,stb[key]))
    file.close()


class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=textwrap.dedent('''\
         \n
         The program has two uses related to scaffold to bin (.stb) files.
         .stb files should be tab-separated, with no header, and two columns: scaffold and bin

         Use 1) Pass a list of genomes to generate a .stb file.

         Example:
         parse_stb.py --reverse -f dereplicate_genomes/* -o representitve_genomes.stb

         Use 2) Pass a single .fasta file and a scaffold to bin file (.stb) to generate a number of
         fasta files based on the .stb file.

         Example:
         parse_stb.py -f concat_genomes.fasta -s scaffold_to_bin.tsv -o genomeList_1
         '''))

    parser.add_argument('-s','--stb',help='scaffold to bin file')
    parser.add_argument('-f','--fasta',help='fasta file to extract scaffolds from. Will treat as compressed if ends in .gz. This can also be a single text file with a genome on each line',nargs='*')
    parser.add_argument('-o','--output',help='output base name', default = '')

    parser.add_argument('--reverse',help='generate a stb from a list of genomes',\
                        action = "store_true")
    args = parser.parse_args()

    if args.reverse == False:
        if len(args.fasta) != 1:
            print("must give one and only one fasta file")
            sys.exit()
        extract_bins(args.fasta[0], args.stb, args.output)

    if args.reverse == True:
        stb = gen_stb(args.fasta)
        print_stb(stb, args.output)
```

Copy and paste the code into a `parse_stb.py` file created in the console and compile it with the code provided above. Take consideration of were the final bins are located and save the table as `representative_genomes.stb` in a `15-instrain` folder.

The output file should look similar to this.

```vim
SRR8555091_bin.1.permissive_0   SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta
SRR8555091_bin.1.permissive_1   SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta
SRR8555091_bin.1.permissive_2   SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta
SRR8555091_bin.1.permissive_3   SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta
SRR8555091_bin.1.permissive_4   SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta
SRR8555091_bin.1.permissive_5   SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta
SRR8555091_bin.1.permissive_6   SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta
SRR8555091_bin.1.permissive_7   SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta
SRR8555091_bin.1.permissive_8   SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta
SRR8555091_bin.1.permissive_9   SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta
```

After this, the code for inStrain can be evaluated. The call command for inStrain is `inStrain` and, as mentioned above, the module of interest is `profile`.

```bash
(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/instrain-1.9.0.sif inStrain profile -h
usage: inStrain profile [-o OUTPUT] [--use_full_fasta_header] [--force_compress] [-p PROCESSES] [-d] [-h] [--version]
                        [-l MIN_READ_ANI] [--min_mapq MIN_MAPQ] [--max_insert_relative MAX_INSERT_RELATIVE]
                        [--min_insert MIN_INSERT] [--pairing_filter {paired_only,non_discordant,all_reads}]
                        [--priority_reads PRIORITY_READS] [--maximum_reads MAXIMUM_READS] [--detailed_mapping_info]
                        [-c MIN_COV] [-f MIN_FREQ] [-fdr FDR] [-g GENE_FILE] [-s [STB [STB ...]]] [--mm_level]
                        [--skip_mm_profiling] [--database_mode] [--min_scaffold_reads MIN_SCAFFOLD_READS]
                        [--min_genome_coverage MIN_GENOME_COVERAGE] [--min_snp MIN_SNP] [--store_everything]
                        [--scaffolds_to_profile SCAFFOLDS_TO_PROFILE] [--rarefied_coverage RAREFIED_COVERAGE]
                        [--window_length WINDOW_LENGTH] [--skip_genome_wide] [--skip_plot_generation]
                        bam fasta

REQUIRED:
  bam                   Sorted .bam file
  fasta                 Fasta file the bam is mapped to

I/O PARAMETERS:
  -o OUTPUT, --output OUTPUT
                        Output prefix (default: inStrain)
  --use_full_fasta_header
                        Instead of using the fasta ID (space in header before space), use the full header. Needed for some
                        mapping tools (including bbMap) (default: False)
  --force_compress      Force compression of all output files (default: False)

SYSTEM PARAMETERS:
  -p PROCESSES, --processes PROCESSES
                        Number of processes to use (default: 6)
  -d, --debug           Make extra debugging output (default: False)
  -h, --help            show this help message and exit
  --version             show programs version number and exit

READ FILTERING OPTIONS:
  -l MIN_READ_ANI, --min_read_ani MIN_READ_ANI
                        Minimum percent identity of read pairs to consensus to use the reads. Must be >, not >= (default: 0.95)
  --min_mapq MIN_MAPQ   Minimum mapq score of EITHER read in a pair to use that pair. Must be >, not >= (default: -1)
  --max_insert_relative MAX_INSERT_RELATIVE
                        Multiplier to determine maximum insert size between two reads - default is to use 3x median insert
                        size. Must be >, not >= (default: 3)
  --min_insert MIN_INSERT
                        Minimum insert size between two reads - default is 50 bp. If two reads are 50bp each and overlap
                        completely, their insert will be 50. Must be >, not >= (default: 50)
  --pairing_filter {paired_only,non_discordant,all_reads}
                        How should paired reads be handled?
                        paired_only = Only paired reads are retained
                        non_discordant = Keep all paired reads and singleton reads that map to a single scaffold
                        all_reads = Keep all reads regardless of pairing status (NOT RECOMMENDED; See documentation for deatils)
                         (default: paired_only)
  --priority_reads PRIORITY_READS
                        The location of a list of reads that should be retained regardless of pairing status (for example long
                        reads or merged reads). This can be a .fastq file or text file with list of read names (will assume
                        file is compressed if ends in .gz (default: None)
  --maximum_reads MAXIMUM_READS
                        Maximum number of reads. Requires sambamba to do the subsetting. (default: None)

READ OUTPUT OPTIONS:
  --detailed_mapping_info
                        Make a detailed read report indicating deatils about each individual mapped read (default: False)

VARIANT CALLING OPTIONS:
  -c MIN_COV, --min_cov MIN_COV
                        Minimum coverage to call an variant (default: 5)
  -f MIN_FREQ, --min_freq MIN_FREQ
                        Minimum SNP frequency to confirm a SNV (both this AND the FDR snp count cutoff must be true to call a
                        SNP). (default: 0.05)
  -fdr FDR, --fdr FDR   SNP false discovery rate- based on simulation data with a 0.1 percent error rate (Q30) (default: 1e-06)

GENE PROFILING OPTIONS:
  -g GENE_FILE, --gene_file GENE_FILE
                        Path to prodigal .fna genes file. If file ends in .gb or .gbk, will treat as a genbank file
                        (EXPERIMENTAL; the name of the gene must be in the gene qualifier) (default: None)

GENOME WIDE OPTIONS:
  -s [STB [STB ...]], --stb [STB [STB ...]]
                        Scaffold to bin. This can be a file with each line listing a scaffold and a bin name, tab-seperated.
                        This can also be a space-seperated list of .fasta files, with one genome per .fasta file. If nothing is
                        provided, all scaffolds will be treated as belonging to the same genome (default: [])

READ ANI OPTIONS:
  --mm_level            Create output files on the mm level (see documentation for info) (default: False)
  --skip_mm_profiling   Dont perform analysis on an mm level; saves RAM and time; impacts plots and raw_data (default: False)

PROFILE OPTIONS:
  --database_mode       Set a number of parameters to values appropriate for mapping to a large fasta file. Will set:
                        --min_read_ani 0.92 --skip_mm_profiling --min_genome_coverage 1 (default: False)
  --min_scaffold_reads MIN_SCAFFOLD_READS
                        Minimum number of reads mapping to a scaffold to proceed with profiling it (default: 1)
  --min_genome_coverage MIN_GENOME_COVERAGE
                        Minimum number of reads mapping to a genome to proceed with profiling it. MUST profile .stb if this is
                        set (default: 0)
  --min_snp MIN_SNP     Absolute minimum number of reads connecting two SNPs to calculate LD between them. (default: 20)
  --store_everything    Store intermediate dictionaries in the pickle file; will result in significantly more RAM and disk
                        usage (default: False)
  --scaffolds_to_profile SCAFFOLDS_TO_PROFILE
                        Path to a file containing a list of scaffolds to profile- if provided will ONLY profile those scaffolds
                        (default: None)
  --rarefied_coverage RAREFIED_COVERAGE
                        When calculating nucleotide diversity, also calculate a rarefied version with this much coverage
                        (default: 50)
  --window_length WINDOW_LENGTH
                        Break scaffolds into windows of this length when profiling (default: 10000)

OTHER  OPTIONS:
  --skip_genome_wide    Do not generate tables that consider groups of scaffolds belonging to genomes (default: False)
  --skip_plot_generation
                        Do not make plots (default: False)
```

Write a code that performs the abundance estimation taking in consideration the `.stb` file, 64 threads, and stores the output in a `15-instrain/$sample` folder. Take in consideration that the `for` loop will be needed. In addition, state the `--database_mode` and `--skip_plot_generation` flags, as the index is provided and diversity graph will be later analyzed.

> If the job takes a while running, the results in the common repository `/home/public/met-workshop`

inStrain created different folders within the output directory. However, the `output` folder is the one containing the tables with the abundance information.

```bash
(base) [dorian.rojas@accessnode test]$ ls 15-instrain/SRR9988205/output/
SRR9988205_genome_info.tsv  SRR9988205_mapping_info.tsv   SRR9988205_SNVs.tsv
SRR9988205_linkage.tsv      SRR9988205_scaffold_info.tsv
```

The `scaffold_info.tsv` have general information about the provided scaffolds (e.g. coverage, breadth, nucleotide diversity, etc...). The `mapping_info.tsv` have the number of reads that mapped to each scaffold, and additional metrics. The `SNVs.tsv` describes mutations found in the different scaffolds. The `linkage.tsv` have the SNPs linkage (how likely two mutations are presented together). Finally, the `genome_info.tsv` have the metrics of interest for the abundance estimation.

```vim
genome  coverage        breadth nucl_diversity  length  true_scaffolds  detected_scaffolds      coverage_median coverage_std    coverage_SEM     breadth_minCov  breadth_expected        nucl_diversity_rarefied conANI_reference        popANI_reference        iRep     iRep_GC_corrected       linked_SNV_count        SNV_distance_mean       r2_mean d_prime_mean    consensus_divergent_sitespopulation_divergent_sites      SNS_count       SNV_count       filtered_read_pair_count        reads_unfiltered_pairs  reads_mean_PID   divergent_site_count    reads_unfiltered_reads
SRR9988205_bin.1.permissive_filtered_kept_contigs.fasta 5.496498350383906       0.994117863265243       0.0008198111233099      3781789  781     781     5       2.6698752312885534      0.001402173454713       0.6071142520114158      0.9921982564903672      0.0      0.9999629787393434      0.999998257823028               False   203     26.88177339901478       0.627394623981047       0.9875701459034792       85      4       3       459     92969   94802   0.9982983789351916      462     191829
SRR9988205_bin.2.permissive_filtered_kept_contigs.fasta 9.223122918012026       0.9982518369063236      0.0005069402542589      4456678  545     545     9       4.186153449138965       0.0020076430886163      0.8682996617660059      0.9997095321144326      0.0      0.9999790683179094      0.9999994831683434      1.9121057585254528      False   525     30.02285714285714       0.5179457671047509       0.9642752701766494      81      2       2       676     190590  194320  0.998575959217324       678     390820
```

In this sense, some important concepts need to be understand to fully interpret this table.

Concept|Definition
----|----
Coverage|Average number of reads mapping to a region, a mesure of depth
Breadth|Percentage of bases in a region covered by at least one read
Expected breadth|Breadth ecpected if reads were evenly distribuited

In general, the breadth would allows us to determinate whether or not a scaffold is considered part of the genome. For instance, if the breadth is lower than the expected breadth, the reads are aligned to a single part of the region (e.g. transposon).

This table with metrics will be used for further analysis.

## Ecological indexes estimation and metadata analysis (R command line)

The output abundance tables have to be combined with the taxonomic classification prior the analysis of the microbial diversity. This would allow to have a possible interpretation based on the species identified in the MAGs.

For this, several python codes that allows the data analysis are required. First, the `metadata_tax.py` code will create a file that combines the bins name and the taxonomic classification from GTDB-tk.

> Most of these scripts are an modified version of a code written by the M.Sc. Maria Alejandra Soto, a secondee from the project, for the demostrative purposes of the workshop.

```python
import pandas as pd

example_tax_classification = pd.read_csv('/home/dorian.rojas/test/11-gtdbtk/gtdbtk.bac120.summary.tsv', delimiter='\t')
example_tax_classification['user_genome'] = example_tax_classification['user_genome'] + '_filtered_kept_contigs.fasta'
example = example_tax_classification[['user_genome', 'classification']]

# Saving as a new file in the metadata folder 
example.to_csv('/home/dorian.rojas/test/16-metadata/tax.csv', index =False)
```

The end table should look similar to this.

```vim
user_genome,classification
SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta,d__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Treponematales;f__Treponemataceae;g__Treponema_D;s__Treponema_D sp900541945
SRR8555091_bin.2.strict_filtered_kept_contigs.fasta,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__UBA932;g__Cryptobacteroides;s__Cryptobacteroides sp000434935
SRR9988205_bin.1.permissive_filtered_kept_contigs.fasta,d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
SRR9988205_bin.2.permissive_filtered_kept_contigs.fasta,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis
```

Subsequently, the code `metadata_checkm.py` generates a table that combines the taxonomic classification with the quality metrics from CheckM.

```python
import pandas as pd

# gtdbtk classification
sample_tax = pd.read_csv('/home/dorian.rojas/test/11-gtdbtk/gtdbtk.bac120.summary.tsv', sep = '\t')
sample_tax['user_genome'] = sample_tax['user_genome'] + '_filtered_kept_contigs.fasta'

# checkm2
sample_checkm = pd.read_csv('/home/dorian.rojas/test/10-final-MAGs/combined_gunc_and_checkm2.tsv', sep = '\t')
sample_checkm = sample_checkm.rename(columns = {'genome':'user_genome'}) 
sample_checkm['Genome_type'] = 'MAG'

# Combining classification and checkm2 info 
sample_tax = sample_tax.merge(sample_checkm[['user_genome','Genome_type','Completeness','Contamination']], on = 'user_genome', how = 'left')
sample_tax.to_csv('/home/dorian.rojas/test/16-metadata/tax_checkm.csv', index=False)
```

The table is similar to the previously created with adding the CheckM results and a new column named `Genome_type` that indicates the bins represents a MAGs.

```vim
user_genome,classification,closest_genome_reference,closest_genome_reference_radius,closest_genome_taxonomy,closest_genome_ani,closest_genome_af,closest_placement_reference,closest_placement_radius,closest_placement_taxonomy,closest_placement_ani,closest_placement_af,pplacer_taxonomy,classification_method,note,"other_related_references(genome_id,species_name,radius,ANI,AF)",msa_percent,translation_table,red_value,warnings,Genome_type,Completeness,Contamination
SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta,d__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Treponematales;f__Treponemataceae;g__Treponema_D;s__Treponema_D sp900541945,GCA_022010055.1,95.0,d__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Treponematales;f__Treponemataceae;g__Treponema_D;s__Treponema_D sp900541945,98.12,0.752,GCA_022010055.1,95.0,d__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Treponematales;f__Treponemataceae;g__Treponema_D;s__Treponema_D sp900541945,98.12,0.752,d__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Treponematales;f__Treponemataceae;g__Treponema_D;s__,taxonomic classification defined by topology and ANI,topological placement and ANI have congruent species assignments,"GCA_900317625.1, s__Treponema_D sp900317625, 95.0, 85.39, 0.268; GCA_002478955.1, s__Treponema_D sp002478955, 95.0, 85.14, 0.279; GCA_934718065.1, s__Treponema_D sp934718065, 95.0, 84.66, 0.211",79.66,11,,,MAG,86.0,1.87
SRR8555091_bin.2.strict_filtered_kept_contigs.fasta,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__UBA932;g__Cryptobacteroides;s__Cryptobacteroides sp000434935,GCA_000434935.1,95.0,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__UBA932;g__Cryptobacteroides;s__Cryptobacteroides sp000434935,96.83,0.773,GCA_000434935.1,95.0,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__UBA932;g__Cryptobacteroides;s__Cryptobacteroides sp000434935,96.83,0.773,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__UBA932;g__Cryptobacteroides;s__,taxonomic classification defined by topology and ANI,topological placement and ANI have congruent species assignments,,72.67,11,,,MAG,78.16,2.29
SRR9988205_bin.1.permissive_filtered_kept_contigs.fasta,d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli,GCF_003697165.2,95.0,d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli,98.39,0.921,GCF_003697165.2,95.0,d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli,98.39,0.921,d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__,taxonomic classification defined by topology and ANI,topological placement and ANI have congruent species assignments,"GCF_000194175.1, s__Escherichia coli_F, 95.0, 95.66, 0.841; GCF_002965065.1, s__Escherichia sp002965065, 95.0, 94.29, 0.69; GCF_004211955.1, s__Escherichia sp004211955, 95.0, 93.15, 0.728; GCF_005843885.1, s__Escherichia sp005843885, 95.0, 92.99, 0.746; GCF_029876145.1, s__Escherichia ruysiae, 95.0, 92.88, 0.744; GCF_000026225.1, s__Escherichia fergusonii, 95.0, 92.75, 0.565; GCF_011881725.1, s__Escherichia coli_E, 95.0, 92.4, 0.713; GCF_014836715.1, s__Escherichia whittamii, 95.0, 92.08, 0.741; GCF_002900365.1, s__Escherichia marmotae, 95.0, 91.17, 0.711; GCF_000759775.1, s__Escherichia albertii, 95.0, 90.64, 0.611",86.0,11,,,MAG,82.96,0.45
SRR9988205_bin.2.permissive_filtered_kept_contigs.fasta,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis,GCF_000025985.1,95.0,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis,99.01,0.901,GCF_000025985.1,95.0,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis,99.01,0.901,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__,taxonomic classification defined by topology and ANI,topological placement and ANI have congruent species assignments,"GCA_947646015.1, s__Bacteroides sp947646015, 95.0, 97.98, 0.349; GCF_019583405.1, s__Bacteroides fragilis_B, 95.0, 91.9, 0.692; GCF_014288095.1, s__Bacteroides hominis, 95.0, 88.11, 0.582",90.21,11,,,MAG,95.32,0.27
```

Now, it is required to use these files to create table that combines both the abundance estimated from inStrain and the taxonomic classification. First, create a directory named `16-metadata/instrain_copies` and copy the `genome_info.tsv` generated from inStrain into this folder. This will allow to edit their contents without affecting the original files.

```bash
(base) [dorian.rojas@accessnode test]$ mkdir 16-metadata/instrain_copies/
(base) [dorian.rojas@accessnode test]$ cp 15-instrain/*/output/*_genome_info.tsv 16-metadata/instrain_copies/
(base) [dorian.rojas@accessnode test]$ ll -h 16-metadata/*
-rw-rw-r-- 1 dorian.rojas dorian.rojas 4,4K feb 26 12:29 16-metadata/tax_checkm.csv
-rw-rw-r-- 1 dorian.rojas dorian.rojas  742 feb 26 12:39 16-metadata/tax.csv

16-metadata/instrain_copies:
total 8,0K
-rw-rw-r-- 1 dorian.rojas dorian.rojas 1,3K feb 26 12:42 SRR8555091_genome_info.tsv
-rw-rw-r-- 1 dorian.rojas dorian.rojas 1,2K feb 26 12:42 SRR9988205_genome_info.tsv
```

The python code `abundance_tax_filter.py` will generate one abundance table per sample with that include the taxonomic classification into the inStrain results.

This code could allow the filtering of the identified bins based in the breadth value. As mentiones above, certain value could represent that these bins are not to be considered bacteria but rather certain genetic elements or mismapping. Therefore, depending on the study objective, researchers aim to remove bins with low abundance. There is no 'standard' for the breadth value, so it is strongly dependent on the overall results of the metagenomic analysis. However, in this case, considering the small number of bins analyzed, this filtering step won't be performed.

```python
import pandas as pd
import os

## Extract taxonomy information for the GTDBtk results

metadata_file_paths = [
    '/home/dorian.rojas/test/16-metadata/tax_checkm.csv',
    '/home/dorian.rojas/test/16-metadata/tax.csv'
]

# Step 1: Create a dictionary of user_genome and taxonomic classification
classification_dict = {}

for file in metadata_file_paths:
    taxonomy_data = pd.read_csv(file)
    taxonomy_data['user_genome'] = taxonomy_data['user_genome'].str.replace('_filtered_kept_contigs.fasta', '', regex=False)
    classification_dict.update(dict(zip(taxonomy_data['user_genome'], taxonomy_data['classification'])))

## Add taxonomy information to abundance estimation results (InStrain output) and filter rows with genome coverage lower than 0.5 (optional)

# List of .csv files to process (original inStrain outputs are in another folder (/home/dorian.rojas/test/15-instrain). The tables that will be updated are copies)
abundance_tables_parent_dir = '/home/dorian.rojas/test/16-metadata/instrain_copies/' 
output_dir = '/home/dorian.rojas/test/16-metadata/'  

#Step 2: Iterate through .tsv files in the parent directory and update them 
for root, _, files in os.walk(abundance_tables_parent_dir):
    for file in files:
        if file.endswith('.tsv'):
            tsv_path = os.path.join(root, file)
            abundance_tsv_data = pd.read_csv(tsv_path, sep='\t')
            
            # Clean user_genome values in the .tsv file
            abundance_tsv_data = abundance_tsv_data.rename(columns = {'genome':'user_genome'}) 
            abundance_tsv_data['user_genome'] = abundance_tsv_data['user_genome'].str.replace('_filtered_kept_contigs.fasta', '', regex=False)
            abundance_tsv_data['user_genome'] = abundance_tsv_data['user_genome'].str.replace('.fasta', '', regex=False)
            
            # Map the classification column to the user_genome column
            abundance_tsv_data['classification'] = abundance_tsv_data['user_genome'].map(classification_dict)
            
            # Remove rows where the 'breadth' column is lower than 0.5 (optional)
            # abundance_tsv_data =  abundance_tsv_data[abundance_tsv_data['breadth'] >= 0.5]

            # Save the updated .tsv file
            output_path = os.path.join(output_dir, file)
            abundance_tsv_data.to_csv(output_path, sep='\t', index=False)
            print(f"Updated file saved to: {output_path}")
```

The output from this python should look similar to the inStrain with an additional column named `classification` at the end.

```bash
(base) [dorian.rojas@accessnode test]$ ls 16-metadata/
instrain_copies  SRR8555091_genome_info.tsv  SRR9988205_genome_info.tsv  tax_checkm.csv  tax.csv
(base) [dorian.rojas@accessnode test]$ cat 16-metadata/SRR9988205_genome_info.tsv
user_genome     coverage        breadth nucl_diversity  length  true_scaffolds  detected_scaffolds      coverage_median     coverage_std    coverage_SEM    breadth_minCov  breadth_expected        nucl_diversity_rarefied     conANI_reference        popANI_reference        iRep    iRep_GC_corrected       linked_SNV_count   SNV_distance_mean        r2_mean d_prime_mean    consensus_divergent_sites       population_divergent_sites SNS_count        SNV_count       filtered_read_pair_count        reads_unfiltered_pairs  reads_mean_PID  divergent_site_count        reads_unfiltered_reads  classification
SRR9988205_bin.1.permissive     5.496498350383906       0.994117863265243       0.0008198111233099      3781789     781     781     5       2.6698752312885534      0.001402173454713       0.6071142520114158      0.9921982564903672  0.0     0.9999629787393434      0.999998257823028               False   203     26.88177339901478   0.627394623981047       0.9875701459034792      85      4       3       459     92969   94802   0.9982983789351916  462     191829  d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
SRR9988205_bin.2.permissive     9.223122918012026       0.9982518369063236      0.0005069402542589      4456678     545     545     9       4.186153449138965       0.0020076430886163      0.8682996617660059      0.9997095321144326  0.0     0.9999790683179094      0.9999994831683434      1.9121057585254528      False   52530.02285714285714        0.5179457671047509      0.9642752701766494      81      2       2       676     190590      194320  0.998575959217324       678     390820  d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis
```

To continue, it is require to consolidate the abundances of each bin into a single table. The following code (`consolidate_abundance.py`) creates an abundance table with each MAGs as columns and sample as rows. Notice the `filtered_read_pair_count` is the extracted column used to evaluate the abundance of each bin.

```python
import pandas as pd
import glob
import os

# Define the path to your inStrain otput (genome_info.tsv files)
file_paths = glob.glob('/home/dorian.rojas/test/16-metadata/*_genome_info.tsv')

# Create an empty DataFrame to store all the abundance results
final_df = pd.DataFrame()

# Process each file
for file_path in file_paths:
    file_name = os.path.basename(file_path).replace('_genome_info.tsv', '')
    df = pd.read_csv(file_path, sep='\t')

# Create a dictionary with 'filtered_read_pair_count' values, even if they are 0
    data = {
        'file_name': [file_name],
        **{genome: [filtered_count] for genome, filtered_count in zip(df['user_genome'], df['filtered_read_pair_count'])}
    }

    # Convert the dictionary to a DataFrame and concatenate it to the final DataFrame
    result = pd.DataFrame(data)
    final_df = pd.concat([final_df, result], ignore_index=True)

# Fill NaN values with 0 to ensure all genomes have a value
final_df.fillna(0, inplace=True)

# Save the combined DataFrame to a CSV file
final_df.to_csv('/home/dorian.rojas/test/16-metadata/final_combined_abundance_output.csv', index=False)
print("Combined CSV file 'final_combined_abundance_output.csv' has been created.")
```

In this case, due to the small number of bins and samples, the output is short.

```vim
file_name,SRR8555091_bin.1.permissive,SRR8555091_bin.2.strict,SRR9988205_bin.1.permissive,SRR9988205_bin.2.permissive
SRR8555091,55394.0,56537.0,0.0,0.0
SRR9988205,0.0,0.0,92969.0,190590.0
```

Finally, prior the diversity analysis, a final table with the input for the R packages needs to be creates. For this, a final python code will combine both the `tax.csv` and the `tax_checkm.csv` into an acceptable format for the R packages.

Prior to run this code, a `MAGs_list.txt` file needs to be created. This file is made using the command `awk '{print $2}' path/to/gtdb-input.txt > output/MAGs_list.txt`. For instance:

```bash
(base) [dorian.rojas@accessnode test]$ awk '{print $2}' 10-final-MAGs/gtdb-input.txt > 16-metadata/MAGs_list.txt
(base) [dorian.rojas@accessnode test]$ cat 16-metadata/MAGs_list.txt
SRR8555091_bin.1.permissive
SRR8555091_bin.2.strict
SRR9988205_bin.1.permissive
SRR9988205_bin.2.permissive
```

This file contains the genome name of the final bins present in the pangenome. Finally, the `metadata_for_R.py` code will generate the required table.

```python
import pandas as pd

# File paths
metadata_file_paths = [
    '/home/dorian.rojas/test/16-metadata/tax_checkm.csv',
    '/home/dorian.rojas/test/16-metadata/tax.csv'
]

genome_list_path = '/home/dorian.rojas/test/16-metadata/MAGs_list.txt'  
output_csv_path = '/home/dorian.rojas/test/16-metadata/all_MAGs_taxonomy.csv' 

# Step 1: Read the genome list
genome_list = []
with open(genome_list_path, 'r') as file:
    genome_list = [line.strip() for line in file]

# Step 2: Create a dictionary of user_genome -> classification
classification_dict = {}

for file in metadata_file_paths:
    taxonomy_data = pd.read_csv(file)
    taxonomy_data['user_genome'] = taxonomy_data['user_genome'].str.replace('_filtered_kept_contigs.fasta', '', regex=False)
    taxonomy_data['user_genome'] = taxonomy_data['user_genome'].str.replace('.fa', '', regex=False)
    classification_dict.update(dict(zip(taxonomy_data['user_genome'], taxonomy_data['classification'])))

# Step 3: create a dictionary to store the genome's name (User_Genome) and the extracted taxonomic levels 
columns = ['User_Genome', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
output_data = []

for genome in genome_list:
    classification = classification_dict.get(genome, None)
    if classification:
        taxonomy_levels = classification.split(';')
        taxonomy_dict = {
            'User_Genome': genome,
            'Domain': next((x.replace('d__', '') for x in taxonomy_levels if x.startswith('d__')), ''),
            'Phylum': next((x.replace('p__', '') for x in taxonomy_levels if x.startswith('p__')), ''),
            'Class': next((x.replace('c__', '') for x in taxonomy_levels if x.startswith('c__')), ''),
            'Order': next((x.replace('o__', '') for x in taxonomy_levels if x.startswith('o__')), ''),
            'Family': next((x.replace('f__', '') for x in taxonomy_levels if x.startswith('f__')), ''),
            'Genus': next((x.replace('g__', '') for x in taxonomy_levels if x.startswith('g__')), ''),
            'Species': next((x.replace('s__', '') for x in taxonomy_levels if x.startswith('s__')), '')
        }
        output_data.append(taxonomy_dict)  # Add the dictionary to the output_data list

# Step 4: Save the output to a .csv file
output_df = pd.DataFrame(output_data, columns=columns)  # converts the list of dictionaries (output_data) into a DataFrame.
output_df.to_csv(output_csv_path, index=False)

print(f"Taxonomy data has been saved to {output_csv_path}.")
```

The output is named `all_MAGs_taxonomy.csv` and contains the column names refering to the sample and each taxonomic classification level and rows as the data respective for each column.

```vim
User_Genome,Domain,Phylum,Class,Order,Family,Genus,Species
SRR8555091_bin.1.permissive,Bacteria,Spirochaetota,Spirochaetia,Treponematales,Treponemataceae,Treponema_D,Treponema_D sp900541945
SRR8555091_bin.2.strict,Bacteria,Bacteroidota,Bacteroidia,Bacteroidales,UBA932,Cryptobacteroides,Cryptobacteroides sp000434935
SRR9988205_bin.1.permissive,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia coli
SRR9988205_bin.2.permissive,Bacteria,Bacteroidota,Bacteroidia,Bacteroidales,Bacteroidaceae,Bacteroides,Bacteroides fragilis
```

The results from these data maniputations will be used to conduct the diversity analysis in the R code.

### Coding the phyloseq and microviz R packages

The diversity analysis of the microbiome requires three different files, an Operational Taxonomic Units (OTUs) table, a taxonomy table, and a metadata table. From the previous analysis, the `final_combined_abundance_output.csv` represents the OTUs table with the bins and respective abundance. Meanwhile, the `all_MAGs_taxonomy.csv` is the tax table.

The metadata table is a desired feature of most metagenomic analysis for the interpretation of microbial composition. Commonly, these tables would be composed of a column with the name of each of the samples and information about the individual (e.g. age, weigth, BMI, height, country). These serves as categorical factors to evaluate the microbial composotion and the abundances of each bin.

This information is missing or it is poor in our samples. Hence, this part of the workshop is demonstrative based on a similar analysis conducted with a greater amount of samples and extensive metadata.

## Code templates

**`bowtie-index.slurm` file:**

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="index"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL

cd /home/dorian.rojas/test

CTN_PATH=/opt/ohpc/pub/containers/BIO/


mkdir -p 14-bowtie/1-index


#Bowtie2 does not recognize the use of \. Command has to be written in a single line
$CTN_PATH/bowtie2-2.5.4.sif bowtie2-build --threads 64 10-final-MAGs/bins/pangenome.fasta 14-bowtie/1-index/pangenome


date
time
```

**`bowtie.slurm` file:**

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="mapping"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL
start=`date +%s`

cd /home/dorian.rojas/test

CTN_PATH=/opt/ohpc/pub/containers/BIO/
output=/home/dorian.rojas/test/14-bowtie
index=/home/dorian.rojas/test/14-bowtie/1-index

mkdir -p 14-bowtie/

for sample in $@; do

#Bowtie2 does not recognize the use of \. Command has to be written in a single line
#Mapping
$CTN_PATH/bowtie2-2.5.4.sif bowtie2 --threads 64 -x $index/pangenome -1 1-data/${sample}_sub_1.fastq -2 1-data/${sample}_sub_2.fastq -S $output/${sample}.sam

#SAM to BAM
$CTN_PATH/samtools-1.21.sif samtools view -F 4 -bS $output/${sample}.sam > $output/${sample}-raw.bam

#Sort and index
$CTN_PATH/samtools-1.21.sif samtools sort -@ 64 $output/${sample}-raw.bam -o $output/${sample}.sorted.bam
$CTN_PATH/samtools-1.21.sif samtools index $output/${sample}.sorted.bam

rm $output/${sample}-raw.bam

done

date
end=`date +%s`
runtime=$((end-start))
```

**`parse_stb.py` file:**

```python
#!/usr/bin/env python

import os
import sys
import gzip
import textwrap
import argparse

def extract_bins(fasta, stb_file, out_base):
    if (out_base != '') & (out_base[-1] != '/'):
        out_base += '_'

    # Make scaffold to bin dictionary
    stb = {}
    with open(stb_file) as stb_reader:
        for line in stb_reader:
            line = line.strip()

            if line.startswith('#') or line.startswith('scaffold_name'):
                continue

            scaffold, bin = line.split('\t')[:2]
            stb[scaffold.strip()] = bin.strip()

    # Parse the FASTA file and write to bins
    if fasta.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'

    with open_func(fasta, mode) as handle:
        record_id = ''
        record_seq = ''
        for line in handle:
            if line.startswith('>'):
                if record_id and record_id in stb:
                    bin = stb[record_id]
                    with open(f"{out_base}{bin}.fa", 'a') as nw:
                        nw.write(f">{record_id}\n{record_seq}\n")
                record_id = line[1:].strip().split(' ')[0]
                record_seq = ''
            else:
                record_seq += line.strip()

        # Write the last record
        if record_id and record_id in stb:
            bin = stb[record_id]
            with open(f"{out_base}{bin}.fa", 'a') as nw:
                nw.write(f">{record_id}\n{record_seq}\n")

import os
import gzip

def gen_stb(fastas):
    if ((len(fastas) == 1) & (not fastas[0].endswith('.gz'))):
        # See if this is a text file, not a fasta file
        text_list = True
        genomes = []
        with open(fastas[0], 'r') as o:
            for line in o.readlines():
                if line.startswith('>'):
                    text_list = False
                    break
                else:
                    genomes.append(line.strip())
        if text_list:
            print("Treating .fasta input as list")
            fastas = genomes

    stb = {}
    for fasta in fastas:
        bin = os.path.basename(fasta)
        if fasta.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        with open_func(fasta, mode) as handle:
            for line in handle:
                if line.startswith('>'):
                    id = line[1:].strip().split(' ')[0]
                    stb[id] = bin

    return stb


def print_stb(stb, output):
    file = open(output,'w')
    for key in sorted(stb, key=stb.get):
        file.write("{0}\t{1}\n".format(key,stb[key]))
    file.close()


class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=textwrap.dedent('''\
         \n
         The program has two uses related to scaffold to bin (.stb) files.
         .stb files should be tab-separated, with no header, and two columns: scaffold and bin

         Use 1) Pass a list of genomes to generate a .stb file.

         Example:
         parse_stb.py --reverse -f dereplicate_genomes/* -o representitve_genomes.stb

         Use 2) Pass a single .fasta file and a scaffold to bin file (.stb) to generate a number of
         fasta files based on the .stb file.

         Example:
         parse_stb.py -f concat_genomes.fasta -s scaffold_to_bin.tsv -o genomeList_1
         '''))

    parser.add_argument('-s','--stb',help='scaffold to bin file')
    parser.add_argument('-f','--fasta',help='fasta file to extract scaffolds from. Will treat as compressed if ends in .gz. This can also be a single text file with a genome on each line',nargs='*')
    parser.add_argument('-o','--output',help='output base name', default = '')

    parser.add_argument('--reverse',help='generate a stb from a list of genomes',\
                        action = "store_true")
    args = parser.parse_args()

    if args.reverse == False:
        if len(args.fasta) != 1:
            print("must give one and only one fasta file")
            sys.exit()
        extract_bins(args.fasta[0], args.stb, args.output)

    if args.reverse == True:
        stb = gen_stb(args.fasta)
        print_stb(stb, args.output)
```

**`instrain.slurm` file:**

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="instrain"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL
start=`date +%s`

cd /home/dorian.rojas/test

CTN_PATH=/opt/ohpc/pub/containers/BIO/
bam=/home/dorian.rojas/test/14-bowtie/
pangenome=/home/dorian.rojas/test/10-final-MAGs/bins
output=/home/dorian.rojas/test/15-instrain

mkdir -p 15-instrain/$sample

for sample in $@; do

$CTN_PATH/instrain-1.9.0.sif inStrain profile $bam/${sample}.sorted.bam \
        $pangenome/pangenome.fasta -o $output/$sample -p 64 \
        -s $output/representative_genomes.stb \
        --database_mode --skip_plot_generation

done

date
end=`date +%s`
runtime=$((end-start))
echo 'Runtime was ' $runtime
```

**`metadata_tax.py` file:**

```python
import pandas as pd

example_tax_classification = pd.read_csv('/home/dorian.rojas/test/11-gtdbtk/gtdbtk.bac120.summary.tsv', delimiter='\t')
example_tax_classification['user_genome'] = example_tax_classification['user_genome'] + '_filtered_kept_contigs.fasta'
example = example_tax_classification[['user_genome', 'classification']]

# Saving as a new file in the metadata folder 
example.to_csv('/home/dorian.rojas/test/16-metadata/tax.csv', index =False)
```

**`metadata_checkm.py` file:**

```python
import pandas as pd

# gtdbtk classification
sample_tax = pd.read_csv('/home/dorian.rojas/test/11-gtdbtk/gtdbtk.bac120.summary.tsv', sep = '\t')
sample_tax['user_genome'] = sample_tax['user_genome'] + '_filtered_kept_contigs.fasta'

# checkm2
sample_checkm = pd.read_csv('/home/dorian.rojas/test/10-final-MAGs/combined_gunc_and_checkm2.tsv', sep = '\t')
sample_checkm = sample_checkm.rename(columns = {'genome':'user_genome'}) 
sample_checkm['Genome_type'] = 'MAG'

# Combining classification and checkm2 info 
sample_tax = sample_tax.merge(sample_checkm[['user_genome','Genome_type','Completeness','Contamination']], on = 'user_genome', how = 'left')
sample_tax.to_csv('/home/dorian.rojas/test/16-metadata/tax_checkm.csv', index=False)
```

**`abundance_tax_filter.py` file:**

```python
import pandas as pd
import os

## Extract taxonomy information for the GTDBtk results

metadata_file_paths = [
    '/home/dorian.rojas/test/16-metadata/tax_checkm.csv',
    '/home/dorian.rojas/test/16-metadata/tax.csv'
]

# Step 1: Create a dictionary of user_genome and taxonomic classification
classification_dict = {}

for file in metadata_file_paths:
    taxonomy_data = pd.read_csv(file)
    taxonomy_data['user_genome'] = taxonomy_data['user_genome'].str.replace('_filtered_kept_contigs.fasta', '', regex=False)
    classification_dict.update(dict(zip(taxonomy_data['user_genome'], taxonomy_data['classification'])))

## Add taxonomy information to abundance estimation results (InStrain output) and filter rows with genome coverage lower than 0.5 (optional)

# List of .csv files to process (original inStrain outputs are in another folder (/home/dorian.rojas/test/15-instrain). The tables that will be updated are copies)
abundance_tables_parent_dir = '/home/dorian.rojas/test/16-metadata/instrain_copies/' 
output_dir = '/home/dorian.rojas/test/16-metadata/'  

#Step 2: Iterate through .tsv files in the parent directory and update them 
for root, _, files in os.walk(abundance_tables_parent_dir):
    for file in files:
        if file.endswith('.tsv'):
            tsv_path = os.path.join(root, file)
            abundance_tsv_data = pd.read_csv(tsv_path, sep='\t')
            
            # Clean user_genome values in the .tsv file
            abundance_tsv_data = abundance_tsv_data.rename(columns = {'genome':'user_genome'}) 
            abundance_tsv_data['user_genome'] = abundance_tsv_data['user_genome'].str.replace('_filtered_kept_contigs.fasta', '', regex=False)
            abundance_tsv_data['user_genome'] = abundance_tsv_data['user_genome'].str.replace('.fasta', '', regex=False)
            
            # Map the classification column to the user_genome column
            abundance_tsv_data['classification'] = abundance_tsv_data['user_genome'].map(classification_dict)
            
            # Remove rows where the 'breadth' column is lower than 0.5 (optional)
            # abundance_tsv_data =  abundance_tsv_data[abundance_tsv_data['breadth'] >= 0.5]

            # Save the updated .tsv file
            output_path = os.path.join(output_dir, file)
            abundance_tsv_data.to_csv(output_path, sep='\t', index=False)
            print(f"Updated file saved to: {output_path}")
```

**`consolidate_abundance.py` file:**

```python
import pandas as pd
import glob
import os

# Define the path to your inStrain otput (genome_info.tsv files)
file_paths = glob.glob('/home/dorian.rojas/test/16-metadata/*_genome_info.tsv')

# Create an empty DataFrame to store all the abundance results
final_df = pd.DataFrame()

# Process each file
for file_path in file_paths:
    file_name = os.path.basename(file_path).replace('_genome_info.tsv', '')
    df = pd.read_csv(file_path, sep='\t')

# Create a dictionary with 'filtered_read_pair_count' values, even if they are 0
    data = {
        'file_name': [file_name],
        **{genome: [filtered_count] for genome, filtered_count in zip(df['user_genome'], df['filtered_read_pair_count'])}
    }

    # Convert the dictionary to a DataFrame and concatenate it to the final DataFrame
    result = pd.DataFrame(data)
    final_df = pd.concat([final_df, result], ignore_index=True)

# Fill NaN values with 0 to ensure all genomes have a value
final_df.fillna(0, inplace=True)

# Save the combined DataFrame to a CSV file
final_df.to_csv('/home/dorian.rojas/test/16-metadata/final_combined_abundance_output.csv', index=False)
print("Combined CSV file 'final_combined_abundance_output.csv' has been created.")
```

**`metadata_for_R.py` file:**

```python

```
