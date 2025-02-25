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

[]

### Coding the phyloseq and microviz R packages

[]

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
