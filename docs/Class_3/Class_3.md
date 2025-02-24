# Class 3: MAGs quality evaluation

- - - -

## MAGs quality parameters (MDMcleaner, GUNC, CheckM2)

The Minimum Information about a MAG (MIMAG) required for their classification have been defined ([Bowers et al. 2017](https://www.nature.com/articles/nbt.3893)). The complete parameters to delimit MAGs include assembly quality, genome completeness, and measures of contamination. In addition, these characteristics favor the determination of MAGs as finished and draft genomes.

Criterion|Assembly annotation|Completion|Contamination
-------|-------------------|----------|-------------
Finished|Single contigous sequencing without gaps or ambiguities with consensus error rate equivalent to Q50 or better|>90%|<5%
High-quality draft|Gaps in repetitive regions. Presence of 23S, 16S, and 5S rRNA gene, or >18 tRNAs|>90%|<5%
Medium-quality draft|Many fragments without gene annotation, but basic assembly stats|>50%|<10%
Low-quality draft|Many fragments without gene annotation, but basic assembly stats|<50%|<10%

>The basic assembly stats of MAGs are similar to those of genomic assemblies: N50, L50, largest contig, number of contigs, total size, reads percentage mapped to assembly, and number of predicted genes ([Bowers et al. 2017](https://www.nature.com/articles/nbt.3893)). However, other params might be reported in addition.

The quality of reassembled bins were evaluated using the MetaWRAP wapper. Each of the modules included the quality report by CheckM. After the publication of MIMAG, CheckM became the _de facto_ gold standard for determining genome quality. However, it holds certain challenges. For instance, highly fragmented MAGs might not present a recognizable conserved marked genes. This difficults the discrimination between correctly assigned contigs and contamination. In this regard, some tools provide more accurate stats than CheckM. Thus, the use of various quality control softwares is recommended to ensure correct and complete scores.

To complement the quality assessment, this class focus on running several tools for quality check. These software include: MDMcleaner, GUNC, and CheckM2.

[MDMcleaner](https://github.com/KIT-IBG-5/mdmcleaner) (Microbial Dark Matter cleaner) is a python workflow that allows to dectect and remove contamination with minimum error. Properly, it performs a quality check as well as improving the bins quality by removing incorrectly assigned reads. For this, marker genes are extracted from the contigs (e.g. small subunit rRNA genes, large subunit rRNA genes, universal protein coding marker genes, total coding sequences, and tRNA genes). Non-coding sequences are discarded due to their mostly abundance in eukaryotic contaminants. The recollected genes are aligned against a comprised database based on GTBD, SILVA, and RefSeq up-to-date databases. This results in a set of genomes per species that avoids overrespresentation of taxa. The preliminar taxonomy is assignated using a Least Common Ancestor (LCA) approach on the respective blast-hits per genes. Later, the final classification is defined in a similar way, but using contigs blast-hits. The over-classification is avoided by pruning those ranks that are not supported by the common 16S rRNA and protein coding genes cut-off values.

> This LCA classification might have ambiguities and contamination due to some high taxonomic ranks (e.g. domain, phylum). Therefore, an additional workflow analyzes this conflict and performs a downstream evaluation for reference databases contaminations.

Later, a weighted majority consensus of the overall genome classification is extracted and each contig is assigned a trustworthiness score. This score is based on deviation from the taxonomic assignation, marker genes, and corresponding alignment identities that were involved. It ranges from lowest trustworthiness `0` to highest trustworthiness `10`. The workflow output is a detailed report on the contigs scores, and individual `fasta` files with the separated sequences. These files indicate those kept, deleted, and possibly-further-evaluated contigs.

Briefly, MDMcleaner allows a MAG evaluation based on varios databases of conserved marker genes from both cultured and uncultured taxa. The further analysis uses aligment scores from individual genes and general contigs for the calculation of a trustworthiness score that is used to discriminate the sequences. The strategy has high sensitivity to contaminants in fragmented genomes and underrepresented taxa ([Vollmers et al. 2022](https://academic.oup.com/nar/article/50/13/e76/6583244?login=false)).

[GUNC](https://github.com/grp-bork/gunc) (Genome UNClutterer) is a tool that allows the detection of chrimerism and contamination from prokaryotic genomes. These result from misassembly and mis-binning, although the latter is the most common source of errors. For this, GUNC identifies genes using prodigal and mapped through diamond to the GUNC reference database (based on proGenomes). This is used to calculate the GUNC scores and generate interactive Sankey diagrams that favors the visualization of taxonomic composition. GUNC evaluates the chimerism and representation in the difference taxonomic levels to estimate Clase Separations Scores (CSS) and Reference Representation Scores (RRS). CSS are high in case of gene classification to different groups and RRS are high in case of close mapping to the GUNC reference.

Overall, GUNC uses a entropy-based measurement of lineage homogeneity to calculate chimerism. The output files are divided in two `.tsv` files that summarises the scores with the highest CSS and all the taxonomic levels. If the genome has a CSS value <= 0.45, it passes the GUNC analysis indicating that might not be chimeric (or it cannot be detected). This value is based on a benchmarking study using simulated genomes ([Orakov et al. 2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02393-0)).

[CheckM2](https://github.com/chklovski/CheckM2) is an update from CheckM and works in very similar ways. However, CheckM2 analysis involves a machine learning model universally trained to assess the completeness and contamination of bins ([Chklosvski et al. 2023](https://www.nature.com/articles/s41592-023-01940-w)). CheckM2 is highly accurate on microorganisms with reduced genomes or with only a few genomic representative. In this regard, to predict contamination and completion two machine learning models are employed. The general gradient boost model is efficient under poorly represented microbes in the GenBank or RefSeq. Although, this approach is applied for all inputs regardless of taxonomic novelty. Moreover, the specific neural network model is accurate when predicting the score of microorganisms closely related to the training set. To determine the appropiate model, a cosine similarity calculation is used by deafult. However, a particular model or the result from both can be indicated throught command options.

In this regard, CheckM2 does a similar job to CheckM but with improved and more accurate results.

### Running mdmcleaner

> Prior running the MAGs processing, the MDMcleaner database must be downloaded and directory set in the configuration. This might take a while (>13 hours). Thus, it was previously installed and set for usage.

The workflow was designed to accept both Single Amplified Genomes (SAGs) and MAGs. Therefore, MDMcleaner holds multiple commands. MDMcleaner has to be individually installed under the base mamba environment. For this, create a `.sh` file similar to the one used to install metaWRAP. This can be named as `MDMcleaner_install.sh`. Copy and paste the code below and run using `sh MDMcleaner_install.sh`. Make sure to be in the `(base)` environment prior to compile the shell script (use `conda deactivate` if you have any environment activated).

```vim
#!/bin/bash
conda create -n mdmcleaner
conda install -c bioconda -n mdmcleaner
conda activate mdmcleaner
```

This short script will create a new environment named `mdmcleaner` (`conda create...`) and install the MDMcleaner software extracted from the bioconda channel (`conda install...`). Finally the environment will be automatically activated. For any further use, by running `conda activate mdmcleaner` in the console, this environment will work correctly.

Run the basic mdmclenaer command and check the options.

```bash
(base) [dorian.rojas@accessnode ~]$ mdmcleaner -h
usage: mdmcleaner [-h] {clean,makedb,get_markers,completeness,acc2taxpath,refdb_contams,set_configs,show_configs,check_dependencies,version} ...

MDMcleaner pipeline v0.8.7 for decontaminating and classifying microbial dark matter MAGs and SAGs

positional arguments:
  {clean,makedb,get_markers,completeness,acc2taxpath,refdb_contams,set_configs,show_configs,check_dependencies,version}
    clean               classify and filter contigs from microbial dark matter MAGs and SAGs
    makedb              Download and create MDMcleaner database
    get_markers         extracts protein coding and/or rRNA gene sequences from input genome(s)
    completeness        estimate completeness (roughly based on presence of universally required tRNA types). Results are printed directly to stdout
    acc2taxpath         Get full taxonomic path assorciated with a specific acession number
    refdb_contams       EXPERIMENTAL: evaluate potentiel refDB-contaminations
    set_configs         setting or changing settings in config files
    show_configs        check settings in config files (with highest ranking source for each setting)
    check_dependencies  checks if all dependencies for MDMcleaner are being met
    version             show version info and exit

options:
  -h, --help            show this help message and exit
```

The command `set_configs` is required to indicate the path towards the database. The dabatase is downloaded in the public folder `/home/public/DB` and it is named `gtdb`. However, the tools only required the folder CONTAINING the database. This needs to be set each time the software is used within the for loop before the actual command. Check the options.

```bash
(base) [dorian.rojas@accessnode test]$ mdmcleaner set_configs -h
usage: mdmcleaner set_configs [-h] [-s {local,global}] [--blastp BLASTP] [--blastn BLASTN]
                              [--diamond DIAMOND] [--barrnap BARRNAP] [--hmmsearch HMMSEARCH]
                              [--aragorn ARAGORN] [--db_basedir DB_BASEDIR] [--threads THREADS]

options:
  -h, --help            show this help message and exit
  -s {local,global}, --scope {local,global}
                        change settings in local or global config file. 'global' likely require admin
                        privileges. 'local' will modify or create a mdmcleaner.config file in the current
                        working directory. default = 'local'
  --blastp BLASTP       path to blastp binaries (if not in PATH)
  --blastn BLASTN       path to blastn binaries (if not in PATH)
  --diamond DIAMOND     path to diamond binaries (if not in PATH)
  --barrnap BARRNAP     path to barrnap binaries (if not in PATH)
  --hmmsearch HMMSEARCH
                        path to hmmsearch binaries (if not in PATH)
  --aragorn ARAGORN     path to aragorn binaries (if not in PATH)
  --db_basedir DB_BASEDIR
                        path to basedirectory for reference database
  --threads THREADS     threads to use by default
```

Here, the option `--db_basedir` is the one of interest for setting the database path.

The command `clean` is the one performing the MAGs analysis. It has a simple usage.

```bash
(base) [dorian.rojas@accessnode ~]$ mdmcleaner clean -h
usage: mdmcleaner clean [-h] -i INPUT_FASTAS [INPUT_FASTAS ...] [-o OUTPUT_FOLDER] [--outblacklist OUTBLACKLIST] [-v] [-c CONFIGFILE] [-t THREADS] [-f]
                        [--overview_files_basename OVERVIEW_BASENAME] [-b BLACKLISTFILE] [--no_filterfasta] [--ignore_default_blacklist] [--fast_run]

options:
  -h, --help            show this help message and exit
  -i INPUT_FASTAS [INPUT_FASTAS ...], --input_fastas INPUT_FASTAS [INPUT_FASTAS ...]
                        input fastas of genomes and/or bins
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        output-folder for MDMcleaner results. Default = 'mdmcleaner_output'
  --outblacklist OUTBLACKLIST
                        Outputfile for new blacklist additions. If a preexisting file is selected, additions will be appended to end of that file
  -v, --version         show programs version number and exit
  -c CONFIGFILE, --config CONFIGFILE
                        provide a local config file with basic settings (such as the location of database-files). default: looks for config files named 'mdmcleaner.config' in
                        current working directory. settings in the local config file will override settings in the global config file
                        '/home/dorian.rojas/bin/miniforge3/lib/python3.12/site-packages/mdmcleaner/mdmcleaner.config'
  -t THREADS, --threads THREADS
                        Number of threads to use. (default = use setting from config file)
  -f, --force           Force reclassification of pre-existing blast-results
  --overview_files_basename OVERVIEW_BASENAME
                        basename for overviewfiles (default="overview"
  -b BLACKLISTFILE, --blacklistfile BLACKLISTFILE
                        File listing reference-DB sequence-names that should be ignored during blast-analyses (e.g. known refDB-contaminations...
  --no_filterfasta      Do not write filtered contigs to final output fastas (Default = False)
  --ignore_default_blacklist
                        Ignore the default blacklist (Default = False)
  --fast_run            skip detailed analyses of potential reference ambiguities (runs may be faster but also classification may be less exact, and potential reference database
                        contaminations will not be verified)
```

MDMcleaner outputs overview files of the complete process and the `mdmcleane.config` file where the database and other commands are indicated. These are created in the working directory. Therefore, to avoid each of the samples to overwrite the results, the working directory must change. This would also require to direct the `clean` command to the complete path of the input files. You can set this with a variable to make the code shorter. The for loop should look something simialr to the example below:

```vim
cd /home/dorian.rojas/test
WD=/home/dorian.rojas/test

for sample in $@; do

mkdir 7-mdmcleaner/${sample}_sub
cd 7-mdmcleaner/${sample}_sub

mdmcleaner set_configs --db_basedir /home/public/DB
mdmcleaner clean -i $WD/6-reassemble/${sample}_sub/reassembled_bins/*.fa \
        -o $WD/7-mdmcleaner/${sample}_sub/ -t 64

cd ../../

done

date
time
```

MDMcleaner can analyze single genomes and the input command recognizes a single file. Hence, a wildcard is required (see `*.fa` in the example).

Write a code that would clean the reassembled bins and save the results under a folder named `7-mdmcleaner`, use the maximum of threads possible (64).

> If the job takes a while running, the results in the common repository `/home/public/met-workshop`

As mentioned above, the output from MDMcleaner are very straightforward. It involves a subfolder for each analyzed bin, the overview files, configs file, and several documents with in each bin subfolder.

```bash
(base) [dorian.rojas@accessnode test]$ ls 7-mdmcleaner/SRR8555091_sub/*
7-mdmcleaner/SRR8555091_sub/mdmcleaner.config
7-mdmcleaner/SRR8555091_sub/new_blacklist_additions.tsv
7-mdmcleaner/SRR8555091_sub/overview_all_before_cleanup.tsv
7-mdmcleaner/SRR8555091_sub/overview_errorlist.txt

7-mdmcleaner/SRR8555091_sub/bin.1.permissive:
arch_pfam.hmm.domtblout
arch_tigr.hmm.domtblout
bact_pfam.hmm.domtblout
bact_tigr.hmm.domtblout
bin.1.permissive_filtered_kept_contigs.fasta.gz
bin.1.permissive_filtered_need_evaluation_high.fasta.gz
bin.1.permissive_filtered_need_evaluation_low.fasta.gz
bin.1.permissive_filtered_removed.fasta.gz
bin.1.permissive_kronainput.tsv
bin.1.permissive_rRNA_tsu_rRNA.fasta
bin.1.permissive_totalprots.faa
bin.1.permissive_totalprots_vs_gtdbplus_protdb.dmnd.blast.tsv
bin.1.permissive_tRNAs.fasta.gz
bindata_progress.pickle
bindata_trna_progress.json.gz
fullcontiginfos_beforecleanup.tsv
nucblasts.json.gz
prok_cog.hmm.domtblout
prok_pfam.hmm.domtblout
prok_tigr.hmm.domtblout
protblasts.json.gz

7-mdmcleaner/SRR8555091_sub/bin.2.strict:
arch_pfam.hmm.domtblout
arch_tigr.hmm.domtblout
bact_pfam.hmm.domtblout
bact_tigr.hmm.domtblout
bin.2.strict_filtered_kept_contigs.fasta.gz
bin.2.strict_filtered_need_evaluation_high.fasta.gz
bin.2.strict_filtered_need_evaluation_low.fasta.gz
bin.2.strict_filtered_removed.fasta.gz
bin.2.strict_kronainput.tsv
bin.2.strict_totalprots.faa
bin.2.strict_totalprots_vs_gtdbplus_protdb.dmnd.blast.tsv
bin.2.strict_tRNAs.fasta.gz
bindata_progress.pickle
bindata_trna_progress.json.gz
fullcontiginfos_beforecleanup.tsv
nucblasts.json.gz
prok_cog.hmm.domtblout
prok_pfam.hmm.domtblout
prok_tigr.hmm.domtblout
protblasts.json.gz
```

The `fullcontiginfos_beforecleanup.tsv` file contains the information reagrding the quality of each contig and the indications of the software program. Here is where the category for each contigs (keep, remove...) is presented next to the index value.

Moreover, the `<bin_name>_filtered_kept_contigs.fasta.gz` is were all the good quality contigs are presented. This is the file that will be used for further analysis.

### Running GUNC

GUNC is installed as a singularity container. Therefore, it is requried to call the complete path to the containers folder.

```bash
(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/gunc-1.0.6.sif gunc -h
usage: gunc [-h] [-v]  ...

Tool for detection of chimerism and contamination in prokaryotic genomes.

options:
  -h, --help     show this help message and exit
  -v, --version  Print version number and exit.

GUNC subcommands:

    run          Run chimerism detection.
    download_db  Download GUNC db.
    merge_checkm
                 Merge GUNC and CheckM outputs.
    plot         Create interactive visualisation.
    summarise    Reevaluate results using a different contamination_portion cutoff.
```

The software holds different commands including one for downloading the database. As priorly mentioned, the database by default is progenomes, which is already downloaded in the computer. Check the options and basicun command to perform the MAGs analysis.

```bash
(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/gunc-1.0.6.sif gunc run -h
usage: gunc run [-h] [-r] (-i  | -f  | -d ) [-e] [-g] [-t] [-o] [--temp_dir] [--sensitive]
                [--detailed_output] [--contig_taxonomy_output] [--use_species_level] [--min_mapped_genes]
                [-v]

options:
  -h, --help            show this help message and exit
  -r , --db_file        DiamondDB reference file. Default: GUNC_DB envvar
  -i , --input_fasta    Input file in FASTA format.
  -f , --input_file     File with paths to FASTA format files.
  -d , --input_dir      Input dir with files in FASTA format.
  -e , --file_suffix    Suffix of files in input_dir. Default: .fa
  -g, --gene_calls      Input files are FASTA faa format. Default: False
  -t , --threads        number of CPU threads. Default: 4
  -o , --out_dir        Output dir. Default: cwd
  --temp_dir            Directory to store temp files. Default: cwd
  --sensitive           Run with high sensitivity. Default: False
  --detailed_output     Output scores for every taxlevel. Default: False
  --contig_taxonomy_output
                        Output assignments for each contig. Default: False
  --use_species_level   Allow species level to be picked as maxCSS. Default: False
  --min_mapped_genes    Dont calculate GUNC score if number of mapped genes is below this value. Default:
                        11
  -v, --verbose         Verbose output for debugging
```

Notice the default for GUNC is the suffix `.fa`. Therefore, `-e` command is required to indicate the mdmcleaner results files.

Moreover, the database for the gunc software needs to be indicated in the command through the `-r` flag. This db is presented in the path `/home/public/DB/gunc_db_progenomes2.1.dmnd`.

Additionally, the results from mdmcleaner are presented in different subfolders for each bins. These need to be copied into a single subfolder to be used as the input for the gunc command (`-d` flag). In this subfolder, the names of each bins require to have the the prefix with the sample code. For this, the following code must be added after before running.

```vim
#Adding sample name prefix
        for file in 7-mdmcleaner/${sample}_sub/bins/*.fasta.gz; do
                mv "$file" "$(dirname "$file")/${sample}_$(basename "$file")"
        done
```

Code for a simple analysis of the bins that were kept after the mdmcleaner for each bin. Save the results in a `8-gunc` folder and set the threads to the maximum available per node.

By the end of the code, create a new folder in the `7-mdmcleaner` directory called `bins` and copy all the `*fasta.gz` files that were analyzed. For this, include something similar to this after running mdmcleaner.

```vim
#Creating a final folder with all the bins
mkdir -p 7-mdmcleaner/bins
cp 7-mdmcleaner/${sample}_sub/bins/*.fasta.gz 7-mdmcleaner/bins
```

This code is more complex as it includes the addition of `*.fasta.gz` files to a new folder and their rename. Hence, remember the template for this `.slurm` file is in the end of this document.

> If the job takes a while running, the results in the common repository `/home/public/met-workshop`

The output from GUNC should look something like this:

```bash
(base) [dorian.rojas@accessnode test]$ ls 8-gunc/*/
8-gunc/SRR8555091_sub/:
diamond_output  gene_calls  GUNC.progenomes_2.1.maxCSS_level.tsv

8-gunc/SRR9988205_sub/:
diamond_output  gene_calls  GUNC.progenomes_2.1.maxCSS_level.tsv
```

The `GUNC.progenomes_2.1.maxCSS_level.tsv` contains the table with the GUNC analysis results. Here, the last columnn `pass.GUNC` indicates whether the bins passed (`True`) or failed (`False`) the GUNC analysis.

```vim
genome  n_genes_called  n_genes_mapped  n_contigs       taxonomic_level proportion_genes_retained_in_major_clades   genes_retained_index    clade_separation_score  contamination_portion   n_effective_surplus_clades mean_hit_identity        reference_representation_score  pass.GUNC
SRR9988205_bin.1.permissive_filtered_kept_contigs.fasta 3983    3951    781     genus   0.93    0.93    0.060.14    0.31    0.99    0.92    True
SRR9988205_bin.2.permissive_filtered_kept_contigs.fasta 3737    3591    544     kingdom 1.0     0.96    0.00.0      0.0     0.96    0.92    True
```

Basicallly, those bins with the category `True` are the one further selected for following analysis.

### Running CheckM2

CheckM2 is installed as a container and has multiple commands.

```bash
(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/checkm2-1.0.2.sif checkm2 -h
                  ____ _               _    __  __ ____
                 / ___| |__   ___  ___| | _|  \/  |___ \
                | |   | '_ \ / _ \/ __| |/ / |\/| | __) |
                | |___| | | |  __/ (__|   <| |  | |/ __/
                 \____|_| |_|\___|\___|_|\_\_|  |_|_____|

                ...::: CheckM2 v1.0.2 :::...

  General usage:
    predict         -> Predict the completeness and contamination of genome bins in a folder. Example usage:


        checkm2 predict --threads 30 --input <folder_with_bins> --output-directory <output_folder>
    testrun         -> Runs Checkm2 on internal test genomes to ensure it runs without errors. Example usage:

         checkm2 testrun --threads 10
    database        -> Download and set up required CheckM2 DIAMOND database for annotation

  Use checkm2 <command> -h for command-specific help.
```

Check the `predict` command to evaluate completion and contamination of bins.

```bash
(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/checkm2-1.0.2.sif checkm2 predict -h
usage: checkm2 predict [-h] [--debug] [--version] [--quiet] [--lowmem] --input INPUT [INPUT ...]
                       --output-directory OUTPUT_DIRECTORY [--general] [--specific] [--allmodels]
                       [--genes] [-x EXTENSION] [--tmpdir TMPDIR] [--force] [--resume]
                       [--threads num_threads] [--stdout] [--remove_intermediates] [--ttable ttable]
                       [--database_path DATABASE_PATH] [--dbg_cos] [--dbg_vectors]

Predict the completeness and contamination of genome bins in a folder. Example usage:

        checkm2 predict --threads 30 --input <folder_with_bins> --output-directory <output_folder>

optional arguments:
  -h, --help            show this help message and exit
  --debug               output debug information
  --version             output version information and quit
  --quiet               only output errors
  --lowmem              Low memory mode. Reduces DIAMOND blocksize to significantly reduce RAM usage at the expense of longer runtime

required arguments:
  --input INPUT [INPUT ...], -i INPUT [INPUT ...]
                        Path to folder containing MAGs or list of MAGS to be analyzed
  --output-directory OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY, -o OUTPUT_DIRECTORY
                        Path output to folder

additional arguments:
  --general             Force the use of the general quality prediction model (gradient boost)
  --specific            Force the use of the specific quality prediction model (neural network)
  --allmodels           Output quality prediction for both models for each genome.
  --genes               Treat input files as protein files. [Default: False]
  -x EXTENSION, --extension EXTENSION
                        Extension of input files. [Default: .fna]
  --tmpdir TMPDIR       specify an alternative directory for temporary files
  --force               overwrite output directory [default: do not overwrite]
  --resume              Reuse Prodigal and DIAMOND results found in output directory [default: not set]
  --threads num_threads, -t num_threads
                        number of CPUS to use [default: 1]
  --stdout              Print results to stdout [default: write to file]
  --remove_intermediates
                        Remove all intermediate files (protein files, diamond output) [default: dont]
  --ttable ttable       Provide a specific progidal translation table for bins [default: automatically determine either 11 or 4]
  --database_path DATABASE_PATH
                        Provide a location for the CheckM2 database for a given predict run [default: use either internal path set via <checkm2 database> or CHECKM2DB environmental variable]
  --dbg_cos             DEBUG: write cosine similarity values to file [default: don't]
  --dbg_vectors         DEBUG: dump pickled feature vectors to file [default: don't]
```

Do not get overwhelmed by the large help message, the command is very simple, check the usage example first. In this case the default extension for the files is `.fna`. Hence, the `.gz` of the mdmcleaner results should be specificed in the CheckM2 command with the `-x` flags.

The database for CheckM2 must be indicated in the predict command. The path is `/home/dorian.rojas/DB/CheckM2_database/uniref100.KO.1.dmnd`.

Code for a prediction analysis of the same bins used for the GUNC assessment and set the threats to the maximum to ensure a faster running. Save the results in a folder named `9-checkm2`.

> If the job takes a while running, the results in the common repository `/home/public/met-workshop`

Evaluate the results folder.

```bash
(base) [dorian.rojas@accessnode test]$ ls 9-checkm2/*/
9-checkm2/SRR8555091/:
checkm2.log  diamond_output  protein_files  quality_report.tsv

9-checkm2/SRR9988205/:
checkm2.log  diamond_output  protein_files  quality_report.tsv
```

The tab separated table `quality_report.tsv` contains the quality evaluation of the bins.

```vim
Name    Completeness    Contamination   Completeness_Model_Used Translation_Table_Used  Coding_Density  Contig_N50  Average_Gene_Length     Genome_Size     GC_Content      Total_Coding_Sequences  Total_Contigs   Max_Contig_Length   Additional_Notes
SRR9988205_bin.1.permissive_filtered_kept_contigs.fasta 82.96   0.45    Neural Network (Specific Model) 11 0.896    5827    284.7303144654088       3781789 0.51    3975    781     37930   None
SRR9988205_bin.2.permissive_filtered_kept_contigs.fasta 95.32   0.27    Neural Network (Specific Model) 11 0.904    12541   360.1825247922809       4456678 0.44    3731    545     73642   None
```

## Data analysis

To further complete the analyss of this data, it is recommended to combine both the GUNC and CheckM2 results. This will allow an easier identification and analysis of those bins that passed the GUNC and have the quality parameters of medium-quality bins.

For this, a python code was created under the name `gtdb-input.py`. The code is presented below. The directory paths and output file path need to be modified prior its use. Copy the python file from the common public folder or copy and paste the code below in a text file in the console.

Prior to running, it is recommended to create a new directory named `10-final_mags`. Here is where all the output table should be store.

> This code is a modification from the template provided by the M.Sc Maria Alejandra Soto, a secondee from the project

```python
import pandas as pd
import os
import re

# Indicating directories for GUNC
PD = '/home/dorian.rojas/test/8-gunc'
OF = '/home/dorian.rojas/test/8-gunc/combined_gunc.tsv'

# List of dataframes
dataframes = []

# Explore PD
for root, dirs, files in os.walk(PD):
    for file in files:
        if file == 'GUNC.progenomes_2.1.maxCSS_level.tsv':
            file_path = os.path.join(root, file)
            # Add to list
            df = pd.read_csv(file_path, sep='\t')
            dataframes.append(df)

# Concat
combined_df = pd.concat(dataframes, ignore_index=True)

# Save
combined_df.to_csv(OF, sep='\t', index=False)

print(f'Combined gunc table saved as: {OF}')


# Indicating directories for CheckM2
PDir = '/home/dorian.rojas/test/9-checkm2'
Ofile = '/home/dorian.rojas/test/9-checkm2/combined_checkm2.tsv'

# Dataframes
Dframes = []

# Explore Pdir
for root, dirs, files in os.walk(PDir):
    for file in files:
        if file == 'quality_report.tsv':
            file_path = os.path.join(root, file)
            # Add to list
            df = pd.read_csv(file_path, sep = '\t')
            Dframes.append(df)

# Concat
combined_df = pd.concat(Dframes, ignore_index=True)

# Save
combined_df.to_csv(Ofile, sep = '\t', index = False)

print(f'Combined checkm table saved as: {Ofile}')


# Add GUNC to CheckM
pd.set_option('display.max_colwidth', None)

# Read CheckM2
checkm2_all = pd.read_csv('/home/dorian.rojas/test/9-checkm2/combined_checkm2.tsv', sep='\t')
# Rename 'Name' column to match GUNC table
checkm2_all = checkm2_all.rename(columns={'Name':'genome'})

# Read GUNC, 
gunc_all = pd.read_csv('/home/dorian.rojas/test/8-gunc/combined_gunc.tsv', sep='\t')
# Removing '.' trailing
gunc_all['genome'] = gunc_all['genome'].str.replace(r'_filtered_kept_contigs\.fasta\.$', '_filtered_kept_contigs.fasta', regex=True) 


# Add pass.GUNC column. FO = Final Output
FO = '/home/dorian.rojas/test/10-final-MAGs/combined_gunc_and_checkm2.tsv'

combined_gunc_and_checkm2 = checkm2_all.merge(gunc_all[['genome', 'pass.GUNC']], on='genome', how='left')
combined_gunc_and_checkm2.to_csv(FO, sep='\t', index=False)

print(f'Final output table was saved at: {FO}')


# Remove rows of bins that failed GUNC and do not have the medium-quality
pass_output = '/home/dorian.rojas/test/10-final-MAGs/final_mags.tsv'

final_mags = combined_gunc_and_checkm2.drop(combined_gunc_and_checkm2[combined_gunc_and_checkm2['pass.GUNC'] == False].index)

final_mags = final_mags[
    (final_mags['Completeness'] >= 50) &
    (final_mags['Contamination'] <= 5)
]

final_mags.to_csv(pass_output, index = False, sep='\t')

print(f'Final medium-quality MAGs that passed GUNC analysis table was saved at: {pass_output}')


# Create gtdbtk input
mags_set = '/home/dorian.rojas/test/7-mdmcleaner/bins'
kept_mags = pass_output
path_and_name = '/home/dorian.rojas/test/10-final-MAGs/gtdb-input.txt'

# Read table
df = pd.read_csv(kept_mags, sep='\t')
# Rename
mags_ids = set(df['genome'].str.replace(r'_filtered_kept_contigs\.fasta$', '', regex=True))
# Exctrac mag ID
pattern = re.compile(r'(.+?)_filtered_kept_contigs\.fasta\.gz$')

# Write ouput file
with open(path_and_name, 'w') as out_file:
    # Tranverse files in folder
    for file_name in os.listdir(mags_set):
        # Check correct extension
        if file_name.endswith('_filtered_kept_contigs.fasta.gz'):
            # Use pattern
            match = pattern.match(file_name)
            if match:
                mag_id = match.group(1)
                # Check if in the list 
                if mag_id in mags_ids:
                    # Get path
                    file_path = os.path.join(mags_set, file_name)
                    # Write ID
                    out_file.write(f'{file_path}\t{mag_id}\n')

print(f'MAGs path and name file saved as: {path_and_name}')
```

To run the python code use the command `python <file_name>`. If there is an error regarding the installation of pandas run `pip install pandas` in the console and try rerunning the python code.

The results from this python code are three files. These correspond to the concatenation of the GUNC and CheckM2 results (`combined_gunc_and_checkm2.tsv`), the collection of those bins that passed the GUNC analysis and have medium-quality statistics (`final_mags.tsv`), and a table that will be used for further analysis (`gtdb-input.txt`).

```vim
genome  Completeness    Contamination   Completeness_Model_Used Translation_Table_Used  Coding_Density  Contig_N50  Average_Gene_Length     Genome_Size     GC_Content      Total_Coding_Sequences  Total_Contigs   Max_Contig_Length   Additional_Notes        pass.GUNC
SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta 86.0    1.87    Neural Network (Specific Model) 11 0.915    5037    327.78  2146512 0.42    2000    487     22248           True
SRR8555091_bin.2.strict_filtered_kept_contigs.fasta     78.16   2.29    Neural Network (Specific Model) 11 0.934    7463    359.9509202453988       2633626 0.49    2282    443     34713           True
SRR9988205_bin.1.permissive_filtered_kept_contigs.fasta 82.96   0.45    Neural Network (Specific Model) 11 0.896    5827    284.7303144654088       3781789 0.51    3975    781     37930           True
SRR9988205_bin.2.permissive_filtered_kept_contigs.fasta 95.32   0.27    Neural Network (Specific Model) 11 0.904    12541   360.1825247922809       4456678 0.44    3731    545     73642           True
```

```vim
genome  Completeness    Contamination   Completeness_Model_Used Translation_Table_Used  Coding_Density  Contig_N50  Average_Gene_Length     Genome_Size     GC_Content      Total_Coding_Sequences  Total_Contigs   Max_Contig_Length   Additional_Notes        pass.GUNC
SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta 86.0    1.87    Neural Network (Specific Model) 11 0.915    5037    327.78  2146512 0.42    2000    487     22248           True
SRR8555091_bin.2.strict_filtered_kept_contigs.fasta     78.16   2.29    Neural Network (Specific Model) 11 0.934    7463    359.9509202453988       2633626 0.49    2282    443     34713           True
SRR9988205_bin.1.permissive_filtered_kept_contigs.fasta 82.96   0.45    Neural Network (Specific Model) 11 0.896    5827    284.7303144654088       3781789 0.51    3975    781     37930           True
SRR9988205_bin.2.permissive_filtered_kept_contigs.fasta 95.32   0.27    Neural Network (Specific Model) 11 0.904    12541   360.1825247922809       4456678 0.44    3731    545     73642           True
```

It is noticeable that the combined and final table are the same. This is expected as all the resulting bins passed the required parameters. Therefore, a similar table is the output.

Finally, the last output file from this code is what's called a "manisfest" or "batch file" for some software. These correspond to a `.tsv` table with the full path of the `.fasta`  files and their respective name.

```vim
/home/dorian.rojas/test/7-mdmcleaner/bins/SRR8555091_bin.1.permissive_filtered_kept_contigs.fasta.gz    SRR8555091_bin.1.permissive
/home/dorian.rojas/test/7-mdmcleaner/bins/SRR8555091_bin.2.strict_filtered_kept_contigs.fasta.gz        SRR8555091_bin.2.strict
/home/dorian.rojas/test/7-mdmcleaner/bins/SRR9988205_bin.1.permissive_filtered_kept_contigs.fasta.gz    SRR9988205_bin.1.permissive
/home/dorian.rojas/test/7-mdmcleaner/bins/SRR9988205_bin.2.permissive_filtered_kept_contigs.fasta.gz    SRR9988205_bin.2.permissive
```

Each software that allows the utlization of a batch file have a specific style for the document. The output from this python file is specific for the GTDB-tk tools used for taxonomic annotations. This is the reason behind the name `gtdb-input.txt`.

## Code examples

**`mdmcleaner.slurm` file:**

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="mdmcleaner"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL

cd /home/dorian.rojas/test
WD=/home/dorian.rojas/test

for sample in $@; do

mkdir -p 7-mdmcleaner/${sample}_sub
cd 7-mdmcleaner/${sample}_sub

mdmcleaner set_configs --db_basedir /home/public/DB
mdmcleaner clean -i $WD/6-reassemble/${sample}_sub/reassembled_bins/*.fa \
        -o $WD/7-mdmcleaner/${sample}_sub/ -t 64

cd ../../

done

date
time
```

**`gunc.slurm` file:**

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="gunc"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL

cd /home/dorian.rojas/test

CTN_PATH=/opt/ohpc/pub/containers/BIO/
$CTN_PATH/gunc-1.0.6.sif gunc -v

for sample in $@; do

mkdir -p 8-gunc/${sample}_sub

#Moving the mdmcleaner files into a singlular subfolder
mkdir -p 7-mdmcleaner/${sample}_sub/bins
cp 7-mdmcleaner/${sample}_sub/bin.*/*_kept_contigs.fasta.gz 7-mdmcleaner/${sample}_sub/bins

        #Adding sample name prefix
        for file in 7-mdmcleaner/${sample}_sub/bins/*.fasta.gz; do
                mv "$file" "$(dirname "$file")/${sample}_$(basename "$file")"
        done

#Running gunc
$CTN_PATH/gunc-1.0.6.sif gunc run -d 7-mdmcleaner/${sample}_sub/bins \
       -o 8-gunc/${sample}_sub/ -e .gz -t 64 -r /home/dorian.rojas/DB/gunc_db_progenomes2.1.dmnd


#Creating a final folder with all the bins
mkdir -p 7-mdmcleaner/bins
cp 7-mdmcleaner/${sample}_sub/bins/*.fasta.gz 7-mdmcleaner/bins

done

date
time
```

**`checkm2.slurm` file:**

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="checkm2"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL

cd /home/dorian.rojas/test

CTN_PATH=/opt/ohpc/pub/containers/BIO/
DBPATH=/home/dorian.rojas/DB/CheckM2_database/uniref100.KO.1.dmnd

for sample in $@; do

mkdir -p 9-checkm2/${sample}

#export CHECKM2DB="/home/public/met-workshop/uniref100.KO.1.dmnd"

$CTN_PATH/checkm2-1.0.2.sif checkm2 predict --threads $SLURM_NTASKS \
         --input 7-mdmcleaner/${sample}_sub/bins --output-directory 9-checkm2/$sample/ \
        -x .gz --database_path $DBPATH

done

date
time
```

**`gtdb-input.py` file:**

```python
import pandas as pd
import os
import re

# Indicating directories for GUNC
PD = '/home/dorian.rojas/test/8-gunc'
OF = '/home/dorian.rojas/test/8-gunc/combined_gunc.tsv'

# List of dataframes
dataframes = []

# Explore PD
for root, dirs, files in os.walk(PD):
    for file in files:
        if file == 'GUNC.progenomes_2.1.maxCSS_level.tsv':
            file_path = os.path.join(root, file)
            # Add to list
            df = pd.read_csv(file_path, sep='\t')
            dataframes.append(df)

# Concat
combined_df = pd.concat(dataframes, ignore_index=True)

# Save
combined_df.to_csv(OF, sep='\t', index=False)

print(f'Combined gunc table saved as: {OF}')


# Indicating directories for CheckM2
PDir = '/home/dorian.rojas/test/9-checkm2'
Ofile = '/home/dorian.rojas/test/9-checkm2/combined_checkm2.tsv'

# Dataframes
Dframes = []

# Explore Pdir
for root, dirs, files in os.walk(PDir):
    for file in files:
        if file == 'quality_report.tsv':
            file_path = os.path.join(root, file)
            # Add to list
            df = pd.read_csv(file_path, sep = '\t')
            Dframes.append(df)

# Concat
combined_df = pd.concat(Dframes, ignore_index=True)

# Save
combined_df.to_csv(Ofile, sep = '\t', index = False)

print(f'Combined checkm table saved as: {Ofile}')


# Add GUNC to CheckM
pd.set_option('display.max_colwidth', None)

# Read CheckM2
checkm2_all = pd.read_csv('/home/dorian.rojas/test/9-checkm2/combined_checkm2.tsv', sep='\t')
# Rename 'Name' column to match GUNC table
checkm2_all = checkm2_all.rename(columns={'Name':'genome'})

# Read GUNC,
gunc_all = pd.read_csv('/home/dorian.rojas/test/8-gunc/combined_gunc.tsv', sep='\t')
# Removing '.' trailing
gunc_all['genome'] = gunc_all['genome'].str.replace(r'_filtered_kept_contigs\.fasta\.$', '_filtered_kept_contigs.fasta', regex=True)


# Add pass.GUNC column. FO = Final Output
FO = '/home/dorian.rojas/test/10-final-MAGs/combined_gunc_and_checkm2.tsv'

combined_gunc_and_checkm2 = checkm2_all.merge(gunc_all[['genome', 'pass.GUNC']], on='genome', how='left')
combined_gunc_and_checkm2.to_csv(FO, sep='\t', index=False)

print(f'Final output table was saved at: {FO}')


# Remove rows of bins that failed GUNC and do not have the medium-quality
pass_output = '/home/dorian.rojas/test/10-final-MAGs/final_mags.tsv'

final_mags = combined_gunc_and_checkm2.drop(combined_gunc_and_checkm2[combined_gunc_and_checkm2['pass.GUNC'] == False].index)

final_mags = final_mags[
    (final_mags['Completeness'] >= 50) &
    (final_mags['Contamination'] <= 5)
]

final_mags.to_csv(pass_output, index = False, sep='\t')

print(f'Final medium-quality MAGs that passed GUNC analysis table was saved at: {pass_output}')


# Create gtdbtk input
mags_set = '/home/dorian.rojas/test/7-mdmcleaner/bins'
kept_mags = pass_output
path_and_name = '/home/dorian.rojas/test/10-final-MAGs/gtdb-input.txt'

# Read table
df = pd.read_csv(kept_mags, sep='\t')
# Rename
mags_ids = set(df['genome'].str.replace(r'_filtered_kept_contigs\.fasta$', '', regex=True))
# Exctrac mag ID
pattern = re.compile(r'(.+?)_filtered_kept_contigs\.fasta\.gz$')

# Write ouput file
with open(path_and_name, 'w') as out_file:
    # Tranverse files in folder
    for file_name in os.listdir(mags_set):
        # Check correct extension
        if file_name.endswith('_filtered_kept_contigs.fasta.gz'):
            # Use pattern
            match = pattern.match(file_name)
            if match:
                mag_id = match.group(1)
                # Check if in the list
                if mag_id in mags_ids:
                    # Get path
                    file_path = os.path.join(mags_set, file_name)
                    # Write ID
                    out_file.write(f'{file_path}\t{mag_id}\n')
```
