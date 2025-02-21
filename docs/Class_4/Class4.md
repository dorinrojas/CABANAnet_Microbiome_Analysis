# Class 4: Functional and taxonomical annotation of MAGs

- - - -

## Taxonomical annotatiion of MAGs (GTDBtk)

Until this section of the workshop, the output are solemnly fasta sequences. Clean, filtered, and perfect fasta sequences that are collected in bins. Each of these bins represent, in theory, one species.

The following steps towards the interpretation of this information is the annotation. Particularly, the taxonomical annotation to have an overall idea of which bacteria is presented in the analyzed microbiome. This is strongly relevant for the analysis as knowing the microbial profile is one of the ultimate goals of this workshop.

To perform this taxonomical annotation, multiple tools are available. For instance, [PhyloPhlAn](https://github.com/biobakery/phylophlan) and [MiGA](https://github.com/bio-miga/miga) (Microbial Genome Atlas) are based on the NCBI taxonomy database and have been widely used. More compherensive tools, such as the Genome Taxonomy Database toolkit (GTDB-tk) employ the GTDB for identification.

Commonly, the taxonomic annotation is performed based on the aligments of the marker genes, such as the 16S rRNA gene, or aligment of database to the input bins. However, this is labour intensive and computational demanding. Therefore, other approach like the analysis based on statistics like percentage identity or Average Nucleotide Identity (ANI) have been implemented.

The GTDB-tk performs the taxonomic classification by placing bins into protein referecnes tree that are domain-specific and uses the criteria from Relative Evolutionary Divergence and ANI statistics.

For this, the toolkit uses the reference trees, alignments, and taxonomies available in the GTDB website. To make the initial placement, the toolkit identifies the genes using Prodigal and selects marker genes using HMMER. The input genomes are assigned to those domain-specific trees based on the highest proportion of identifies genes. The genes are aligned and concatenated into a single files that is used for placement using pplacer [(Chaumeil et al. 2020)](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182?login=false).

Then, through the placements, RED and ANI values a taxonomy is assigned. RED values are used to resolve when rank assignments are ambiguous. Meanwhile, species are determined usinf the ANI values calculated based on the FastANI software. The criteria for genus is >65% and the species indicator is 95%. If these parameters are not met, the bins is classified as a novel genus or species, respectively [(Chaumeil et al. 2020)](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182?login=false).

Overall, GTDB-tk is currently widely used for taxonomic annotation as the GTDB has gained popularity in the metagenomics research community. The workshop employs the second version of the software (GTDB-tk v2). This is very similar to the previous version. The main difference relays on the division of reference tree into class-level trees rather than domain-specific. This allows a more memory efficient approach for the taxonomy assignation process, as the GTDB has significantly grown in the last years [(Chaumeil et al. 2022)](https://academic.oup.com/bioinformatics/article/38/23/5315/6758240).

### Running gtdbtk

GTDB-tk is installed through a singularity container. The software presents multiple options as it represents a toolkit.

```bash
(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/gtdbtk-2.4.0.sif gtdbtk -h

              ...::: GTDB-Tk v2.4.0 :::...

  Workflows:
    classify_wf -> Classify genomes by placement in GTDB reference tree
                     (ani_screening -> identify -> align -> classify)
    de_novo_wf  -> Infer de novo tree and decorate with GTDB taxonomy
                     (identify -> align -> infer -> root -> decorate)

  Methods:
    identify -> Identify marker genes in genome
    align    -> Create multiple sequence alignment
    classify -> Determine taxonomic classification of genomes
    infer    -> Infer tree from multiple sequence alignment
    root     -> Root tree using an outgroup
    decorate -> Decorate tree with GTDB taxonomy

  Tools:
    infer_ranks        -> Establish taxonomic ranks of internal nodes using RED
    ani_rep            -> Calculates ANI to GTDB representative genomes
    trim_msa           -> Trim an untrimmed MSA file based on a mask
    export_msa         -> Export the untrimmed archaeal or bacterial MSA file
    remove_labels      -> Remove labels (bootstrap values, node labels) from an Newick tree
    convert_to_itol    -> Convert a GTDB-Tk Newick tree to an iTOL tree
    convert_to_species -> Convert GTDB genome IDs to GTDB species names


  Testing:
    test          -> Validate the classify_wf pipeline with 3 archaeal genomes
    check_install -> Verify third party programs and GTDB reference package

  Use: gtdbtk <command> -h for command specific help
```

Similar to CheckM, this toolkit have integrated complete workflows that involves several of the methods of the tools. In this regard, the `classify_wf` options is the required for the taxonomical annotation of the bins.

This module includes the `identify`, `align`, and `classify` modules of the complete toolkit. The `identify` module predicts the marker genes using Prodigal based on the 11 translation table. Moreover, the `align` module creates the multiple sequence aligment for the placement of the bins into the class-level trees. The marker genes used for this alignment are the AR53/BAC120 genes. Finaly, the `classify` module determines the taxonomic classification of the metagenomics bins.

```bash
(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/gtdbtk-2.4.0.sif gtdbtk classify_wf -h
usage: gtdbtk classify_wf (--genome_dir GENOME_DIR | --batchfile BATCHFILE) --out_dir OUT_DIR
                          (--skip_ani_screen | --mash_db MASH_DB) [--no_mash] [--mash_k MASH_K]
                          [--mash_s MASH_S] [--mash_v MASH_V] [--mash_max_distance MASH_MAX_DISTANCE] [-f]
                          [-x EXTENSION] [--min_perc_aa MIN_PERC_AA] [--prefix PREFIX] [--genes]
                          [--cpus CPUS] [--pplacer_cpus PPLACER_CPUS] [--force]
                          [--scratch_dir SCRATCH_DIR] [--write_single_copy_genes] [--keep_intermediates]
                          [--min_af MIN_AF] [--tmpdir TMPDIR] [--debug] [-h]

mutually exclusive required arguments:
  --genome_dir GENOME_DIR
                        directory containing genome files in FASTA format
  --batchfile BATCHFILE
                        path to file describing genomes - tab separated in 2 or 3 columns (FASTA file,
                        genome ID, translation table [optional])

required named arguments:
  --out_dir OUT_DIR     directory to output files

mutually exclusive required arguments:
  --skip_ani_screen     Skip the ani_screening step to classify genomes using mash and skani. (default:
                        False)
  --mash_db MASH_DB     path to save/read (if exists) the Mash reference sketch database (.msh)

optional Mash arguments:
  --no_mash             skip pre-filtering of genomes using Mash (default: False)
  --mash_k MASH_K       k-mer size [1-32] (default: 16)
  --mash_s MASH_S       maximum number of non-redundant hashes (default: 5000)
  --mash_v MASH_V       maximum p-value to keep [0-1] (default: 1.0)
  --mash_max_distance MASH_MAX_DISTANCE
                        Maximum Mash distance to select a potential GTDB genome as representative of a
                        user genome. (default: 0.15)

optional arguments:
  -f, --full_tree       use the unsplit bacterial tree for the classify step; this is the original GTDB-Tk
                        approach (version < 2) and requires more than 320 GB of RAM to load the reference
                        tree (default: False)
  -x, --extension EXTENSION
                        extension of files to process, gz = gzipped (default: fna)
  --min_perc_aa MIN_PERC_AA
                        exclude genomes that do not have at least this percentage of AA in the MSA
                        (inclusive bound) (default: 10)
  --prefix PREFIX       prefix for all output files (default: gtdbtk)
  --genes               indicates input files contain predicted proteins as amino acids (skip gene
                        calling).Warning: This flag will skip the ANI comparison steps (ani_screen and
                        classification). (default: False)
  --cpus CPUS           number of CPUs to use (default: 1)
  --pplacer_cpus PPLACER_CPUS
                        number of CPUs to use during pplacer placement
  --force               continue processing if an error occurs on a single genome (default: False)
  --scratch_dir SCRATCH_DIR
                        reduce pplacer memory usage by writing to disk (slower).
  --write_single_copy_genes
                        output unaligned single-copy marker genes (default: False)
  --keep_intermediates  keep intermediate files in the final directory (default: False)
  --min_af MIN_AF       minimum alignment fraction to assign genome to a species cluster (default: 0.5)
  --tmpdir TMPDIR       specify alternative directory for temporary files (default: /tmp)
  --debug               create intermediate files for debugging purposes (default: False)
  -h, --help            show help message
```

Code for the `classify_wf` workflow for the indicating the batch file created with the python script as the input and a novel `11-gtdbtk` as the outpur directory. Notice the default extension of the software is `.fna`. Therefore, set the extension to `.gz` using the `-x` flag. In this case, threats are specified using the `--cpus` option. Also, due to the time consuming nature of the process, add the `-skip_ani_screen` function.

The database for the gtdbtk must be indicated by adding `export GTDBTK_DATA_PATH="/home/public/DB/release220/"` prior the tool command. `export` is a base command that sets an environmental variable. Most variables are different depending on the running software; therefore, this information has to be retrieved from the source code or the github repository.

This analysis needs to be run only once in the batch file. As the input is not divided by the samples and the batch file includes samples from all the bins. This allows the `.slurm` to be slightly different. For instance, the `for` loop can be removed from the command the output is set to the general directory, avoiding the creation of sample-specific subfolders.

Finally, the `.slurm` file is indicated using the `sbatch` command in the console. The `batch.sh` is not required.

> If the job takes a while running, the results in the common repository `/home/public/met-workshop`

Evaluate the results folder.

```bash
(base) [dorian.rojas@accessnode test]$ ls 11-gtdbtk/
align  classify  gtdbtk.bac120.summary.tsv  gtdbtk.json  gtdbtk.log  gtdbtk.warnings.log  identify
```

There are individuals folders with the results for each GTDB-tk module. The `gtdbtk.bac120.summary.tsv` contains the information regarding the taxonomic classification for the respective marker gene. This tab-separared table is divides in several columns.

```vim
user_genome     classification  closest_genome_reference        closest_genome_reference_radius closest_genome_taxonomy     closest_genome_ani      closest_genome_af       closest_placement_reference     closest_placement_radius    closest_placement_taxonomy      closest_placement_ani   closest_placement_af    pplacer_taxonomy    classification_method   note    other_related_references(genome_id,species_name,radius,ANI,AF)  msa_percent translation_table       red_value       warnings
SRR8555091_bin.1.permissive     d__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Treponematales;f__Treponemataceae;g__Treponema_D;s__Treponema_D sp900541945 GCA_022010055.1 95.0    d__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Treponematales;f__Treponemataceae;g__Treponema_D;s__Treponema_D sp900541945 98.12   0.752   GCA_022010055.1     95.0    d__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Treponematales;f__Treponemataceae;g__Treponema_D;s__Treponema_D sp900541945 98.12   0.752   d__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Treponematales;f__Treponemataceae;g__Treponema_D;s__        taxonomic classification defined by topology and ANItopological placement and ANI have congruent species assignments        GCA_900317625.1, s__Treponema_D sp900317625, 95.0, 85.39, 0.268; GCA_002478955.1, s__Treponema_D sp002478955, 95.0, 85.14, 0.279; GCA_934718065.1, s__Treponema_D sp934718065, 95.0, 84.66, 0.211       79.66   11      N/A     N/A
SRR8555091_bin.2.strict d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__UBA932;g__Cryptobacteroides;s__Cryptobacteroides sp000434935 GCA_000434935.1 95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__UBA932;g__Cryptobacteroides;s__Cryptobacteroides sp000434935 96.83   0.773   GCA_000434935.1     95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__UBA932;g__Cryptobacteroides;s__Cryptobacteroides sp000434935 96.83   0.773   d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__UBA932;g__Cryptobacteroides;s__      taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments    N/A     72.67   11      N/A     N/A
SRR9988205_bin.1.permissive     d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli   GCF_003697165.2 95.0    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli   98.39       0.921   GCF_003697165.2 95.0    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli   98.39   0.921   d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__   taxonomic classification defined by topology and ANI        topological placement and ANI have congruent species assignments   GCF_000194175.1, s__Escherichia coli_F, 95.0, 95.66, 0.841; GCF_002965065.1, s__Escherichia sp002965065, 95.0, 94.29, 0.69; GCF_004211955.1, s__Escherichia sp004211955, 95.0, 93.15, 0.728; GCF_005843885.1, s__Escherichia sp005843885, 95.0, 92.99, 0.746; GCF_029876145.1, s__Escherichia ruysiae, 95.0, 92.88, 0.744; GCF_000026225.1, s__Escherichia fergusonii, 95.0, 92.75, 0.565; GCF_011881725.1, s__Escherichia coli_E, 95.0, 92.4, 0.713; GCF_014836715.1, s__Escherichia whittamii, 95.0, 92.08, 0.741; GCF_002900365.1, s__Escherichia marmotae, 95.0, 91.17, 0.711; GCF_000759775.1, s__Escherichia albertii, 95.0, 90.64, 0.611      86.0    11      N/AN/A
SRR9988205_bin.2.permissive     d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis        GCF_000025985.1 95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis        99.01   0.901   GCF_000025985.1     95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis        99.01   0.901   d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__    taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments    GCA_947646015.1, s__Bacteroides sp947646015, 95.0, 97.98, 0.349; GCF_019583405.1, s__Bacteroides fragilis_B, 95.0, 91.9, 0.692; GCF_014288095.1, s__Bacteroides hominis, 95.0, 88.11, 0.582     90.21   11      N/A     N/A
(base) [dorian.rojas@accessnode test]$
```

The results here showed the classification for all the four analyzed bins. Two of these were identified to the species level (*Escherichia coli* and *Bacteroides fragilis*), while the other remained as genera (*Treponema* sp. and *Cryptobacteroides* sp.).

These results provide a taxonomical interpretation to the sequences. However, further analysis is still required. For instance, the abundance estimation, sequencing depth, metabolic potential, and other clinical aspects (Antibiotics Resistance Genes, Mobile Genentic Elements) can still be analyzed.

Prior to these, the dereplication of the MAGs can be performed.

## MAGs deprelication (dRep)

MAGs dereplication is an optional process in the metagenomics pipeline analysis. This part of the workshop is going to be demonstrative as the dereplication process can't be correctly performed in the working dataset.

Dereplication consists in analyzing the bins to remove all those sequences that have high similarity among them. These similarities are calculated based on different methods; for instance, ANI values. As mentioned above, this step is optional and in some cases it is skipped to further obtained the desired results (e.g. strain-level studies).

This is performed in order to deal with computational demanding analysis, diversity issues (various genomes from the same species), and unspecific mapping while determining abundance and sequencing depth. The final output is a list of so called representative genomes.

There are several tools that can be used to performed the dereplication process. However, the most widely used is [dRep](https://github.com/MrOlm/drep). dRep employs a bi-phasic approach to deprelicate genomes by using gANI and Mash software. gANI is a very efficient method for dereplication; however, it would take a long time to analyze large genome sets. Moreover, Mash is time-efficient, but compromises the accuracy of the dereplication when completeness decreases [(Olm et al. 2017)](https://academic.oup.com/ismej/article/11/12/2864/7537826).

In this sense, the bi-phasic approach first divides genomes into clusters using Mash. Then, these are compaired in a pairwise manner using gANI generating secondary clusters of near-identical genomes that are dereplicated [(Olm et al. 2017)](https://academic.oup.com/ismej/article/11/12/2864/7537826).

> Notice this approach does not considers the taxonomical classification of the bins and it is merely based on statistics. This allows the analysis to be performed either before or after the taxonomical assigment with GTDB-tk.

Here, a basic dereplication process for an example dataset of bins is conducted.

### Running dRep

dRep is installed as a container. The first help message is a menu to select the compare or dereplicate functions of the tools. In this case, the dereplication is our option of interest.

```bash
(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/drep-3.5.0.sif dRep -h

                ...::: dRep v3.5.0 :::...

  Matt Olm. MIT License. Banfield Lab, UC Berkeley. 2017 (last updated 2024)

  See https://drep.readthedocs.io/en/latest/index.html for documentation
  Choose one of the operations below for more detailed help.

  Example: dRep dereplicate -h

  Commands:
    compare            -> Compare and cluster a set of genomes
    dereplicate        -> De-replicate a set of genomes
    check_dependencies -> Check which dependencies are properly installed


(base) [dorian.rojas@accessnode test]$ /opt/ohpc/pub/containers/BIO/drep-3.5.0.sif dRep dereplicate -h
usage: dRep dereplicate [-p PROCESSORS] [-d] [-h] [-g [GENOMES ...]] [-l LENGTH] [-comp COMPLETENESS]
                        [-con CONTAMINATION] [--ignoreGenomeQuality] [--genomeInfo GENOMEINFO]
                        [--checkM_method {lineage_wf,taxonomy_wf}] [--set_recursion SET_RECURSION]
                        [--checkm_group_size CHECKM_GROUP_SIZE]
                        [--S_algorithm {ANIn,goANI,fastANI,skani,ANImf,gANI}] [-ms MASH_SKETCH]
                        [--SkipMash] [--SkipSecondary] [--skani_extra SKANI_EXTRA]
                        [--n_PRESET {normal,tight}] [-pa P_ANI] [-sa S_ANI] [-nc COV_THRESH]
                        [-cm {total,larger}]
                        [--clusterAlg {complete,single,average,ward,centroid,median,weighted}]
                        [--multiround_primary_clustering] [--primary_chunksize PRIMARY_CHUNKSIZE]
                        [--greedy_secondary_clustering] [--run_tertiary_clustering]
                        [-comW COMPLETENESS_WEIGHT] [-conW CONTAMINATION_WEIGHT]
                        [-strW STRAIN_HETEROGENEITY_WEIGHT] [-N50W N50_WEIGHT] [-sizeW SIZE_WEIGHT]
                        [-centW CENTRALITY_WEIGHT] [-extraW EXTRA_WEIGHT_TABLE] [--gen_warnings]
                        [--warn_dist WARN_DIST] [--warn_sim WARN_SIM] [--warn_aln WARN_ALN] [--skip_plots]
                        work_directory

positional arguments:
  work_directory        Directory where data and output are stored
                        *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

SYSTEM PARAMETERS:
  -p PROCESSORS, --processors PROCESSORS
                        threads (default: 6)
  -d, --debug           make extra debugging output (default: False)
  -h, --help            show this help message and exit

GENOME INPUT:
  -g [GENOMES ...], --genomes [GENOMES ...]
                        genomes to filter in .fasta format. Not necessary if Bdb or Wdb already exist. Can
                        also input a text file with paths to genomes, which results in fewer OS issues
                        than wildcard expansion (default: None)

GENOME FILTERING OPTIONS:
  -l LENGTH, --length LENGTH
                        Minimum genome length (default: 50000)
  -comp COMPLETENESS, --completeness COMPLETENESS
                        Minimum genome completeness (default: 75)
  -con CONTAMINATION, --contamination CONTAMINATION
                        Maximum genome contamination (default: 25)

GENOME QUALITY ASSESSMENT OPTIONS:
  --ignoreGenomeQuality
                        Dont run checkM or do any quality filtering. NOT RECOMMENDED! This is useful for
                        use with bacteriophages or eukaryotes or things where checkM scoring does not
                        work. Will only choose genomes based on length and N50 (default: False)
  --genomeInfo GENOMEINFO
                        location of .csv file containing quality information on the genomes. Must contain:
                        ["genome"(basename of .fasta file of that genome), "completeness"(0-100 value for
                        completeness of the genome), "contamination"(0-100 value of the contamination of
                        the genome)] (default: None)
  --checkM_method {lineage_wf,taxonomy_wf}
                        Either lineage_wf (more accurate) or taxonomy_wf (faster) (default: lineage_wf)
  --set_recursion SET_RECURSION
                        Increases the python recursion limit. NOT RECOMMENDED unless checkM is crashing
                        due to recursion issues. Recommended to set to 2000 if needed, but setting this
                        could crash python (default: 0)
  --checkm_group_size CHECKM_GROUP_SIZE
                        The number of genomes passed to checkM at a time. Increasing this increases RAM
                        but makes checkM faster (default: 2000)

GENOME COMPARISON OPTIONS:
  --S_algorithm {ANIn,goANI,fastANI,skani,ANImf,gANI}
                        Algorithm for secondary clustering comaprisons:
                        fastANI = Kmer-based approach; very fast
                        skani = Even faster Kmer-based approacht
                        ANImf   = (DEFAULT) Align whole genomes with nucmer; filter alignment; compare aligned regions
                        ANIn    = Align whole genomes with nucmer; compare aligned regions
                        gANI    = Identify and align ORFs; compare aligned ORFS
                        goANI   = Open source version of gANI; requires nsmimscan
                         (default: fastANI)
  -ms MASH_SKETCH, --MASH_sketch MASH_SKETCH
                        MASH sketch size (default: 1000)
  --SkipMash            Skip MASH clustering, just do secondary clustering on all genomes (default: False)
  --SkipSecondary       Skip secondary clustering, just perform MASH clustering (default: False)
  --skani_extra SKANI_EXTRA
                        Extra arguments to pass to skani triangle (default: )
  --n_PRESET {normal,tight}
                        Presets to pass to nucmer
                        tight   = only align highly conserved regions
                        normal  = default ANIn parameters (default: normal)

GENOME CLUSTERING OPTIONS:
  -pa P_ANI, --P_ani P_ANI
                        ANI threshold to form primary (MASH) clusters (default: 0.9)
  -sa S_ANI, --S_ani S_ANI
                        ANI threshold to form secondary clusters (default: 0.95)
  -nc COV_THRESH, --cov_thresh COV_THRESH
                        Minmum level of overlap between genomes when doing secondary comparisons (default:
                        0.1)
  -cm {total,larger}, --coverage_method {total,larger}
                        Method to calculate coverage of an alignment
                        (for ANIn/ANImf only; gANI and fastANI can only do larger method)
                        total   = 2*(aligned length) / (sum of total genome lengths)
                        larger  = max((aligned length / genome 1), (aligned_length / genome2))
                         (default: larger)
  --clusterAlg {complete,single,average,ward,centroid,median,weighted}
                        Algorithm used to cluster genomes (passed to scipy.cluster.hierarchy.linkage
                        (default: average)

GREEDY CLUSTERING OPTIONS
These decrease RAM use and runtime at the expense of a minor loss in accuracy.
Recommended when clustering 5000+ genomes:
  --multiround_primary_clustering
                        Cluster each primary clunk separately and merge at the end with single linkage.
                        Decreases RAM usage and increases speed, and the cost of a minor loss in precision
                        and the inability to plot primary_clustering_dendrograms. Especially helpful when
                        clustering 5000+ genomes. Will be done with single linkage clustering (default:
                        False)
  --primary_chunksize PRIMARY_CHUNKSIZE
                        Impacts multiround_primary_clustering. If you have more than this many genomes,
                        process them in chunks of this size. (default: 5000)
  --greedy_secondary_clustering
                        Use a heuristic to avoid pair-wise comparisons when doing secondary clustering.
                        Will be done with single linkage clustering. Only works for fastANI S_algorithm
                        option at the moment (default: False)
  --run_tertiary_clustering
                        Run an additional round of clustering on the final genome set. This is especially
                        useful when greedy clustering is performed and/or to handle cases where similar
                        genomes end up in different primary clusters. Only works with dereplicate, not
                        compare. (default: False)

SCORING CRITERIA
Based off of the formula:
A*Completeness - B*Contamination + C*(Contamination * (strain_heterogeneity/100)) + D*log(N50) + E*log(size) + F*(centrality - S_ani)

A = completeness_weight; B = contamination_weight; C = strain_heterogeneity_weight; D = N50_weight; E = size_weight; F = cent_weight:
  -comW COMPLETENESS_WEIGHT, --completeness_weight COMPLETENESS_WEIGHT
                        completeness weight (default: 1)
  -conW CONTAMINATION_WEIGHT, --contamination_weight CONTAMINATION_WEIGHT
                        contamination weight (default: 5)
  -strW STRAIN_HETEROGENEITY_WEIGHT, --strain_heterogeneity_weight STRAIN_HETEROGENEITY_WEIGHT
                        strain heterogeneity weight (default: 1)
  -N50W N50_WEIGHT, --N50_weight N50_WEIGHT
                        weight of log(genome N50) (default: 0.5)
  -sizeW SIZE_WEIGHT, --size_weight SIZE_WEIGHT
                        weight of log(genome size) (default: 0)
  -centW CENTRALITY_WEIGHT, --centrality_weight CENTRALITY_WEIGHT
                        Weight of (centrality - S_ani) (default: 1)
  -extraW EXTRA_WEIGHT_TABLE, --extra_weight_table EXTRA_WEIGHT_TABLE
                        Path to a tab-separated file with two-columns, no headers, listing genome and
                        extra score to apply to that genome (default: None)

WARNINGS:
  --gen_warnings        Generate warnings (default: False)
  --warn_dist WARN_DIST
                        How far from the threshold to throw cluster warnings (default: 0.25)
  --warn_sim WARN_SIM   Similarity threshold for warnings between dereplicated genomes (default: 0.98)
  --warn_aln WARN_ALN   Minimum aligned fraction for warnings between dereplicated genomes (ANIn)
                        (default: 0.25)

ANALYZE:
  --skip_plots          Dont make plots (default: False)

Example: dRep dereplicate output_dir/ -g /path/to/genomes/*.fasta
```

This software is very complete for different kind of analysis. Therefore, it contains a great quantity of functions that could be selected in the command. However, the main code is composed of the required flags `dRep dereplicate outout_directory -g path/to/genomes/*.fasta`.

The software requires the integration of CheckM data about completeness and contamination which could favor the accuracy of the process. This can be created through a simple python script.

> This code is an modified example provided by the M.Sc Maria Alejandra Soto, a secondee from the project, for the demostrative purposes of the workshop.

```python
import glob
import pandas as pd
from tqdm.auto import tqdm
import re
import os
import shutil


# Extract list of genomes (MAGs) not classified at the species level
ale_tax_classification = pd.read_csv('/home/dorian.rojas/1-test-alejandra/0-data/test_drep/mags_taxonomy_gtdb.csv')
ale_tax_classification['classification'] = ale_tax_classification[['domain','phylum','class','order','family','genus','species']].agg(';'.join, axis = 1)
ale_tax_classification['user_genome'] = ale_tax_classification['user_genome'] + '_filtered_kept_contigs.fasta.gz'
ale_for_drep = ale_tax_classification[['user_genome', 'classification']]
 
example_tax_classification = pd.read_csv('/home/dorian.rojas/1-test-alejandra/10-gtdbtk/gtdbtk.bac120.summary.tsv', delimiter='\t')
example_tax_classification['user_genome'] = example_tax_classification['user_genome'] + '_filtered_kept_contigs.fasta.gz'
example_for_drep = example_tax_classification[['user_genome', 'classification']]

# Saving as a new file in the drep folder 
all_for_drep = pd.concat([example_for_drep, ale_for_drep], ignore_index=True)
all_for_drep.to_csv('/home/dorian.rojas/1-test-alejandra/11-dRep/for_drep.csv', index =False)

###################################### optional ###########################################
# Create a directory with the genomes files present in the previously generated list
directories_to_search = ['/home/dorian.rojas/1-test-alejandra/6-mdmcleaner/bins', '/home/dorian.rojas/1-test-alejandra/0-data/test_drep/data']

# Load the list of genomes to analize 
for_drep = pd.read_csv('/home/dorian.rojas/1-test-alejandra/11-dRep/for_drep.csv')

# Create an output directory for the copied files
  # I copied the files for organization purposes (notice it was created in the 9-final-MAGs folder). However, it is also recommended in case not all files were created by your username (e.g. public databases), this avoids permission issues. 
   # This script can be modified so it does not copy the files but instead only creates the .csv file.
output_directory = '/home/dorian.rojas/1-test-alejandra/9-final-MAGs/all_genomes'
os.makedirs(output_directory, exist_ok=True)

i = 0
# Iterate through the directories
for directory in directories_to_search:
    # Recursively walk through the directory and its subdirectories
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith('_filtered_kept_contigs.fasta.gz'):
                # Check if the filename (with extension) is in the 'user_genome' column
                if filename in for_drep['user_genome'].values:
                    i += 1
                    print(f'File number {i} found: {filename}!')
                    filepath = os.path.join(root, filename)
                    # Copy the file to the output directory
                    shutil.copy(filepath, os.path.join(output_directory, filename))

# Create a CSV file with the paths of copied files
copied_files_df = pd.DataFrame({'file_path': [os.path.join(output_directory, f) for f in os.listdir(output_directory)]})
copied_files_df.to_csv('/home/dorian.rojas/1-test-alejandra/11-dRep/paths_genomes_for_drep.csv', index=False)


genomes_to_find = pd.read_csv('/home/dorian.rojas/1-test-alejandra/11-dRep/for_drep.csv')

# Changing column names for unification (genome, completness, contamination)
checkm_sample = pd.read_csv('/home/dorian.rojas/1-test-alejandra/9-final-MAGs/final_mags.tsv', sep = '\t')
checkm_sample['genome'] = checkm_sample['genome'].str.replace(r'_filtered_kept_contigs\.fasta$', '_filtered_kept_contigs.fasta.gz', regex=True)
checkm_sample = checkm_sample.rename(columns = {'genome':'user_genome'}) 
sample_checkm2_for_drep = checkm_sample[checkm_sample['user_genome'].isin(genomes_to_find['user_genome'])]

ale_checkm = pd.read_csv('/home/dorian.rojas/1-test-alejandra/0-data/test_drep/mags_checkm2.csv')
ale_checkm['genome'] = ale_checkm['genome'].str.replace(r'_filtered_kept_contigs\.fasta$', '_filtered_kept_contigs.fasta.gz', regex=True)
ale_checkm = ale_checkm.rename(columns = {'genome':'user_genome'})
ale_checkm_for_drep = ale_checkm[ale_checkm['user_genome'].isin(genomes_to_find['user_genome'])]

# Creating consolidated file with all the genomes and their quality metrics
all_checkm2_for_drep = pd.concat([ale_checkm_for_drep[['user_genome', 'Completeness', 'Contamination']], sample_checkm2_for_drep[['user_genome','Completeness','Contamination']]], ignore_index = True)

# Renaming columns to be identified by dRep
all_checkm2_for_drep = all_checkm2_for_drep.rename(columns = {'user_genome':'genome','Completeness':'completeness','Contamination':'contamination'})
all_checkm2_for_drep.to_csv('/home/dorian.rojas/1-test-alejandra/11-dRep/all_checkm2_for_drep.csv', index=False)
```

The final output are several files in the `11-dRep` folder. The `all_checkm2_for_drep.csv` and `paths_genomes_for_drep.csv` are the inputs for the dRep analysis, it contains the quality information from CheckM2 and the path to the `.fasta` files of each bin.

```bash
(base) [dorian.rojas@accessnode 1-test-alejandra]$ ls 11-dRep/
all_checkm2_for_drep.csv  for_drep.csv  paths_genomes_for_drep.csv

(base) [dorian.rojas@accessnode 1-test-alejandra]$ head 11-dRep/all_checkm2_for_drep.csv 
genome,completeness,contamination
49647_1#2_clean_bin.19.orig_filtered_kept_contigs.fasta.gz,88.21,1.12
```

Finally, the `drep.slurm` file can be completed. For this analysis we stabilish various flags. These include: `-g` to indicate the path files, `-p` for the threats, `-comp` for minimum completeness of 50, `-con` for maximum contamination of 5, `--S_ani` for 0.95 ANI threshold, `--cov_thresh` for 0.30 coverage threshold, `--S_algorithm` for ANI algorithm to use, and the `--genomeInfo` for the CheckM2 information.

> Notice the `--S_algorithm` selects the fastANI approach althought the ideal of dRep is the use of gANI. However, the fastANI is faster and works best for the demonstrative purpose of the workshop.

Similarly to the GTDB-tk code, this ignores the `batch.sh` file and the `for` loop and it is sent directly with all the samples in the `paths_genomes_for_drep.csv` file using `sbatch drep.slurm`. The final example of the slurm file used to run the code for dereplication is presented below.

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="drep"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL

cd /home/dorian.rojas/1-test-alejandra

CTN_PATH=/opt/ohpc/pub/containers/BIO/
$CTN_PATH/drep-3.5.0.sif dRep --version

mkdir -p 11-dRep/

$CTN_PATH/drep-3.5.0.sif dRep dereplicate ./11-dRep -g  11-dRep/paths_genomes_for_drep.csv \
        -p $SLURM_NTASKS -comp 50 -con 5 --genomeInfo 11-dRep/all_checkm2_for_drep.csv \
        --S_ani 0.95 --cov_thresh 0.30 --S_algorithm fastANI

date
time
```

The output from the dRep analysis are several files and folders within the output folder.

```bash
(base) [dorian.rojas@accessnode 1-test-alejandra]$ ls 11-dRep/
all_checkm2_for_drep.csv  data_tables           figures       log
data                      dereplicated_genomes  for_drep.csv  paths_genomes_for_drep.csv
```

Notice the `.csv` file were created by the python script. Hence, the output files from dRep are within the new folders.

The folder figures contain different images in `.pdf` format that favor the understand of the similarity relationship between the analyzed genomes. These correspond, for example, to the dendogram of the primary and secondary cluster. Below the example fo the dendrogram of the primary clusterng is presented, the dotted line represents the ANI threshold of 0.95,

![Primary Clustering Dendogram](Primary_clustering_dendrogram.jpg)

Analyze the rest of the images and try to understand thei relevance for the analysis. In addition, an image representing the quality score of the cluster is presented with the selected (representative) genome marked with an *. This graph is presented in the `Cluster_scoring.pdf` file. This allows the user to easily identify that dRep selected the best genome

![Cluster Scoring](Cluster_scoring.png)

Moreover, the `data_tables` folders contains the information generated to select the representatives genomes. The `log` folder contains the log files. The `data` folder includes the data greated during the workflo2 (e.g. ANI values). Finally, the `dereplicated_genomes` holds the sequences files of the representative genomes.

```bash
(base) [dorian.rojas@accessnode 11-dRep]$ ls data_tables/
Bdb.csv  Cdb.csv  genomeInfo.csv  genomeInformation.csv  Mdb.csv  Ndb.csv  Sdb.csv  Wdb.csv  Widb.csv

(base) [dorian.rojas@accessnode 11-dRep]$ ls log/
cluster_arguments.json  logger.log

(base) [dorian.rojas@accessnode 11-dRep]$ ls dereplicated_genomes/
49647_1#100_clean_bin.11_filtered_kept_contigs.fasta.gz
49647_1#100_clean_bin.14_filtered_kept_contigs.fasta.gz
49647_1#100_clean_bin.18_filtered_kept_contigs.fasta.gz
49647_1#100_clean_bin.1_filtered_kept_contigs.fasta.gz
[...]
```

The final results from this `dereplicated_genomes` are the one used for the abundance estimation.

## Functional annotation of MAGs (ABRicate, eggNOG-mapper)

Further analysis is performed using a pangenome. A pangenome is a file with multiple sequences representing the genomes from the sample. Basically, it is a concatenation of the bins into a singular file. Although, the individual files can be used in the analysis, the utilization of a pangenome could avoid several coding issues, faciliting the overall pipeline. This pangenome can be created through different approaches. Here, the provide a simple bash code using the command `cat`.

>Notice the following sections of the workshop are hands-on. Remember where are your files of interest located in your console.

```bash
#!/bin/bash

cd /home/dorian.rojas/test
CTN_PATH=/opt/ohpc/pub/containers/BIO/

# Create a general directpry
mkdir -p 10-final-MAGs/bins

# Select each of the bins that passed the GUNC analysis
for file_name in $(awk '{print $2}' 10-final-MAGs/gtdb-input.txt); do
  # Create a variable for the name of the (de)compressed bins
  compressed="${file_name}_filtered_kept_contigs.fasta.gz"
  newname="${file_name}_filtered_kept_contigs.fasta"

  # Decompress the bins
  gzip -dc 7-mdmcleaner/bins/$compressed > 7-mdmcleaner/bins/$newname

  # Rename the contigs using the bbmap toolkit and output in general directory
  $CTN_PATH/bbmap-39.10.sif rename.sh in=7-mdmcleaner/bins/$newname \
  out=10-final-MAGs/bins/$newname prefix=$file_name
done 

# Concatenate decompressed renamed files
cat 10-final-MAGs/bins/*.fasta > 10-final-MAGs/bins/pangenome.fasta
```

Copy the code into a `pangenome.sh` file. Ensure to change the path where the data is located. This script will decompress the `.fasta.gz` files and rename each contigs with the sample name. This step is crucial as `cat` will combine all contigs regardless of the their names. For this, the toolkit BBMap and module `rename.sh` are used. Finally, all the decompressed renamed files are combined into a single pangenome fasta.

This final file can be use to facilitate the command of other analysis. During this workshop, the functional annotation of clinical features and metabolic pathways are going to be performed. First the annotation of clinical features refers to the identification of genomic characteristics that can be associated to a clinical risk. In this regard, mobile genetic elements and antibiotic resistance genes can be of interest.

Several tools have been developed for this approach. For instance, the Comprehensive Antibiotic Resistance Database (CARD) developed the [Resistance Gene Identifer](https://github.com/arpcard/rgi) (RGI) that allows the annotation against this database. Moreover, the [AMRFinderPlus](https://github.com/ncbi/amr) performs a similar tasks againt the National Center for Biotechnology Information (NCBI) database. The large quantity of databases for antibiotic resistance genes and other features resulted in the need of a more general software that could conduct the identification against multiple databases. As a solution to this issue, ABRicate was proposed.

[ABRicate](https://github.com/tseemann/abricate/tree/master) is a complete tool to perform the annotation of different features taking in consideration multiple databases. Among the predetermine reference, the NCBI, CARD, Resfinder, ARG-ANNOT, MEGARES, EcOH, PlasmidFinder, VFDB, and Ecoli_VF are included. The tools also allows the user to create they custom databases for more specific analysis. To perform the annotation, ABRicate employs the BLAST+ software.

> The authors defined the etymology as "Anti-Biotic Resistance" in the form of an English verb to represent the software taking action against the antibiotic resistance. They also stated that "[it] is unlikely to receive an infamous [JABBA AWARD](http://www.acgt.me/blog/2014/12/1/time-for-a-new-jabba-award-for-just-another-bogus-bioinformatics-acronym) (Just Another Bogus Bioinformatics Acronym Award).

### Running abricate

This software is installed as a container and is simple to use. The only issue is that in order run multiple databases at the same time, a bash `for` loop must be used as the toolf does not accept various databases simultaneously.

```bash
```

Code for a simple annotation using all the available databases by the software.

The output from the analysis is a `.tsv` table with several columns that indicate the name of the file, contigs in which the hit was detected, gene, coverage, database, among others. Notice the use of the pangenome file for this analysis facilitates the running code. However, the rename of the contigs headers was required to further associated the identified hits to a particular bin.

### Running eggnog-mapper

[]

## Code templates

**`gtdbtk.slurm` file:**

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="gtdbtk"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL

cd /home/dorian.rojas/test

CTN_PATH=/opt/ohpc/pub/containers/BIO/
DBPATH=/home/dorian.rojas/DB/CheckM2_database/uniref100.KO.1.dmnd

mkdir -p 11-gtdbtk/

export GTDBTK_DATA_PATH="/home/dorian.rojas/DB/release220/"

$CTN_PATH/gtdbtk-2.4.0.sif gtdbtk classify_wf --cpus $SLURM_NTASKS \
         --batchfile 10-final-MAGs/gtdb-input.txt --out_dir 11-gtdbtk/ \
        -x .gz --skip_ani_screen


date
time
```

**`drep-genomeinfo.py` file:**

```python
import glob
import pandas as pd
from tqdm.auto import tqdm
import re
import os
import shutil


# Extract list of genomes (MAGs) not classified at the species level
ale_tax_classification = pd.read_csv('/home/dorian.rojas/1-test-alejandra/0-data/test_drep/mags_taxonomy_gtdb.csv')
ale_tax_classification['classification'] = ale_tax_classification[['domain','phylum','class','order','family','genus','species']].agg(';'.join, axis = 1)
ale_tax_classification['user_genome'] = ale_tax_classification['user_genome'] + '_filtered_kept_contigs.fasta.gz'
ale_for_drep = ale_tax_classification[['user_genome', 'classification']]
 
example_tax_classification = pd.read_csv('/home/dorian.rojas/1-test-alejandra/10-gtdbtk/gtdbtk.bac120.summary.tsv', delimiter='\t')
example_tax_classification['user_genome'] = example_tax_classification['user_genome'] + '_filtered_kept_contigs.fasta.gz'
example_for_drep = example_tax_classification[['user_genome', 'classification']]

# Saving as a new file in the drep folder 
all_for_drep = pd.concat([example_for_drep, ale_for_drep], ignore_index=True)
all_for_drep.to_csv('/home/dorian.rojas/1-test-alejandra/11-dRep/for_drep.csv', index =False)

###################################### optional ###########################################
# Create a directory with the genomes files present in the previously generated list
directories_to_search = ['/home/dorian.rojas/1-test-alejandra/6-mdmcleaner/bins', '/home/dorian.rojas/1-test-alejandra/0-data/test_drep/data']

# Load the list of genomes to analize 
for_drep = pd.read_csv('/home/dorian.rojas/1-test-alejandra/11-dRep/for_drep.csv')

# Create an output directory for the copied files
  # I copied the files for organization purposes (notice it was created in the 9-final-MAGs folder). However, it is also recommended in case not all files were created by your username (e.g. public databases), this avoids permission issues. 
   # This script can be modified so it does not copy the files but instead only creates the .csv file.
output_directory = '/home/dorian.rojas/1-test-alejandra/9-final-MAGs/all_genomes'
os.makedirs(output_directory, exist_ok=True)

i = 0
# Iterate through the directories
for directory in directories_to_search:
    # Recursively walk through the directory and its subdirectories
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith('_filtered_kept_contigs.fasta.gz'):
                # Check if the filename (with extension) is in the 'user_genome' column
                if filename in for_drep['user_genome'].values:
                    i += 1
                    print(f'File number {i} found: {filename}!')
                    filepath = os.path.join(root, filename)
                    # Copy the file to the output directory
                    shutil.copy(filepath, os.path.join(output_directory, filename))

# Create a CSV file with the paths of copied files
copied_files_df = pd.DataFrame({'file_path': [os.path.join(output_directory, f) for f in os.listdir(output_directory)]})
copied_files_df.to_csv('/home/dorian.rojas/1-test-alejandra/11-dRep/paths_genomes_for_drep.csv', index=False)


genomes_to_find = pd.read_csv('/home/dorian.rojas/1-test-alejandra/11-dRep/for_drep.csv')

# Changing column names for unification (genome, completness, contamination)
checkm_sample = pd.read_csv('/home/dorian.rojas/1-test-alejandra/9-final-MAGs/final_mags.tsv', sep = '\t')
checkm_sample['genome'] = checkm_sample['genome'].str.replace(r'_filtered_kept_contigs\.fasta$', '_filtered_kept_contigs.fasta.gz', regex=True)
checkm_sample = checkm_sample.rename(columns = {'genome':'user_genome'}) 
sample_checkm2_for_drep = checkm_sample[checkm_sample['user_genome'].isin(genomes_to_find['user_genome'])]

ale_checkm = pd.read_csv('/home/dorian.rojas/1-test-alejandra/0-data/test_drep/mags_checkm2.csv')
ale_checkm['genome'] = ale_checkm['genome'].str.replace(r'_filtered_kept_contigs\.fasta$', '_filtered_kept_contigs.fasta.gz', regex=True)
ale_checkm = ale_checkm.rename(columns = {'genome':'user_genome'})
ale_checkm_for_drep = ale_checkm[ale_checkm['user_genome'].isin(genomes_to_find['user_genome'])]

# Creating consolidated file with all the genomes and their quality metrics
all_checkm2_for_drep = pd.concat([ale_checkm_for_drep[['user_genome', 'Completeness', 'Contamination']], sample_checkm2_for_drep[['user_genome','Completeness','Contamination']]], ignore_index = True)

# Renaming columns to be identified by dRep
all_checkm2_for_drep = all_checkm2_for_drep.rename(columns = {'user_genome':'genome','Completeness':'completeness','Contamination':'contamination'})
all_checkm2_for_drep.to_csv('/home/dorian.rojas/1-test-alejandra/11-dRep/all_checkm2_for_drep.csv', index=False)
```

**`drep.slurm` file:**

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="drep"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL

cd /home/dorian.rojas/1-test-alejandra

CTN_PATH=/opt/ohpc/pub/containers/BIO/
$CTN_PATH/drep-3.5.0.sif dRep --version

mkdir -p 11-dRep/

$CTN_PATH/drep-3.5.0.sif dRep dereplicate ./11-dRep -g  11-dRep/paths_genomes_for_drep.csv \
        -p $SLURM_NTASKS -comp 50 -con 5 --genomeInfo 11-dRep/all_checkm2_for_drep.csv \
        --S_ani 0.95 --cov_thresh 0.30 --S_algorithm fastANI

date
time
```

**`pangenome.sh` file:**

```bash
#!/bin/bash

cd /home/dorian.rojas/test
CTN_PATH=/opt/ohpc/pub/containers/BIO/

# Create a general directpry
mkdir -p 10-final-MAGs/bins

# Select each of the bins that passed the GUNC analysis
for file_name in $(awk '{print $2}' 10-final-MAGs/gtdb-input.txt); do
  # Create a variable for the name of the (de)compressed bins
  compressed="${file_name}_filtered_kept_contigs.fasta.gz"
  newname="${file_name}_filtered_kept_contigs.fasta"

  # Decompress the bins
  gzip -dc 7-mdmcleaner/bins/$compressed > 7-mdmcleaner/bins/$newname

  # Rename the contigs using the bbmap toolkit and output in general directory
  $CTN_PATH/bbmap-39.10.sif rename.sh in=7-mdmcleaner/bins/$newname \
  out=10-final-MAGs/bins/$newname prefix=$file_name
done 

# Concatenate decompressed renamed files
cat 10-final-MAGs/bins/*.fasta > 10-final-MAGs/bins/pangenome.fasta
```

**`abricate.slurm` file:**


**`eggnog-mapper.slurm` file:**

