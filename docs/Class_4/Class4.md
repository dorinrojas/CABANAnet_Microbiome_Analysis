# Class 4: Functional and taxonomical annotation of MAGs

- - - -

## Taxonomical annotation of MAGs (GTDBtk)

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

The database for the gtdbtk must be indicated by adding `export GTDBTK_DATA_PATH="/home/public/met-workshop/databases/release220/"` prior the tool command. `export` is a base command that sets an environmental variable. Most variables are different depending on the running software; therefore, this information has to be retrieved from the source code or the github repository.

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

> This code is an modified example provided by the M.Sc. Maria Alejandra Soto, a secondee from the project, for the demostrative purposes of the workshop.

```python
import glob
import pandas as pd
from tqdm.auto import tqdm
import re
import os
import shutil


# Extract list of genomes (MAGs)
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
49647_1#102_clean_bin.12_filtered_kept_contigs.fasta.gz,65.76,0.53
49647_1#106_clean_bin.26_24_filtered_kept_contigs.fasta.gz,99.42,2.34
49647_1#102_clean_bin.11_filtered_kept_contigs.fasta.gz,84.66,3.82
49647_1#102_clean_bin.5_filtered_kept_contigs.fasta.gz,91.64,3.67
49647_1#100_clean_bin.5_filtered_kept_contigs.fasta.gz,99.43,0.35
49647_1#104_clean_bin.1_4_filtered_kept_contigs.fasta.gz,93.69,2.12
49647_1#106_clean_bin.10_33_filtered_kept_contigs.fasta.gz,96.18,0.18
49647_1#106_clean_bin.12_23_filtered_kept_contigs.fasta.gz,60.96,1.27
49647_1#102_clean_bin.3_filtered_kept_contigs.fasta.gz,98.84,0.24
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

The final results from this `dereplicated_genomes` are the one used for the abundance estimation and functional annotation.

As previously mentioned, this part of the workshop is demonstrative as the datasets does not allow to perform dereplciation (due to the presence of only a few bins that do not require dereplication). Therefore, the following section will continue to be performed with the working samples.

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

Moreover, another interesting feature to analyze in the microbiome is the potential metabolic processes encoded in the identified bacteria. These are considered 'potential' as the genomic data can not ensure that these are being expressed.

Similarly to the clinical feature, a vast variety of tools with different methodologies have been developed to perform metabolic annotations. For instance, the software [HUMAnN](https://github.com/biobakery/humann) (HMP Unified Metabolic Analysis Network) has been widely used in this regard. This is an efficient tools for profiling the abundace of microbial metabolic pathways and other molecular functions from metagenomics or metatranscriptomics data.

During this workshop, the software [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-Basic_usage) will be used. This software performs annotation based on its own database of precomputed orthologous groups (OG) and phylogenies. Such approach, avoids biases based on taxonomic profiling. For this, fasta files are translated into protens using Prodigal and the output is aligned to the eggNOG database using either HMMER3, DIAMOND, or MMseqs2. The orthology is infered based on the selected taxonomy (eukaryota, bacteria, archaea, virus) and finally OGs are annotated based on different databases (e.g. GeneOntology, Pfam, KEGG...).

>Interestingly, eggNOG-mapper seems to have a higher annotation rate in comparison to HUMAnN, as it allows the identification based on OGs instead of species identification [(Mi et al. 2024)](https://www.sciencedirect.com/science/article/pii/S2667237524003229). This allows the inclusion of novel, uncultured, poorly represented species, among others.

### Running abricate

This software is installed as a container and is simple to use. The only issue is that in order run multiple databases at the same time, a bash `for` loop must be used as the toolf does not accept various databases simultaneously.

>Please ignore a perl error if it appears. This should not affect the pipeline.

```bash
(base) [dorian.rojas@accessnode 1-test-alejandra]$ /opt/ohpc/pub/containers/BIO/abricate-1.0.1.sif abricate -h
SYNOPSIS
  Find and collate amplicons in assembled contigs
AUTHOR
  Torsten Seemann (@torstenseemann)
USAGE
  % abricate --list
  % abricate [options] <contigs.{fasta,gbk,embl}[.gz] ...> > out.tab
  % abricate [options] --fofn fileOfFilenames.txt > out.tab
  % abricate --summary <out1.tab> <out2.tab> <out3.tab> ... > summary.tab
GENERAL
  --help          This help.
  --debug         Verbose debug output.
  --quiet         Quiet mode, no stderr output.
  --version       Print version and exit.
  --check         Check dependencies are installed.
  --threads [N]   Use this many BLAST+ threads [1].
  --fofn [X]      Run on files listed in this file [].
DATABASES
  --setupdb       Format all the BLAST databases.
  --list          List included databases.
  --datadir [X]   Databases folder [/usr/local/db].
  --db [X]        Database to use [ncbi].
OUTPUT
  --noheader      Suppress column header row.
  --csv           Output CSV instead of TSV.
  --nopath        Strip filename paths from FILE column.
FILTERING
  --minid [n.n]   Minimum DNA %identity [80].
  --mincov [n.n]  Minimum DNA %coverage [80].
MODE
  --summary       Summarize multiple reports into a table.
DOCUMENTATION
  https://github.com/tseemann/abricate
```

The available databases to be used with abricate can be explored with the `--list` command.

```bash
(base) [dorian.rojas@accessnode 1-test-alejandra]$ /opt/ohpc/pub/containers/BIO/abricate-1.0.1.sif abricate --list
DATABASE        SEQUENCES       DBTYPE  DATE
argannot        2223    nucl    2024-Dec-15
card    2631    nucl    2024-Dec-15
ecoh    597     nucl    2024-Dec-15
ecoli_vf        2701    nucl    2024-Dec-15
megares 6635    nucl    2024-Dec-15
ncbi    5386    nucl    2024-Dec-15
plasmidfinder   460     nucl    2024-Dec-15
resfinder       3077    nucl    2024-Dec-15
vfdb    2597    nucl    2024-Dec-15
```

Code for a simple annotation using all the available databases by the software. The output should be indicated using `>>`. Hence, it is recommeded to add an additional `abricate summary [...]` after the `for` loop has finished. This will provide a summary of the databases that can be easier to interpret.

> If the job takes a while running, the results in the common repository `/home/public/met-workshop`

The output from the analysis is a `.tab` table with several columns that indicate the name of the file, contigs in which the hit was detected, gene, coverage, database, among others. One file is generated per databaase used, and an additional summary file has the complete information. Notice the use of the pangenome file for this analysis facilitates the running code. However, the rename of the contigs headers was required to further associated the identified hits to a particular bin.

```bash
(base) [dorian.rojas@accessnode test]$ ls 12-abricate/
argannot.tab  card.tab  ecoh.tab  ecoli_vf.tab  megares.tab  ncbi.tab  plasmidfinder.tab  resfinder.tab  summary.tab  vfdb.tab
```

In this case, no annotations were found in the pangenome. This could be expected considering the small bins and subsampled used for this workshop. However, to have an idea of how a positive result from ABRicate looks, find below and example from the example dataset used for the previous dereplication step.

This is the `argannot.tab` file. The `FILE` and `SEQUENCE` indicate the file and contig in which the genomic feature was annotation that starts from base pair shown in the `START` and `END`. The found gene, accession, and product resistance are indicated in the `GENE`, `DATABASE ACCESSION`, and `PRODUCT RESISTANCE` columns, respectively.

```bash
(base) [dorian.rojas@accessnode 1-test-alejandra]$ head 12-abricate/argannot.tab
#FILE   SEQUENCE        START   END     STRAND  GENE    COVERAGE        COVERAGE_MAP    GAPS    %COVERAGE       %IDENTITY       DATABASE ACCESSION       PRODUCT RESISTANCE
11-dRep/for_pangenome/pangenome.fasta   49647_1#102_clean_bin.12_filtered_kept_contigs_c_000000000099   73      1290    +       (MLS)mef(A)      1-1218/1218     =============== 0/0     100.00  95.07   argannot        U70055:314-1531 (MLS)mef(A)
11-dRep/for_pangenome/pangenome.fasta   49647_1#102_clean_bin.12_filtered_kept_contigs_c_000000000113   3848    5053    -       (MLS)Mef(En2)    1-1206/1206     =============== 0/0     100.00  99.83   argannot        NG_047980:101-1306      (MLS)Mef(En2)
11-dRep/for_pangenome/pangenome.fasta   49647_1#102_clean_bin.12_filtered_kept_contigs_c_000000000293   57      1023    -       (Bla)cfxA6       30-996/996      =============== 0/0     97.09   99.79   argannot        GQ342996:798-1793       (Bla)cfxA6
11-dRep/for_pangenome/pangenome.fasta   49647_1#106_clean_bin.26_24_filtered_kept_contigs_c_000000000030        694     2613    +(Tet)tetM       1-1920/1920     =============== 0/0     100.00  97.24   argannot        DQ534550:1451-3370      (Tet)tetM
11-dRep/for_pangenome/pangenome.fasta   NODE_1208_length_21529_cov_6.442349     20882   21376   -       (Tmt)dfrF       1-495/495=============== 0/0     100.00  100.00  argannot        NG_047755:101-595       (Tmt)dfrF
11-dRep/for_pangenome/pangenome.fasta   NODE_1397_length_19353_cov_262.939631   2675    4648    +       (Tet)tetQ       1-1974/1974      =============== 0/0     100.00  97.52   argannot        Z21523:362-2287 (Tet)tetQ
11-dRep/for_pangenome/pangenome.fasta   NODE_17_length_92276_cov_13.958351      5647    6804    -       (Bla)ampH_Ecoli 1-1158/1158      =============== 0/0     100.00  98.70   argannot        AP012030:395554-396711  (Bla)ampH_Ecoli
11-dRep/for_pangenome/pangenome.fasta   NODE_181_length_73190_cov_11.141587     56031   56525   -       (MLS)lnu(C)     1-495/495=============== 0/0     100.00  98.18   argannot        AY928180:1150-1644      (MLS)lnu(C)
11-dRep/for_pangenome/pangenome.fasta   NODE_4_length_196669_cov_13.947424      25855   27159   +       (Bla)AmpC1_Ecoli        1-1305/1305      =============== 0/0     100.00  98.16   argannot        FN649414:2765051-2766355        (Bla)AmpC1_Ecoli
```

The `summary.tab` is similar. In this case, total annotations and all the annotated genes are presented as columns. The rows represent each `.tab` file per database and the coverage found for the respective gene.

```vim
#FILE   NUM_FOUND       (Bla)AmpC1_Ecoli        (Bla)AmpC2_Ecoli        (Bla)Penicillin_Binding_Protein_Ecoli   (Bla)ampH_Ecoli (Bla)cfxA6       (MLS)Mef(En2)   (MLS)lnu(C)     (MLS)mef(A)     (Tet)tetM       (Tet)tetQ       (Tmt)dfrF       ACRA    ACRB    ACRD     ACRE    ACRF    ACRS    AMPH    ANT6    ANT9    ASMA    BACA    BAER    BAES    BCR     BLAEC   CFX     CPXAR   CRP     CTX      CfxA6   DFRF    ECS88_3547      ECs3712 EMRA    EMRB    EMRD    EMRK    EMRR    EMRY    EPTA    EVGS    Escherichia_coli_acrA    Escherichia_coli_ampC   Escherichia_coli_ampC1_beta-lactamase   Escherichia_coli_ampH   Escherichia_coli_emrE   Escherichia_coli_mdfA    GADW    GADX    H-NS    HNS     IncFIA(HI1)_1_HI1       IncFIB(K)_1_Kpn3        IncR_1  IncX1_1 KDPE    LNUA    LNUC     MARA    MARR    MDFA    MDTA    MDTB    MDTC    MDTE    MDTF    MDTG    MDTH    MDTI    MDTJ    MDTK    MDTM    MDTN    MDTO     MDTP    MEFA    MEFE    MPHB    MSBA    MVRC    Mef(En2)        PBP2    PBP4B   PMRF    ROBA    SOXS    TETM    TETQ    UMNK88_238       YOGI    Z0263   Z0265   Z1307   Z2200   Z2201   Z2204   Z2206   aad9    aadE    acrB    acrD    acrE    acrF    acrS     aec17   aec18   aec19   aec22   aec23   aec24   aec25   aec26   aec28   aec29   aec30   aec31   aec32   ant(6)-Ia_3     artj     b2854   b2972   bacA    baeR    baeS    blaEC-18        cadA    cfaA    cfaB    cfaC    cfaD    cfxA6   cfxA6_1 cheA    cheB     cheR    cheW    cheY    cheZ    clpV    cpxA    csgA    csgB    csgC    csgD    csgE    csgF    csgG    dfrF    eaeH    ecpA     ecpB    ecpC    ecpD    ecpR    ehaA    emrA    emrB    emrK    emrR    emrY    entA    entB    entC    entD    entE    entF     entS    epaO    epaP    epaQ    epaR    epaS    eprH    eprI    eprJ    eprK    eptA    espL1   espL3   espL4   espR1   espX1    espX4   espX5   etrA    evgA    evgS    fdeC    fepA    fepB    fepC    fepD    fepE    fepG    fes     fimD    fimF    fimG     fimH    flgA    flgB    flgC    flgD    flgE    flgF    flgG    flgH    flgI    flgJ    flgK    flgL    flgN    flhA    flhB     flhC    flhD    flhE    fliA    fliC-H4 fliE    fliF    fliG    fliH    fliI    fliJ    fliK    fliL    fliM    fliN    fliO     fliP    fliQ    fliR    fliS    fliT    fliY    fliZ    flk     gadW    gadX    gspC    gspD    gspE    gspF    gspG    gspH     gspI    gspJ    gspK    gspL    gspM    hcp     hlyE    hofB    hofC    hofq    ibeC    kdpE    lnu(AN2)        lnu(C)  lnu(C)_1 lnuC    lpfao113        marA    matF    mdf(A)_1        mdtA    mdtB    mdtC    mdtE    mdtF    mdtG    mdtH    mdtM    mdtN     mdtO    mdtP    mef(A)  mef(A)_3        mef(En2)        mel     motA    motB    mphB    msbA    mtfa    nada    nadb    ompA     orgA    orgB    pmrF    ppdD    ppda    ppdb    ppdc    stgC    stgD    tar/cheM        tet(M)  tet(M)_12       tet(Q)  tet(Q)_1 tetM    tetQ    tia     tolC    upaG/ehaG       vgrG    wzm-O9  wzt-O9  yagV/ecpE       yagW/ecpD       yagX/ecpC       yagY/ecpB        yagZ/ecpA       ycbF    ycbQ    ycbR    ycbS    ycbT    ycbU    ycbV    ycfz    ygdb    ygeG    ygeH    yggr    yghg     ykgK/ecpR       yojI
12-abricate/argannot.tab        11      100.00  100.00  100.00  100.00  97.09   100.00  100.00  100.00  100.00  100.00  100.00  ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .
12-abricate/card.tab    52      .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       100.00  .97.09   .       .       .       .       .       .       .       .       .       .       .       100.00  100.00  100.00  100.00  100.00   100.00  .       .       100.00  .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       100.00  .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       100.00  100.00  100.00  100.00  100.00  .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       99.76   100.00  100.00  .       .       .       .       .       ..       .       .       .       .       .       .       .       .       100.00  .       .       .       .       .       .       .100.00  .       .       .       .       .       .       .       100.00  100.00  100.00  100.00  100.00  .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       100.00  .       .       .       ..       .       .       .       100.00  100.00  .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       100.00  100.00  .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       99.26   .       .       .       100.00  .       100.00  .       .100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  .       .       .       100.00  .       .100.00  100.00  .       .       .       .       .       .       100.00  .       .       .       .       .       .       .       ..       .       .       100.00  100.00  .       100.00  .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       100.00
12-abricate/ecoh.tab    3       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       100.00  .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       100.00  100.00  .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .
12-abricate/ecoli_vf.tab        175     .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       100.00  100.00  .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       100.00  .       100.00  100.00  100.00  100.00  100.00  100.00   100.00  .       .       .       .       .       .       .       100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00   100.00  100.00  100.00  100.00  100.00  .       100.00  100.00  100.00  .       .       .       .       100.00  100.00  100.00   100.00  100.00  .       .       100.00  100.00  100.00  100.00  100.00  100.00  100.00  .       99.35   100.00  100.00  100.00   100.00  100.00  100.00  .       100.00  100.00  100.00  100.00  100.00  100.00  100.00  .       .       .       .       .100.00  100.00  100.00  100.00  100.00  99.28   100.00  100.00  100.00  100.00  99.15   99.91   100.00  100.00  99.70   100.00  .100.00  100.00  99.94   100.00  100.00  99.94   100.00  100.00  .       .       .       100.00  100.00  100.00  100.00  100.00  100.00   100.00  100.00  100.00  100.00;100.00   100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00   100.00  100.00  100.00  99.90   100.00  100.00  99.72   100.00  100.00  .       100.00  100.00  100.00  100.00  100.00  100.00   100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  99.90   .       100.00  100.00   100.00  100.00  100.00  100.00  100.00  100.00  99.67   100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00   .       .       .       .       .       100.00  .       100.00  .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       100.00  100.00  .       .       100.00  100.00  100.00  .       100.00   100.00  .       100.00  100.00  100.00  100.00  100.00  100.00  100.00  .       .       .       .       .       .       100.00   .       100.00  100.00  .       .       .       .       .       .       .       100.00  100.00  100.00  100.00  100.00  100.00   100.00  99.87   100.00  98.98   100.00  100.00  100.00  .       .
12-abricate/megares.tab 64      .       .       .       .       .       .       .       .       .       .       .       100.00  100.00   100.00  100.00  100.00  100.00  100.00  100.00  100.00  99.95   99.76   100.00  100.00  100.00  100.00  97.09   100.00;100.00    100.00  99.79;100.00    .       100.00  .       .       100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  ..       .       .       .       .       100.00  100.00  .       100.00  .       .       .       .       99.26   100.00  100.00  100.00   100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00  100.00   99.51   100.00  100.00  100.00  100.00  .       100.00  100.00  100.00  100.00  100.00  100.00  99.85   .       100.00  ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .
12-abricate/ncbi.tab    11      .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .100.00  100.00  .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       100.00  .       .       .       .       .       100.00   .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .100.00  .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       100.00  100.00  .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       99.51   .       100.00  .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       100.00   .       100.00  .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .
12-abricate/plasmidfinder.tab   4       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       100.00  100.00  100.00  100.00  .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .
12-abricate/resfinder.tab       7       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       100.00  .       .       .       .       .       .       .       .       .       .       .       ..       97.09   .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       100.00  .       .       .       .       100.00   .       .       .       .       .       .       .       .       .       .       .       .       100.00  .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       100.00  .       100.00  .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .
12-abricate/vfdb.tab    44      .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       99.78   .       100.00  .       99.04   100.00  ..       .       .       .       .       .       .       .       .       .       .       .       100.00  100.00  100.00  100.00  100.00   99.28   100.00  .       .       .       .       .       .       .       .       .       .       100.00  .       .       83.81    100.00  99.94   100.00  .       .       .       99.76   100.00  100.00  100.00  100.00  .       100.00  100.00  .       100.00   100.00  100.00  .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       100.00  100.00  100.00  100.00  100.00  100.00  100.00   100.00  100.00  100.00  100.00  .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       100.00  .       .       .       .       .       .       .       ..       .       .       .       .       .       .       .       .       .       .       .       .       .       99.74   100.00  100.00   100.00  100.00  .       .       .       .       .       .       .       .       .       .       .       .       .       100.00   .
```

It is important to mention that although the `summary.tab` file might come handy to have a general perspective of what is annotated in the pangenome. This file does not provide an indicative of which sample/bin presents the encoded gene. Therefore, this can not be used to analyze the difference in these genomic features according to the sample. Performing a simple `cat` of all the different `$db.tab` files could be a better approach.

### Running eggnog-mapper

This software is also known as `emapper.py` or just emapper as this is its call command. It is installed as a contains and presents a great variety of options that favor the selection of user-specific preferences for different type of analysis.

```bash
(base) [dorian.rojas@accessnode ~]$ /opt/ohpc/pub/containers/BIO/eggnog-mapper-2.1.12.sif emapper.py -h
usage: emapper.py [-h] [-v] [--list_taxa] [--cpu NUM_CPU] [--mp_start_method {fork,spawn,forkserver}] [--resume] [--override]
                  [-i FASTA_FILE] [--itype {CDS,proteins,genome,metagenome}] [--translate]
                  [--annotate_hits_table SEED_ORTHOLOGS_FILE] [-c FILE] [--data_dir DIR] [--genepred {search,prodigal}]
                  [--trans_table TRANS_TABLE_CODE] [--training_genome FILE] [--training_file FILE]
                  [--allow_overlaps {none,strand,diff_frame,all}] [--overlap_tol FLOAT]
                  [-m {diamond,mmseqs,hmmer,no_search,cache,novel_fams}] [--pident PIDENT] [--query_cover QUERY_COVER]
                  [--subject_cover SUBJECT_COVER] [--evalue EVALUE] [--score SCORE] [--dmnd_algo {auto,0,1,ctg}]
                  [--dmnd_db DMND_DB_FILE]
                  [--sensmode {default,fast,mid-sensitive,sensitive,more-sensitive,very-sensitive,ultra-sensitive}]
                  [--dmnd_iterate {yes,no}] [--matrix {BLOSUM62,BLOSUM90,BLOSUM80,BLOSUM50,BLOSUM45,PAM250,PAM70,PAM30}]
                  [--dmnd_frameshift DMND_FRAMESHIFT] [--gapopen GAPOPEN] [--gapextend GAPEXTEND] [--block_size BLOCK_SIZE]
                  [--index_chunks CHUNKS] [--outfmt_short] [--dmnd_ignore_warnings] [--mmseqs_db MMSEQS_DB_FILE]
                  [--start_sens START_SENS] [--sens_steps SENS_STEPS] [--final_sens FINAL_SENS] [--mmseqs_sub_mat SUBS_MATRIX]
                  [-d HMMER_DB_PREFIX] [--servers_list FILE] [--qtype {hmm,seq}] [--dbtype {hmmdb,seqdb}] [--usemem] [-p PORT]
                  [--end_port PORT] [--num_servers NUM_SERVERS] [--num_workers NUM_WORKERS]
                  [--timeout_load_server TIMEOUT_LOAD_SERVER] [--hmm_maxhits MAXHITS] [--report_no_hits]
                  [--hmm_maxseqlen MAXSEQLEN] [--Z DB_SIZE] [--cut_ga]
                  [--clean_overlaps none|all|clans|hmmsearch_all|hmmsearch_clans] [--no_annot] [--dbmem]
                  [--seed_ortholog_evalue MIN_E-VALUE] [--seed_ortholog_score MIN_SCORE] [--tax_scope TAX_SCOPE]
                  [--tax_scope_mode TAX_SCOPE_MODE] [--target_orthologs {one2one,many2one,one2many,many2many,all}]
                  [--target_taxa LIST_OF_TAX_IDS] [--excluded_taxa LIST_OF_TAX_IDS] [--report_orthologs]
                  [--go_evidence {experimental,non-electronic,all}] [--pfam_realign {none,realign,denovo}] [--md5]
                  [--output FILE_PREFIX] [--output_dir DIR] [--scratch_dir DIR] [--temp_dir DIR] [--no_file_comments]
                  [--decorate_gff DECORATE_GFF] [--decorate_gff_ID_field DECORATE_GFF_ID_FIELD] [--excel]

options:
  -h, --help            show this help message and exit
  -v, --version         show version and exit. (default: False)
  --list_taxa           List taxa available for --tax_scope/--tax_scope_mode, and exit (default: False)

Execution Options:
  --cpu NUM_CPU         Number of CPUs to be used. --cpu 0 to run with all available CPUs. (default: 1)
  --mp_start_method {fork,spawn,forkserver}
                        Sets the python multiprocessing start method. Check
                        https://docs.python.org/3/library/multiprocessing.html. Only use if the default method is not working
                        properly in your OS. (default: spawn)
  --resume              Resumes a previous emapper run, skipping results in existing output files. (default: False)
  --override            Overwrites output files if they exist. By default, execution is aborted if conflicting files are
                        detected. (default: False)

Input Data Options:
  -i FASTA_FILE         Input FASTA file containing query sequences (proteins by default; see --itype and --translate).
                        Required unless -m no_search. (default: None)
  --itype {CDS,proteins,genome,metagenome}
                        Type of data in the input (-i) file. (default: proteins)
  --translate           When --itype CDS, translate CDS to proteins before search. When --itype genome/metagenome and
                        --genepred search, translate predicted CDS from blastx hits to proteins. (default: False)
  --annotate_hits_table SEED_ORTHOLOGS_FILE
                        Annotate TSV formatted table with 4 fields: query, hit, evalue, score. Usually, a .seed_orthologs file
                        from a previous emapper.py run. Requires -m no_search. (default: None)
  -c FILE, --cache FILE
                        File containing annotations and md5 hashes of queries, to be used as cache. Required if -m cache
                        (default: None)
  --data_dir DIR        Path to eggnog-mapper databases. By default, "data/" or the path specified in the environment variable
                        EGGNOG_DATA_DIR. (default: None)

Gene Prediction Options:
  --genepred {search,prodigal}
                        This is applied when --itype genome or --itype metagenome. search: gene prediction is inferred from
                        Diamond/MMseqs2 blastx hits. prodigal: gene prediction is performed using Prodigal. (default: search)
  --trans_table TRANS_TABLE_CODE
                        This option will be used for Prodigal, Diamond or MMseqs2, depending on the mode. For Diamond searches,
                        this option corresponds to the --query-gencode option. For MMseqs2 searches, this option corresponds to
                        the --translation-table option. For Prodigal, this option corresponds to the -g/--trans_table option.
                        It is also used when --translate, check
                        https://biopython.org/docs/1.75/api/Bio.Seq.html#Bio.Seq.Seq.translate. Default is the corresponding
                        programs defaults. (default: None)
  --training_genome FILE
                        A genome to run Prodigal in training mode. If this parameter is used, Prodigal will run in two steps:
                        firstly in training mode, and secondly using the training to analize the emapper input data. See
                        Prodigal documentation about Traning mode for more info. Only used if --genepred prodigal. (default:
                        None)
  --training_file FILE  A training file from Prodigal training mode. If this parameter is used, Prodigal will run using this
                        training file to analyze the emapper input data. Only used if --genepred prodigal. (default: None)
  --allow_overlaps {none,strand,diff_frame,all}
                        When using 'blastx'-based genepred (--genepred search --itype genome/metagenome) this option controls
                        whether overlapping hits are reported or not, or if only those overlapping hits in a different strand
                        or frame are reported. Also, the degree of accepted overlap can be controlled with --overlap_tol.
                        (default: none)
  --overlap_tol FLOAT   This value (0-1) is the proportion such that if (overlap size / hit length) > overlap_tol, hits are
                        considered to overlap. e.g: if overlap_tol is 0.0, any overlap is considered as such. e.g: if
                        overlap_tol is 1.0, one of the hits must overlap entirely to consider that hits do overlap. (default:
                        0.0)

Search Options:
  -m {diamond,mmseqs,hmmer,no_search,cache,novel_fams}
                        diamond: search seed orthologs using diamond (-i is required). mmseqs: search seed orthologs using
                        MMseqs2 (-i is required). hmmer: search seed orthologs using HMMER. (-i is required). no_search: skip
                        seed orthologs search (--annotate_hits_table is required, unless --no_annot). cache: skip seed
                        orthologs search and annotate based on cached results (-i and -c are required).novel_fams: search
                        against the novel families database (-i is required). (default: diamond)

Search filtering common options:
  --pident PIDENT       Report only alignments above or equal to the given percentage of identity (0-100).No effect if -m
                        hmmer. (default: None)
  --query_cover QUERY_COVER
                        Report only alignments above or equal the given percentage of query cover (0-100).No effect if -m
                        hmmer. (default: None)
  --subject_cover SUBJECT_COVER
                        Report only alignments above or equal the given percentage of subject cover (0-100).No effect if -m
                        hmmer. (default: None)
  --evalue EVALUE       Report only alignments below or equal the e-value threshold. (default: 0.001)
  --score SCORE         Report only alignments above or equal the score threshold. (default: None)

Diamond Search Options:
  --dmnd_algo {auto,0,1,ctg}
                        Diamonds --algo option, which can be tuned to search small query sets. By default, it is adjusted
                        automatically. However, the ctg option should be activated manually. If you plan to search a small
                        input set of sequences, use --dmnd_algo ctg to make it faster. (default: auto)
  --dmnd_db DMND_DB_FILE
                        Path to DIAMOND-compatible database (default: None)
  --sensmode {default,fast,mid-sensitive,sensitive,more-sensitive,very-sensitive,ultra-sensitive}
                        Diamonds sensitivity mode. Note that emappers default is sensitive, which is different from diamonds
                        default, which can be activated with --sensmode default. (default: sensitive)
  --dmnd_iterate {yes,no}
                        --dmnd_iterate yes --> activates the --iterate option of diamond for iterative searches, from faster,
                        less sensitive modes, up to the sensitivity specified with --sensmode. Available since diamond 2.0.11.
                        --dmnd_iterate no --> disables the --iterate mode. (default: yes)
  --matrix {BLOSUM62,BLOSUM90,BLOSUM80,BLOSUM50,BLOSUM45,PAM250,PAM70,PAM30}
                        Scoring matrix (default: None)
  --dmnd_frameshift DMND_FRAMESHIFT
                        Diamond --frameshift/-F option. Not used by default. Recommended by diamond: 15. (default: None)
  --gapopen GAPOPEN     Gap open penalty (default: None)
  --gapextend GAPEXTEND
                        Gap extend penalty (default: None)
  --block_size BLOCK_SIZE
                        Diamond -b/--block-size option. Default is the diamonds default. (default: None)
  --index_chunks CHUNKS
                        Diamond -c/--index-chunks option. Default is the diamonds default. (default: None)
  --outfmt_short        Diamond output will include only qseqid sseqid evalue and score. This could help obtain better
                        performance, if also no --pident, --query_cover or --subject_cover thresholds are used. This option is
                        ignored when the diamond search is run in blastx mode for gene prediction (see --genepred). (default:
                        False)
  --dmnd_ignore_warnings
                        Diamond --ignore-warnings option. It avoids Diamond stopping due to warnings (e.g. when a protein
                        contains only ATGC symbols. (default: False)

MMseqs2 Search Options:
  --mmseqs_db MMSEQS_DB_FILE
                        Path to MMseqs2-compatible database (default: None)
  --start_sens START_SENS
                        Starting sensitivity. (default: 3)
  --sens_steps SENS_STEPS
                        Number of sensitivity steps. (default: 3)
  --final_sens FINAL_SENS
                        Final sensititivy step. (default: 7)
  --mmseqs_sub_mat SUBS_MATRIX
                        Matrix to be used for --sub-mat MMseqs2 search option. Default=default used by MMseqs2 (default: None)

HMMER Search Options:
  -d HMMER_DB_PREFIX, --database HMMER_DB_PREFIX
                        specify the target database for sequence searches. Choose among: euk,bact,arch, or a database loaded in
                        a server, db.hmm:host:port (see hmm_server.py) (default: None)
  --servers_list FILE   A FILE with a list of remote hmmpgmd servers. Each row in the file represents a server, in the format
                        'host:port'. If --servers_list is specified, host and port from -d option will be ignored. (default:
                        None)
  --qtype {hmm,seq}     Type of input data (-i). (default: seq)
  --dbtype {hmmdb,seqdb}
                        Type of data in DB (-db). (default: hmmdb)
  --usemem              Use this option to allocate the whole database (-d) in memory using hmmpgmd. If --dbtype hmm, the
                        database must be a hmmpress-ed database. If --dbtype seqdb, the database must be a HMMER-format
                        database created with esl-reformat. Database will be unloaded after execution. Note that this only
                        works for HMMER based searches. To load the eggnog-mapper annotation DB into memory use --dbmem.
                        (default: False)
  -p PORT, --port PORT  Port used to setup HMM server, when --usemem. Also used for --pfam_realign modes. (default: 51700)
  --end_port PORT       Last port to be used to setup HMM server, when --usemem. Also used for --pfam_realign modes. (default:
                        53200)
  --num_servers NUM_SERVERS
                        When using --usemem, specify the number of servers to fire up.Note that cpus specified with --cpu will
                        be distributed among servers and workers. Also used for --pfam_realign modes. (default: 1)
  --num_workers NUM_WORKERS
                        When using --usemem, specify the number of workers per server (--num_servers) to fire up. By default,
                        cpus specified with --cpu will be distributed among servers and workers. Also used for --pfam_realign
                        modes. (default: 1)
  --timeout_load_server TIMEOUT_LOAD_SERVER
                        Number of attempts to load a server on a specific port. If failed, the next numerical port will be
                        tried. (default: 10)
  --hmm_maxhits MAXHITS
                        Max number of hits to report (0 to report all). (default: 1)
  --report_no_hits      Whether queries without hits should be included in the output table. (default: False)
  --hmm_maxseqlen MAXSEQLEN
                        Ignore query sequences larger than `maxseqlen`. (default: 5000)
  --Z DB_SIZE           Fixed database size used in phmmer/hmmscan (allows comparing e-values among databases). (default:
                        40000000)
  --cut_ga              Adds the --cut_ga to hmmer commands (useful for Pfam mappings, for example). See hmmer documentation.
                        (default: False)
  --clean_overlaps none|all|clans|hmmsearch_all|hmmsearch_clans
                        Removes those hits which overlap, keeping only the one with best evalue. Use the "all" and "clans"
                        options when performing a hmmscan type search (i.e. domains are in the database). Use the
                        "hmmsearch_all" and "hmmsearch_clans" options when using a hmmsearch type search (i.e. domains are the
                        queries from -i file). The "clans" and "hmmsearch_clans" and options will only have effect for hits
                        to/from Pfam. (default: None)

Annotation Options:
  --no_annot            Skip functional annotation, reporting only hits. (default: False)
  --dbmem               Use this option to allocate the whole eggnog.db DB in memory. Database will be unloaded after
                        execution. (default: False)
  --seed_ortholog_evalue MIN_E-VALUE
                        Min E-value expected when searching for seed eggNOG ortholog. Queries not having a significant seed
                        orthologs will not be annotated. (default: 0.001)
  --seed_ortholog_score MIN_SCORE
                        Min bit score expected when searching for seed eggNOG ortholog. Queries not having a significant seed
                        orthologs will not be annotated. (default: None)
  --tax_scope TAX_SCOPE
                        Fix the taxonomic scope used for annotation, so only speciation events from a particular clade are used
                        for functional transfer. More specifically, the --tax_scope list is intersected with the seed orthologs
                        clades, and the resulting clades are used for annotation based on --tax_scope_mode. Note that those
                        seed orthologs without clades intersecting with --tax_scope will be filtered out, and wont annotated.
                        Possible arguments for --tax_scope are: 1) A path to a file defined by the user, which contains a list
                        of tax IDs and/or tax names. 2) The name of a pre-configured tax scope, whose source is a file stored
                        within the 'eggnogmapper/annotation/tax_scopes/' directory By default, available ones are: 'auto'
                        ('all'), 'auto_broad' ('all_broad'), 'all_narrow', 'archaea', 'bacteria', 'bacteria_broad',
                        'eukaryota', 'eukaryota_broad' and 'prokaryota_broad'.3) A comma-separated list of taxonomic names
                        and/or taxonomic IDs, sorted by preference. An example of list of tax IDs would be 2759,2157,2,1 for
                        Eukaryota, Archaea, Bacteria and root, in that order of preference. 4) 'none': do not filter out
                        annotations based on taxonomic scope. (default: auto)
  --tax_scope_mode TAX_SCOPE_MODE
                        For a seed ortholog which passed the filter imposed by --tax_scope, --tax_scope_mode controls which
                        specific clade, to which the seed ortholog belongs, will be used for annotation. Options: 1) broadest:
                        use the broadest clade. 2) inner_broadest: use the broadest clade from the intersection with
                        --tax_scope. 3) inner_narrowest: use the narrowest clade from the intersection with --tax_scope. 4)
                        narrowest: use the narrowest clade. 5) A taxonomic scope as in --tax_scope: use this second list to
                        intersect with seed ortholog clades and use the narrowest (as in inner_narrowest) from the intersection
                        to annotate. (default: inner_narrowest)
  --target_orthologs {one2one,many2one,one2many,many2many,all}
                        defines what type of orthologs (in relation to the seed ortholog) should be used for functional
                        transfer (default: all)
  --target_taxa LIST_OF_TAX_IDS
                        Only orthologs from the specified comma-separated list of taxa and all its descendants will be used for
                        annotation transference. By default, all taxa are used. (default: None)
  --excluded_taxa LIST_OF_TAX_IDS
                        Orthologs from the specified comma-separated list of taxa and all its descendants will not be used for
                        annotation transference. By default, no taxa is excluded. (default: None)
  --report_orthologs    Output the list of orthologs found for each query to a .orthologs file (default: False)
  --go_evidence {experimental,non-electronic,all}
                        Defines what type of GO terms should be used for annotation. experimental = Use only terms inferred
                        from experimental evidence. non-electronic = Use only non-electronically curated terms (default: non-
                        electronic)
  --pfam_realign {none,realign,denovo}
                        Realign the queries to the PFAM domains. none = no realignment is performed. PFAM annotation will be
                        that transferred as specify in the --pfam_transfer option. realign = queries will be realigned to the
                        PFAM domains found according to the --pfam_transfer option. denovo = queries will be realigned to the
                        whole PFAM database, ignoring the --pfam_transfer option. Check hmmer options (--num_servers,
                        --num_workers, --port, --end_port) to change how the hmmpgmd server is run. (default: none)
  --md5                 Adds the md5 hash of each query as an additional field in annotations output files. (default: False)

Output options:
  --output FILE_PREFIX, -o FILE_PREFIX
                        base name for output files (default: None)
  --output_dir DIR      Where output files should be written (default: /home/dorian.rojas)
  --scratch_dir DIR     Write output files in a temporary scratch dir, move them to the final output dir when finished. Speed
                        up large computations using network file systems. (default: None)
  --temp_dir DIR        Where temporary files are created. Better if this is a local disk. (default: /home/dorian.rojas)
  --no_file_comments    No header lines nor stats are included in the output files (default: False)
  --decorate_gff DECORATE_GFF
                        Add search hits and/or annotation results to GFF file from gene prediction of a user specified one. no
                        = no GFF decoration at all. GFF file from blastx-based gene prediction will be created anyway. yes =
                        add search hits and/or annotations to GFF file from Prodigal or blastx-based gene prediction. FILE =
                        decorate the specified pre-existing GFF FILE. e.g. --decorage_gff myfile.gff You change the field
                        interpreted as ID of the feature with --decorate_gff_ID_field. (default: no)
  --decorate_gff_ID_field DECORATE_GFF_ID_FIELD
                        Change the field used in GFF files as ID of the feature. (default: ID)
  --excel               Output annotations also in .xlsx format. (default: False)
```

eggNOG-mapper has default settings for accepting protein files. Therefore, the input file type must be set to `--itype metagenome` and `--genepred prodigal` for performing the gene prediction using Progial. Indicate the code to run in all available cpus. The `-o` and `-i` are required and it is recommended to indicate both the `--output_dir` and `--temp_dir` to the same `13-eggnog` folder.

> If the job takes a while running, the results in the common repository `/home/public/met-workshop`

The output from eggNOG-mapper are sveral tab-separated files that present the used orthologs for inference (`workshop.emapper.seed_orthologs`), orthologs hits (`workshop.emapper.hits`), protein prediction with prodigal (`workshop.emapper.genepred.gff` and `workshop.emapper.genepred.fasta`), and the summary table of the metabolic functions annotated `workshop.emapper.annotations`. This last table is the one of interest for this workshop.

```bash
(base) [dorian.rojas@accessnode test]$ ls 13-eggnog/
workshop.emapper.annotations     workshop.emapper.genepred.gff  workshop.emapper.seed_orthologs
workshop.emapper.genepred.fasta  workshop.emapper.hits
```

```bash
(base) [dorian.rojas@accessnode 13-eggnog]$ head workshop.emapper.annotations
## Mon Feb 24 19:04:14 2025
## emapper-2.1.12
## /usr/local/bin/emapper.py --cpu 64 -i 10-final-MAGs/bins/pangenome.fasta --itype metagenome --genepred prodigal -o workshop --output_dir 13-eggnog/ --temp_dir 13-eggnog/
##
#query  seed_ortholog   evalue  score   eggNOG_OGs      max_annot_lvl   COG_category    Description     Preferred_name  GOs     EC       KEGG_ko KEGG_Pathway    KEGG_Module     KEGG_Reaction   KEGG_rclass     BRITE   KEGG_TC CAZy    BiGG_Reaction   PFAMs
SRR8555091_bin.1.permissive_0_1 906968.Trebr_2272       3.66e-36        124.0   COG0100@1|root,COG0100@2|Bacteria,2J7BP@203691|Spirochaetes      203691|Spirochaetes     J       Located on the platform of the 30S subunit, it bridges several disparate RNA helices of the 16S rRNA. Forms part of the Shine-Dalgarno cleft in the 70S ribosome rpsK    -       -       ko:K02948       ko03010,map03010 M00178,M00179   -       -       br01610,ko00000,ko00001,ko00002,ko03011 -       -       -       Ribosomal_S11
SRR8555091_bin.1.permissive_0_2 1124982.MSI_07440       2.34e-66        203.0   COG0099@1|root,COG0099@2|Bacteria,2J7QU@203691|Spirochaetes      203691|Spirochaetes     J       Located at the top of the head of the 30S subunit, it contacts several helices of the 16S rRNA. In the 70S ribosome it contacts the 23S rRNA (bridge B1a) and protein L5 of the 50S subunit (bridge B1b), connecting the 2 subunits       rpsM    -       -       ko:K02952       ko03010,map03010        M00178,M00179   -       -       br01610,ko00000,ko00001,ko00002,ko03011  -       -       -       Ribosomal_S13
SRR8555091_bin.1.permissive_0_3 1124982.MSI_07420       4.86e-253       702.0   COG0201@1|root,COG0201@2|Bacteria,2J5BH@203691|Spirochaetes      203691|Spirochaetes     U       The central subunit of the protein translocation channel SecYEG. Consists of two halves formed by TMs 1-5 and 6-10. These two domains form a lateral gate at the front which open onto the bilayer between TMs 2 and 7, and are clamped together by SecE at the back. The channel is closed by both a pore ring composed of hydrophobic SecY resides and a short helix (helix 2A) on the extracellular side of the membrane which forms a plug. The plug probably moves laterally to allow the channel to open. The ring and the pore may move independently        secY    -       -       ko:K03076       ko02024,ko03060,ko03070,map02024,map03060,map03070       M00335  -       -       ko00000,ko00001,ko00002,ko02044 3.A.5   -       -       SecY
SRR8555091_bin.1.permissive_0_4 877418.ATWV01000004_gene1907    5.16e-81        243.0   COG0200@1|root,COG0200@2|Bacteria,2J7ZW@203691|Spirochaetes      203691|Spirochaetes     J       Binds to the 23S rRNA   rplO    -       -       ko:K02876       ko03010,map03010 M00178,M00179   -       -       br01610,ko00000,ko00001,ko00002,ko03011 -       -       -       Ribosomal_L27A
SRR8555091_bin.1.permissive_0_5 1124982.MSI_07390       9.72e-89        263.0   COG0098@1|root,COG0098@2|Bacteria,2J60F@203691|Spirochaetes      203691|Spirochaetes     J       Located at the back of the 30S subunit body where it stabilizes the conformation of the head with respect to the body    rpsE    -       -       ko:K02988       ko03010,map03010        M00178,M00179   -       -br01610,ko00000,ko00001,ko00002,ko03011 -       -       -       Ribosomal_S5,Ribosomal_S5_C
```

The `KEGG_Pathway` columns indicate IDs for metabolic pathway of the KEGG database. In this regard, the `Description` and `PFAMs` present an overall idea of the function and gene annotated in the pangenome.

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


# Extract list of genomes (MAGs)
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

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="abricate"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL

cd /home/dorian.rojas/test

CTN_PATH=/opt/ohpc/pub/containers/BIO/


mkdir -p 12-abricate/

for db in {argannot,card,ecoh,ecoli_vf,megares,ncbi,plasmidfinder,resfinder,vfdb}; do

$CTN_PATH/abricate-1.0.1.sif abricate --threads 64 --db $db 10-final_MAGs/bins/pangenome.fasta \
        >> 12-abricate/${db}.tab

done

$CTN_PATH/abricate-1.0.1.sif abricate --sumary 12-abricate/*.tab \
        >> 12-abricate/summary.tab

date
time
```

**`eggnog-mapper.slurm` file:**

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="eggnog"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL

cd /home/dorian.rojas/test

CTN_PATH=/opt/ohpc/pub/containers/BIO/


mkdir -p 13-eggnog/

export EGGNOG_DATA_DIR=/home/dorian.rojas/DB/eggNOG_DB/

$CTN_PATH/eggnog-mapper-2.1.12.sif emapper.py --cpu 64 \
        -i 10-final-MAGs/bins/pangenome.fasta \
        --itype metagenome --genepred prodigal \
        -o workshop --output_dir 13-eggnog/ --temp_dir 13-eggnog/


date
time
```

```vim
user_genome     classification  closest_genome_reference        closest_genome_reference_radius closest_genome_taxonomy closest_genome_ani      closest_genome_af       closest_placement_reference closest_placement_radius        closest_placement_taxonomy      closest_placement_ani   closest_placement_af    pplacer_taxonomy        classification_method   noteother_related_references(genome_id,species_name,radius,ANI,AF)  msa_percent     translation_table       red_value       warnings
49647_1#2_clean_bin.10.permissive       d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__Alistipes;s__Alistipes finegoldii       GCF_000265365.1 95.0d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__Alistipes;s__Alistipes finegoldii       98.78   0.925   GCF_000265365.1 95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__Alistipes;s__Alistipes finegoldii   98.78   0.925   d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__Alistipes;s__       taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments        GCF_021204515.1, s__Alistipes sp021204515, 95.0, 93.47, 0.696; GCA_934754975.1, s__Alistipes sp934754975, 95.0, 91.79, 0.566; GCF_025145845.1, s__Alistipes shahii, 95.0, 86.39, 0.225; GCF_025145285.1, s__Alistipes onderdonkii, 95.0, 85.64, 0.287; GCF_902388705.1, s__Alistipes sp902388705, 95.0, 84.9, 0.232; GCF_025145645.1, s__Alistipes senegalensis, 95.0, 84.13, 0.169; GCF_900541585.1, s__Alistipes sp900541585, 95.0, 84.05, 0.17; GCF_900083545.1, s__Alistipes provencensis, 95.0, 83.76, 0.154   56.48   11      N/A     N/A
49647_1#2_clean_bin.11.orig     d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Sutterella;s__Sutterella wadsworthensis_A        GCF_000297775.1     95.0    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Sutterella;s__Sutterella wadsworthensis_A        99.23   0.939       GCF_000297775.1 95.0    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Sutterella;s__Sutterella wadsworthensis_A        99.23       0.939   d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Sutterella;s__   taxonomic classification defined by topology and ANItopological placement and ANI have congruent species assignments        GCA_902493085.1, s__Sutterella sp902493085, 95.0, 88.6, 0.457; GCA_900764215.1, s__Sutterella sp900764215, 95.0, 87.72, 0.527; GCA_934249275.1, s__Sutterella sp934249275, 95.0, 86.22, 0.428       98.37   11      N/A     N/A
49647_1#2_clean_bin.15.strict   d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia;s__Roseburia inulinivorans   GCF_000174195.1 95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia;s__Roseburia inulinivorans       97.6    0.795   GCF_000174195.1 95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia;s__Roseburia inulinivorans       97.6    0.795   d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia;s__      taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments        GCA_019411365.1, s__Roseburia sp019411365, 95.0, 95.66, 0.28; GCA_900542495.1, s__Roseburia sp900542495, 95.0, 90.78, 0.535 40.36   11      N/A     N/A
49647_1#2_clean_bin.16.orig     d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Duodenibacillus;s__Duodenibacillus sp900544255   GCA_900544255.1     95.0    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Duodenibacillus;s__Duodenibacillus sp900544255   97.75   0.849       GCA_900544255.1 95.0    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Duodenibacillus;s__Duodenibacillus sp900544255   97.75       0.849   d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Duodenibacillus;s__      taxonomic classification defined by topology and ANI        topological placement and ANI have congruent species assignments        N/A     69.06   11      N/A     N/A
49647_1#2_clean_bin.19.orig     d__Bacteria;p__Actinomycetota;c__Coriobacteriia;o__Coriobacteriales;f__Coriobacteriaceae;g__Collinsella;s__     N/A     N/A     N/A     N/A     N/AGCA_030826265.1  95.0    d__Bacteria;p__Actinomycetota;c__Coriobacteriia;o__Coriobacteriales;f__Coriobacteriaceae;g__Collinsella;s__Collinsella aerofaciens_AY   93.98   0.835   d__Bacteria;p__Actinomycetota;c__Coriobacteriia;o__Coriobacteriales;f__Coriobacteriaceae;g__Collinsella;s__ taxonomic classification defined by topology and ANI    classification based on placement in class-level tree       GCA_902494085.1, s__Collinsella sp902494085, 95.0, 94.48, 0.769; GCA_902467875.1, s__Collinsella sp902467875, 95.0, 94.39, 0.827; GCA_000434535.1, s__Collinsella sp000434535, 95.0, 94.38, 0.777; GCA_902470895.1, s__Collinsella sp902470895, 95.0, 94.33, 0.799; GCA_958438955.1, s__Collinsella sp958438955, 95.0, 94.28, 0.777; GCA_900543515.1, s__Collinsella sp900543515, 95.0, 94.26, 0.796; GCF_015670795.1, s__Collinsella aerofaciens_Q, 95.0, 94.25, 0.84; GCA_902363125.1, s__Collinsella sp902363125, 95.0, 94.25, 0.832; GCA_018365375.1, s__Collinsella sp018365375, 95.0, 94.25, 0.821; GCA_902494195.1, s__Collinsella sp902494195, 95.0, 94.25, 0.812; GCA_902462065.1, s__Collinsella sp902462065, 95.0, 94.25, 0.81; GCA_902489665.1, s__Collinsella sp902489665, 95.0, 94.25, 0.808; GCA_900760385.1, s__Collinsella sp900760385, 95.0, 94.25, 0.786; GCA_900547345.1, s__Collinsella sp900547345, 95.0, 94.25, 0.776; GCA_014871645.1, s__Collinsella sp014871645, 95.0, 94.24, 0.82; GCF_001405375.1, s__Collinsella aerofaciens_F, 95.0, 94.24, 0.817; GCF_003458415.1, s__Collinsella sp003458415, 95.0, 94.24, 0.805; GCA_958427355.1, s__Collinsella sp958427355, 95.0, 94.24, 0.804; GCA_902495985.1, s__Collinsella sp902495985, 95.0, 94.22, 0.747; GCA_900545835.1, s__Collinsella sp900545835, 95.0, 94.21, 0.778; GCA_902495945.1, s__Collinsella sp902495945, 95.0, 94.18, 0.776; GCA_900545505.1, s__Collinsella sp900545505, 95.0, 94.17, 0.76; GCF_019041915.1, s__Collinsella sp019041915, 95.0, 94.15, 0.833; GCA_902496335.1, s__Collinsella sp902496335, 95.0, 94.15, 0.797; GCA_902470805.1, s__Collinsella sp902470805, 95.0, 94.12, 0.83; GCF_003438495.1, s__Collinsella sp003438495, 95.0, 94.12, 0.822; GCA_902494535.1, s__Collinsella sp902494535, 95.0, 94.12, 0.817; GCA_958407465.1, s__Collinsella sp958407465, 95.0, 94.12, 0.813; GCA_900547765.1, s__Collinsella sp900547765, 95.0, 94.12, 0.802; GCA_902493375.1, s__Collinsella sp902493375, 95.0, 94.12, 0.798; GCA_902497065.1, s__Collinsella sp902497065, 95.0, 94.12, 0.793; GCF_015556945.1, s__Collinsella aerofaciens_AO, 95.0, 94.11, 0.829; GCA_018376265.1, s__Collinsella sp018376265, 95.0, 94.11, 0.826; GCA_958413895.1, s__Collinsella sp958413895, 95.0, 94.11, 0.818; GCA_900544425.1, s__Collinsella sp900544425, 95.0, 94.07, 0.756; GCA_958346165.1, s__Collinsella sp958346165, 95.0, 94.07, 0.687; GCA_902489515.1, s__Collinsella sp902489515, 95.0, 94.03, 0.734; GCA_937937835.1, s__Collinsella sp937937835, 95.0, 94.02, 0.81; GCA_902489965.1, s__Collinsella sp902489965, 95.0, 94.02, 0.802; GCA_900540095.1, s__Collinsella sp900540095, 95.0, 94.02, 0.795; GCA_902493125.1, s__Collinsella sp902493125, 95.0, 94.02, 0.782; GCA_019413955.1, s__Collinsella sp019413955, 95.0, 94.02, 0.77; GCA_902478865.1, s__Collinsella sp902478865, 95.0, 94.02, 0.746; GCA_902478965.1, s__Collinsella sp902478965, 95.0, 94.01, 0.803; GCA_902493985.1, s__Collinsella sp902493985, 95.0, 94.01, 0.799; GCA_900541185.1, s__Collinsella sp900541185, 95.0, 94.01, 0.793; GCA_958421655.1, s__Collinsella sp958421655, 95.0, 94.01, 0.791; GCA_902478825.1, s__Collinsella sp902478825, 95.0, 94.01, 0.76; GCA_938023085.1, s__Collinsella sp938023085, 95.0, 94.0, 0.789; GCF_015558785.1, s__Collinsella aerofaciens_P, 95.07, 93.99, 0.836; GCF_028206995.1, s__Collinsella aerofaciens_AQ, 95.0, 93.99, 0.823; GCA_958429345.1, s__Collinsella sp958429345, 95.0, 93.99, 0.805; GCA_014871655.1, s__Collinsella sp014871655, 95.0, 93.99, 0.8; GCA_902488525.1, s__Collinsella sp902488525, 95.0, 93.99, 0.799; GCA_958407995.1, s__Collinsella sp958407995, 95.0, 93.99, 0.798; GCA_900541065.1, s__Collinsella sp900541065, 95.0, 93.98, 0.837; GCA_018381235.1, s__Collinsella sp018381235, 95.02, 93.98, 0.82; GCA_900541135.1, s__Collinsella sp900541135, 95.0, 93.98, 0.815; GCF_010509075.1, s__Collinsella aerofaciens, 95.0, 93.98, 0.813; GCF_025289055.1, s__Collinsella sp900549025, 95.0, 93.98, 0.813; GCF_018785245.1, s__Collinsella aerofaciens_AD, 95.0, 93.98, 0.81; GCA_902492235.1, s__Collinsella sp902492235, 95.0, 93.98, 0.806; GCA_958352845.1, s__Collinsella sp958352845, 95.0, 93.98, 0.804; GCF_015670745.1, s__Collinsella aerofaciens_T, 95.0, 93.98, 0.803; GCA_900541175.1, s__Collinsella sp900541175, 95.0, 93.98, 0.803; GCA_958405195.1, s__Collinsella sp958405195, 95.0, 93.98, 0.803; GCA_900763515.1, s__Collinsella sp900763515, 95.0, 93.98, 0.795; GCA_902491055.1, s__Collinsella sp902491055, 95.0, 93.98, 0.778; GCA_900758885.1, s__Collinsella sp900758885, 95.0, 93.98, 0.763; GCA_902492545.1, s__Collinsella sp902492545, 95.0, 93.95, 0.765; GCF_000763055.1, s__Collinsella sp000763055, 95.0, 93.92, 0.806; GCA_902496005.1, s__Collinsella sp902496005, 95.0, 93.92, 0.774; GCA_900548565.1, s__Collinsella sp900548565, 95.0, 93.92, 0.764; GCA_900760245.1, s__Collinsella sp900760245, 95.0, 93.92, 0.75; GCA_902494605.1, s__Collinsella sp902494605, 95.0, 93.9, 0.78; GCA_902478995.1, s__Collinsella sp902478995, 95.0, 93.88, 0.773; GCA_900760905.1, s__Collinsella sp900760905, 95.0, 93.87, 0.798; GCA_019414245.1, s__Collinsella sp019414245, 95.0, 93.87, 0.79; GCA_900752015.1, s__Collinsella sp900752015, 95.0, 93.87, 0.784; GCA_014871665.1, s__Collinsella sp014871665, 95.0, 93.87, 0.777; GCA_902469555.1, s__Collinsella sp902469555, 95.0, 93.87, 0.767; GCA_900544235.1, s__Collinsella sp900544235, 95.0, 93.85, 0.791; GCA_902490565.1, s__Collinsella sp902490565, 95.0, 93.85, 0.785; GCA_900544995.1, s__Collinsella sp900544995, 95.0, 93.85, 0.777; GCA_958437935.1, s__Collinsella sp958437935, 95.0, 93.84, 0.817; GCA_902469665.1, s__Collinsella sp902469665, 95.0, 93.84, 0.802; GCA_900539735.1, s__Collinsella sp900539735, 95.0, 93.84, 0.77; GCA_900544865.1, s__Collinsella sp900544865, 95.0, 93.83, 0.797; GCA_900757505.1, s__Collinsella sp900757505, 95.0, 93.83, 0.763; GCA_900552755.1, s__Collinsella sp900552755, 95.0, 93.82, 0.695; GCA_902495155.1, s__Collinsella sp902495155, 95.0, 93.78, 0.783; GCA_905197385.1, s__Collinsella sp905197385, 95.0, 93.77, 0.752; GCA_008671285.1, s__Collinsella sp003487125, 95.0, 93.7, 0.693; GCA_019414325.1, s__Collinsella sp019414325, 95.0, 93.69, 0.805; GCA_958427975.1, s__Collinsella sp958427975, 95.0, 93.57, 0.788; GCA_937950415.1, s__Collinsella sp937950415, 95.0, 93.56, 0.803; GCA_900542945.1, s__Collinsella sp900542945, 95.0, 93.32, 0.777; GCA_022713905.1, s__Collinsella sp022713905, 95.0, 93.29, 0.781; GCA_022720335.1, s__Collinsella sp022720335, 95.0, 93.17, 0.788   91.38       11      0.99658 Genome not assigned to closest species as it falls outside its pre-defined ANI radius
49647_1#2_clean_bin.2.strict    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola massiliensis        GCF_000382445.1 95.0d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola massiliensis        98.9    0.891   GCF_000382445.1 95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola massiliensis    98.9    0.891   d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__    taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments        GCA_900760795.1, s__Phocaeicola sp900760795, 95.0, 93.84, 0.54; GCF_021730445.1, s__Phocaeicola faecalis, 95.0, 90.74, 0.749; GCA_948956705.1, s__Phocaeicola sp948956705, 95.0, 87.88, 0.251   64.07   11      N/A     N/A
49647_1#2_clean_bin.20.orig     d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Alitiscatomonas;s__Alitiscatomonas sp900066535 GCF_003435375.1 95.0d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Alitiscatomonas;s__Alitiscatomonas sp900066535 97.76   0.922   GCF_003435375.1 95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Alitiscatomonas;s__Alitiscatomonas sp900066535     97.76   0.922   d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Alitiscatomonas;s__        taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments    GCA_902496665.1, s__Alitiscatomonas sp902496665, 95.0, 94.34, 0.811; GCF_003480435.1, s__Alitiscatomonas sp900066055, 95.0, 93.59, 0.73; GCA_934370915.1, s__Alitiscatomonas sp934370915, 95.0, 91.99, 0.589; GCF_003478035.1, s__Alitiscatomonas sp900066785, 95.0, 89.84, 0.466   67.81   11      N/A     N/A
49647_1#2_clean_bin.21.orig     d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__Agathobacter rectalis  GCF_000020605.1 95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__Agathobacter rectalis      97.84   0.792   GCF_000020605.1 95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__Agathobacter rectalis      97.84   0.792   d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__   taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments        GCA_902786155.1, s__Agathobacter sp900317585, 95.0, 94.78, 0.632; GCA_937935345.1, s__Agathobacter sp900546625, 95.0, 94.66, 0.728; GCA_934365945.1, s__Agathobacter sp934365945, 95.0, 93.15, 0.574    90.07   11      N/A     N/A
49647_1#2_clean_bin.22.orig     d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Hominicoprocola;s__Hominicoprocola fusiformis        GCA_020687125.1     95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Hominicoprocola;s__Hominicoprocola fusiformis        96.72   0.825   GCA_020687125.1     95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Hominicoprocola;s__Hominicoprocola fusiformis        96.72   0.825   d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Hominicoprocola;s__      taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments    GCA_934835795.1, s__Hominicoprocola sp934835795, 95.0, 96.64, 0.767; GCA_937884915.1, s__Hominicoprocola sp937884915, 95.0, 92.77, 0.652; GCA_900552015.1, s__Hominicoprocola sp900552015, 95.0, 90.82, 0.647; GCA_900550165.1, s__Hominicoprocola sp900550165, 95.0, 89.96, 0.634; GCA_900763705.1, s__Hominicoprocola sp900763705, 95.0, 89.92, 0.684; GCA_003522105.1, s__Hominicoprocola sp003522105, 95.0, 88.43, 0.539; GCA_934597535.1, s__Hominicoprocola sp934597535, 95.0, 88.38, 0.505; GCA_900546295.1, s__Hominicoprocola sp900546295, 95.0, 88.19, 0.577; GCA_934636865.1, s__Hominicoprocola sp934636865, 95.0, 87.44, 0.547; GCA_934561515.1, s__Hominicoprocola sp900556145, 95.0, 87.07, 0.495; GCA_934572445.1, s__Hominicoprocola sp934572445, 95.0, 86.81, 0.394; GCA_002405995.1, s__Hominicoprocola sp002405995, 95.0, 86.66, 0.453; GCA_934572465.1, s__Hominicoprocola sp934572465, 95.0, 85.68, 0.431; GCA_900768135.1, s__Hominicoprocola sp900768135, 95.0, 84.67, 0.344        98.17   11      N/A     N/A
49647_1#2_clean_bin.23.strict   d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Mycoplasmatales;f__UBA3375;g__UBA3375;s__UBA3375 sp900551955 GCA_900551955.1 95.0    d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Mycoplasmatales;f__UBA3375;g__UBA3375;s__UBA3375 sp900551955     97.14   0.647   GCA_900551955.1 95.0    d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Mycoplasmatales;f__UBA3375;g__UBA3375;s__UBA3375 sp900551955     97.14   0.647   d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Mycoplasmatales;f__UBA3375;g__UBA3375;s__    taxonomic classification defined by topology and ANI        topological placement and ANI have congruent species assignments        N/A     51.08   11      N/A     N/A
49647_1#2_clean_bin.25.orig     d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-217;s__CAG-217 sp000436335     GCA_025757685.1 95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-217;s__CAG-217 sp000436335 97.52   0.926   GCA_025757685.1 95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-217;s__CAG-217 sp000436335 97.52   0.926   d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-217;s__    taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments        N/A     95.11   11      N/AN/A
49647_1#2_clean_bin.26.orig     d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Mediterraneibacter;s__Mediterraneibacter faecis        GCA_001312505.1     95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Mediterraneibacter;s__Mediterraneibacter faecis        96.94   0.807   GCA_001312505.1     95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Mediterraneibacter;s__Mediterraneibacter faecis        96.94   0.807   d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Mediterraneibacter;s__     taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments    GCF_934309345.1, s__Mediterraneibacter sp934309345, 95.0, 94.1, 0.769   76.66   11      N/A     Genome has more than 10.8% of markers with multiple hits
49647_1#2_clean_bin.28.orig     d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__Parabacteroides distasonis  GCF_000012845.1 95.0d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__Parabacteroides distasonis  97.3    0.781   GCF_000012845.1 95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__Parabacteroides distasonis      97.3    0.781   d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__        taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments    GCF_004793765.1, s__Parabacteroides distasonis_A, 95.0, 94.2, 0.626; GCF_011038785.1, s__Parabacteroides sp011038785, 95.0, 90.43, 0.569; GCA_028689545.1, s__Parabacteroides sp028689545, 95.0, 90.23, 0.597; GCF_947855115.1, s__Parabacteroides sp947855115, 95.0, 87.43, 0.525  92.85   11      N/A     N/A
49647_1#2_clean_bin.29.strict   d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Holdemanella;s__Holdemanella sp900556915        GCF_024460475.1     95.0    d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Holdemanella;s__Holdemanella sp900556915        96.61   0.755   GCF_024460475.1     95.0    d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Holdemanella;s__Holdemanella sp900556915        96.61   0.755   d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Holdemanella;s__    taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments    GCF_014306155.1, s__Holdemanella hominis, 95.0, 95.67, 0.699; GCF_009696075.1, s__Holdemanella porci, 95.95, 95.66, 0.762; GCA_003436425.1, s__Holdemanella sp003436425, 95.0, 94.27, 0.713; GCF_934772795.1, s__Holdemanella sp934772795, 95.0, 93.21, 0.7; GCF_000156655.1, s__Holdemanella biformis, 95.0, 90.91, 0.615; GCA_938036765.1, s__Holdemanella sp900547815, 95.0, 90.79, 0.658; GCF_900754615.1, s__Holdemanella sp900754615, 95.0, 89.64, 0.602; GCF_934754935.1, s__Holdemanella sp934754935, 95.0, 89.44, 0.65; GCF_003458715.1, s__Holdemanella sp003458715, 95.0, 89.33, 0.608   91.2    11      N/A     N/A
49647_1#2_clean_bin.30.orig     d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella sp900544825   GCA_900544825.1 95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella sp900544825       98.02   0.87    GCA_900544825.1 95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella sp900544825       98.02   0.87    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__     taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments        GCF_900551985.1, s__Prevotella sp900551985, 95.0, 94.61, 0.689; GCA_002224675.1, s__Prevotella copri_A, 95.0, 93.89, 0.632; GCA_022741625.1, s__Prevotella sp022741625, 95.0, 93.53, 0.64; GCF_026095235.1, s__Prevotella sp900557035, 95.0, 93.36, 0.634; GCA_021636625.1, s__Prevotella sp021636625, 95.0, 93.12, 0.568; GCA_900765465.1, s__Prevotella sp900765465, 95.0, 91.17, 0.557; GCA_958369755.1, s__Prevotella sp958369755, 95.0, 89.0, 0.43; GCF_019249655.1, s__Prevotella sp900767615, 95.85, 88.82, 0.431; GCA_900556795.1, s__Prevotella sp900556795, 95.0, 88.69, 0.395; GCA_022732315.1, s__Prevotella sp022732315, 95.0, 88.33, 0.326; GCF_009494395.1, s__Prevotella sp900546535, 95.0, 88.03, 0.374; GCF_009494135.1, s__Prevotella copri_B, 96.82, 87.77, 0.306; GCA_003465495.1, s__Prevotella copri_G, 96.82, 87.32, 0.318; GCA_900557255.1, s__Prevotella sp900557255, 95.0, 87.29, 0.227; GCA_949820605.1, s__Prevotella sp934191715, 96.61, 86.98, 0.29; GCA_003463005.1, s__Prevotella copri_I, 96.73, 86.95, 0.305; GCA_934196725.1, s__Prevotella sp934196725, 95.0, 86.94, 0.33; GCA_026015295.1, s__Prevotella copri_C, 95.51, 86.87, 0.297; GCA_900555035.1, s__Prevotella sp900555035, 95.0, 86.86, 0.378; GCA_026015625.1, s__Prevotella copri_H, 96.58, 86.8, 0.306; GCA_015074785.1, s__Prevotella sp015074785, 96.07, 86.75, 0.284; GCA_900770515.1, s__Prevotella sp900770515, 95.0, 86.63, 0.309; GCF_026015465.1, s__Prevotella copri_E, 96.47, 86.58, 0.276; GCA_021622245.1, s__Prevotella copri_J, 95.95, 86.4, 0.3; GCF_009494835.1, s__Prevotella copri_F, 96.41, 86.38, 0.297; GCA_934202525.1, s__Prevotella sp934202525, 95.0, 86.38, 0.294; GCF_019127135.1, s__Prevotella sp900551275, 95.0, 86.35, 0.276; GCF_025151535.1, s__Prevotella copri, 96.58, 85.98, 0.291; GCA_021531135.1, s__Prevotella sp945863825, 95.0, 84.19, 0.184; GCF_004535825.1, s__Prevotella hominis, 95.0, 83.91, 0.162; GCA_900548535.1, s__Prevotella sp900548535, 95.0, 82.29, 0.193; GCA_002297965.1, s__Prevotella sp002297965, 95.0, 81.79, 0.184     96.39   11      N/A     N/A
49647_1#2_clean_bin.7.strict    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium duncaniae        GCF_010509575.1     95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium duncaniae        97.63   0.841   GCF_010509575.1     95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium duncaniae        97.63   0.841   d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__      taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments    GCA_937915335.1, s__Faecalibacterium sp937915335, 95.0, 88.82, 0.495; GCA_934526125.1, s__Faecalibacterium sp934526125, 95.0, 86.7, 0.415; GCF_019967995.1, s__Faecalibacterium prausnitzii_E, 95.0, 86.58, 0.433; GCF_003287455.1, s__Faecalibacterium hattorii, 95.0, 86.31, 0.423; GCF_003287405.1, s__Faecalibacterium prausnitzii_J, 95.0, 85.72, 0.316; GCA_934547065.1, s__Faecalibacterium sp934547065, 95.0, 85.65, 0.36; GCF_003324185.1, s__Faecalibacterium prausnitzii, 95.0, 85.39, 0.324; GCF_023347295.1, s__Faecalibacterium sp003449675, 95.0, 85.21, 0.315; GCF_023347235.1, s__Faecalibacterium prausnitzii_A, 95.0, 85.17, 0.301; GCA_025757845.1, s__Faecalibacterium sp900539945, 95.0, 85.16, 0.328; GCF_023347535.1, s__Faecalibacterium sp900539885, 95.0, 85.16, 0.316; GCF_003287495.1, s__Faecalibacterium prausnitzii_I, 95.0, 85.1, 0.333; GCA_934503275.1, s__Faecalibacterium sp934503275, 95.0, 85.05, 0.294; GCA_934540955.1, s__Faecalibacterium sp934540955, 95.0, 85.02, 0.286; GCA_934524985.1, s__Faecalibacterium sp934524985, 95.0, 84.97, 0.284; GCA_934561555.1, s__Faecalibacterium sp934561555, 95.0, 84.83, 0.289; GCF_023347255.1, s__Faecalibacterium prausnitzii_D, 95.0, 84.82, 0.291; GCA_900772565.1, s__Faecalibacterium sp900772565, 95.0, 84.81, 0.281; GCF_002549775.1, s__Faecalibacterium prausnitzii_F, 95.0, 84.78, 0.297; GCA_934525845.1, s__Faecalibacterium sp934525845, 95.0, 84.71, 0.292; GCA_902463905.1, s__Faecalibacterium sp902463905, 95.0, 84.58, 0.214; GCA_934539245.1, s__Faecalibacterium sp934539245, 95.0, 84.56, 0.286; GCA_934561115.1, s__Faecalibacterium sp934561115, 95.0, 84.44, 0.252; GCA_014858325.1, s__Faecalibacterium sp014858325, 95.0, 84.43, 0.241; GCF_020687245.1, s__Faecalibacterium longum, 95.0, 84.38, 0.278; GCA_934512105.1, s__Faecalibacterium sp934512105, 95.0, 84.34, 0.245; GCA_934513205.1, s__Faecalibacterium sp934513205, 95.0, 84.26, 0.248; GCA_934532365.1, s__Faecalibacterium sp934532365, 95.0, 84.19, 0.251; GCA_900551435.1, s__Faecalibacterium sp900551435, 95.0, 84.05, 0.292; GCA_934514285.1, s__Faecalibacterium sp934514285, 95.0, 84.05, 0.268; GCA_902463235.1, s__Faecalibacterium sp902463235, 95.0, 84.04, 0.159; GCF_023347275.1, s__Faecalibacterium sp900758465, 95.0, 83.94, 0.262; GCA_900765105.1, s__Faecalibacterium sp900765105, 95.0, 83.94, 0.226; GCA_900765705.1, s__Faecalibacterium sp900765705, 95.0, 83.8, 0.215; GCA_934553285.1, s__Faecalibacterium sp934553285, 95.0, 83.78, 0.243; GCA_934561485.1, s__Faecalibacterium sp934561485, 95.0, 83.68, 0.227; GCA_934525325.1, s__Faecalibacterium sp934525325, 95.0, 83.64, 0.215; GCA_905215595.1, s__Faecalibacterium sp905215595, 95.0, 83.61, 0.207; GCA_934525735.1, s__Faecalibacterium sp934525735, 95.0, 83.51, 0.219; GCA_934525405.1, s__Faecalibacterium sp934525405, 95.0, 83.5, 0.211; GCA_934501115.1, s__Faecalibacterium sp934501115, 95.0, 83.38, 0.229; GCA_934511465.1, s__Faecalibacterium sp934511465, 95.0, 83.36, 0.217; GCA_934533135.1, s__Faecalibacterium sp934533135, 95.0, 83.35, 0.231; GCA_004558805.1, s__Faecalibacterium prausnitzii_M, 95.0, 82.91, 0.204; GCA_002313795.1, s__Faecalibacterium prausnitzii_L, 95.0, 82.39, 0.175; GCA_934518555.1, s__Faecalibacterium sp934518555, 95.0, 82.34, 0.185; GCA_934548045.1, s__Faecalibacterium sp934548045, 95.0, 81.81, 0.187     51.08   11      N/A     N/A
49647_1#2_clean_bin.8.permissive        d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Mediterraneibacter;s__Mediterraneibacter torques       GCF_000153925.1     95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Mediterraneibacter;s__Mediterraneibacter torques       99.62   0.949   GCF_000153925.1     95.0    d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Mediterraneibacter;s__Mediterraneibacter torques       99.62   0.949   d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Mediterraneibacter;s__     taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments    GCA_905209865.1, s__Mediterraneibacter sp900752395, 95.0, 86.14, 0.215  41.75   11      N/A     N/A
49647_1#2_clean_bin.9.permissive        d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Succinivibrionaceae;g__Succinivibrio;s__Succinivibrio sp000431835       GCA_000431835.1 95.0    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Succinivibrionaceae;g__Succinivibrio;s__Succinivibrio sp000431835   97.87       0.88    GCA_000431835.1 95.0    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Succinivibrionaceae;g__Succinivibrio;s__Succinivibrio sp000431835       97.87   0.88    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Succinivibrionaceae;g__Succinivibrio;s__    taxonomic classification defined by topology and ANI        topological placement and ANI have congruent species assignments        GCA_002449335.1, s__Succinivibrio sp002449335, 95.0, 90.62, 0.667; GCA_900315395.1, s__Succinivibrio sp900315395, 95.0, 90.61, 0.624; GCA_934218695.1, s__Succinivibrio sp934218695, 95.0, 87.07, 0.322     81.17   11      N/A     N/A
49647_1#70_clean_bin.1.strict   d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola sp000434735 GCA_000434735.1 95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola sp000434735     97.22   0.829   GCA_000434735.1 95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola sp000434735     97.22   0.829   d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__    taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments        N/A     79.611      N/A     N/A
49647_1#70_clean_bin.2.strict   d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri_A       GCA_002224675.1 95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri_A   96.08   0.76    GCA_002224675.1 95.0    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri_A   96.08   0.76    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__     taxonomic classification defined by topology and ANI    topological placement and ANI have congruent species assignments        GCA_021636625.1, s__Prevotella sp021636625, 95.0, 95.39, 0.687; GCA_022741625.1, s__Prevotella sp022741625, 95.0, 94.81, 0.673; GCF_026095235.1, s__Prevotella sp900557035, 95.0, 94.78, 0.716; GCA_900544825.1, s__Prevotella sp900544825, 95.0, 94.36, 0.709; GCF_900551985.1, s__Prevotella sp900551985, 95.0, 93.14, 0.676; GCA_900765465.1, s__Prevotella sp900765465, 95.0, 92.45, 0.522; GCA_958369755.1, s__Prevotella sp958369755, 95.0, 89.23, 0.309; GCA_900556795.1, s__Prevotella sp900556795, 95.0, 88.5, 0.321; GCF_019249655.1, s__Prevotella sp900767615, 95.85, 87.69, 0.463; GCA_022732315.1, s__Prevotella sp022732315, 95.0, 86.88, 0.345; GCF_009494395.1, s__Prevotella sp900546535, 95.0, 86.79, 0.36; GCA_900557255.1, s__Prevotella sp900557255, 95.0, 86.59, 0.157; GCA_900770515.1, s__Prevotella sp900770515, 95.0, 86.56, 0.234; GCA_003465495.1, s__Prevotella copri_G, 96.82, 86.55, 0.25; GCA_026015625.1, s__Prevotella copri_H, 96.58, 86.54, 0.27; GCF_026015465.1, s__Prevotella copri_E, 96.47, 86.45, 0.238; GCA_900555035.1, s__Prevotella sp900555035, 95.0, 86.33, 0.356; GCA_021622245.1, s__Prevotella copri_J, 95.95, 86.18, 0.245; GCA_026015295.1, s__Prevotella copri_C, 95.51, 86.14, 0.247; GCA_003463005.1, s__Prevotella copri_I, 96.73, 86.09, 0.242; GCF_025151535.1, s__Prevotella copri, 96.58, 86.04, 0.226; GCF_009494835.1, s__Prevotella copri_F, 96.41, 85.94, 0.237; GCA_015074785.1, s__Prevotella sp015074785, 96.07, 85.9, 0.261; GCA_949820605.1, s__Prevotella sp934191715, 96.61, 85.89, 0.247; GCF_009494135.1, s__Prevotella copri_B, 96.82, 85.89, 0.219; GCF_019127135.1, s__Prevotella sp900551275, 95.0, 85.5, 0.264; GCA_934196725.1, s__Prevotella sp934196725, 95.0, 85.21, 0.255; GCA_934202525.1, s__Prevotella sp934202525, 95.0, 85.19, 0.23   63.59   11      N/A     N/A
49647_1#70_clean_bin.6.strict   d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli       GCF_003697165.2     95.0    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli       96.85   0.843   N/A     N/AN/A      N/A     N/A     d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__       taxonomic classification defined by topology and ANI        topological placement and ANI have incongruent species assignments      GCF_003697165.2, s__Escherichia coli, 95.0, 96.85, 0.843; GCF_000194175.1, s__Escherichia coli_F, 95.0, 95.89, 0.887; GCF_002965065.1, s__Escherichia sp002965065, 95.0, 94.52, 0.708; GCF_004211955.1, s__Escherichia sp004211955, 95.0, 92.72, 0.777; GCF_005843885.1, s__Escherichia sp005843885, 95.0, 92.5, 0.796; GCF_029876145.1, s__Escherichia ruysiae, 95.0, 92.2, 0.767; GCF_011881725.1, s__Escherichia coli_E, 95.0, 92.19, 0.791; GCF_014836715.1, s__Escherichia whittamii, 95.0, 91.58, 0.779; GCF_000026225.1, s__Escherichia fergusonii, 95.0, 90.94, 0.564; GCF_002900365.1, s__Escherichia marmotae, 95.0, 90.83, 0.73; GCF_000759775.1, s__Escherichia albertii, 95.0, 90.03, 0.672        98.47   11      N/A     N/A
```
