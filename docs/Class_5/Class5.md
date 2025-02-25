# Class 5: Ecological characterization of the microbiome

- - - -

## Taxon abundance and depth estimation (Bowtie2, InStrain)

[]

### Running bowtie2

[]

### Running instrain

[]

## Ecological indixes estimation and metadata analysis (R command line)

[]

### Coding the phyloseq and microviz R packages

[]

## Code templates

**`raw_reads_paths.txt` command:**

```bash
ls -1 -d /home/dorian.rojas/test/1-data/*_sub_?.fastq > raw_reads_paths.txt
```

**`raw_reads_paths.txt` file**

```vim
/home/dorian.rojas/test/1-data/SRR8555091_sub_1.fastq
/home/dorian.rojas/test/1-data/SRR8555091_sub_2.fastq
/home/dorian.rojas/test/1-data/SRR9988205_sub_1.fastq
/home/dorian.rojas/test/1-data/SRR9988205_sub_2.fastq
```

**`create-manifest.sh` file:**

```bash
#!/bin/bash
fastq_paths="$1"

echo "ID,R1,R2" > manifest.csv  # Headers

while read line; do
    if [[ $line == *"_sub_1.fastq" ]]; then  # Annotates R1
        read1="$line"
        # Annotates R2
        read2="${line%_sub_1.fastq}_sub_2.fastq"

        if [[ ! -f "$read2" ]]; then
            echo "Error: Cannot find corresponding R2 file for $read1" >&2
            continue  # next line in the input file
        fi

        lane_id=$(basename "${line%_sub_1.fastq}")
        echo "$lane_id,$read1,$read2" >> manifest.csv # Adds sample and raw data paths
    fi
done < "$fastq_paths"
```

**`manifest.csv` file:**

```vim
ID,R1,R2
SRR8555091,/home/dorian.rojas/test/1-data/SRR8555091_sub_1.fastq,/home/dorian.rojas/test/1-data/SRR8555091_sub_2.fastq
SRR9988205,/home/dorian.rojas/test/1-data/SRR9988205_sub_1.fastq,/home/dorian.rojas/test/1-data/SRR9988205_sub_2.fastq
```

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

```

**`instrain.slurm` file:**

```vim

```