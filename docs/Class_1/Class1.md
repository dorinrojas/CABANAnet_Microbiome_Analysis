# Class 1: Basic Unix shell Bash commands and workshop methodology

- - - -

## Using Unix through the terminal

Unix shell is a command-line interface (CLI) platform and scripting language initially developed to construct other softwares. For example, the iOS software is written in a derivated of Unix language. This is an efficient universal operating system used in many computers as it allows to run several softwares and small programs.

There are different types of Unix shells. However, the most popular among computers is Bash (Bourne Again SHell). Computational clusters are based on Bash. Therefore, it is important to understand the main commands of the CLI.

We are used to interact with computers in our daily life through a graphical user interface (GUI). Here, we provide instructions by clicking in folders and other directions that are learned intuitively. Contrary to this, the shell lacks the a visual presentation and works as a command line, which might come as a shock for new users.

### Accessing the Unix shell

The shell is a more efficient way to communicate and work with a computer. It can be accessed through the 'terminal' (also called 'Command prompt'). A terminal initiates with the `$` character. This is called a prompt and is an indicator of the shell waiting for an input (the same reason it is also called 'Command prompt'). Prior to the dollar symbol, it is common to find the username of your computer.

```bash
[dorian.rojas@accessnode test]$
```

Contratry to Apple computers, Windows is built under a different system; thus, the commands we are going to explore do not work in the windowns shell. However, the start of a Windows shell looks similar to this:

```ps
C:\Users\rojas>
```

In order to be working on the same operating system, we are going to enter the computational cluster through a SSH client. This is a sofware that enables the connection to remote computers (in this case, our cluster) by applying a secure shell protocol. Most computers have the SSH client installed within the terminal. However, in some cases, Windows users have to use and external client.

> If the SSH client is not installed, we recommend installing the MobaXterm server. Please follow the instructions for downloading the software [here](https://mobaxterm.mobatek.net/download.html).

To connect to the cluster, open a terminal and input the following command:

```ps
ssh [user]@172.16.24.2
```

`[user]` must be changed for the username you were provided with. For example:

```ps
ssh dorian.rojas@172.16.24.2
```

In order to access the remote computer, your confidential password is required. As a precaution, the terminal does not show the password input. Once the password is correctly typed, the welcome message of the HPC cluster of the University of Costa Rica appears in screen.

```bash
Register this system with Red Hat Insights: insights-client --register
Create an account or view all your systems at https://red.ht/insights-dashboard
Last login: __
[dorian.rojas@accessnode ~]$
```

You can open several terminal at the same time. Most people are used to work with two or three simultaneously.

#### Bash command composition

A shell input is composed by the prompt (`$`), the command (e.g. `ls`), flags (or options, e.g.`-h`), and arguments (indications, e.g. a directory, file). Moreover, flags do not necessarily require arguments. These can be placed together in a single `-` (dash). For instance, `ls -lh` holds two flags, the `-l` for listing `-h` for transforming data into human-readble. Here, `ls -l -h` is a synonym of `ls -lh`.

### Basic Bash commands

Here, some main Bash commands and exercises to get comfortable with the 'black screen' are presented. Remember, the Unix command line is case sensitive.

Command|Function|Command|Function
-------|--------|-------|--------
`pwd`|Returns path|`ls`|Lists directories and files
`cd <path>`|Changes dir|`ls -lh`|Lists with human-readable data
`cd ..`|Goes to parent dir|`cp <file> <file_copy>`|Copies a file
`rm <file>`|Removes a file|`cp -r <dir> <dir_copy>`|Copies a dir
`rm -r <dir>`|Removes a dir|`mv <file> <path>`|Moves file to path
`mkdir`|Creates a dir|`mv <file/dir> <new_name>`|Changes name
`nano <file>`|Creates an empty file|`sbatch <file.slurm>`|Runs .slurm file
`squeue`|Shows running jobs|`scancel <jobID>`|Cancels job

Try some commands, input these in your terminal:

```bash
[dorian.rojas@accessnode curso]$ ls
[dorian.rojas@accessnode curso]$ mkdir test
[dorian.rojas@accessnode curso]$ cd test/
[dorian.rojas@accessnode test]$
```

The `ls` command does not output anything. This is because the account or folder is empty. The following command creates a new directory named 'test' (`code`) and, subsequently, it moves to that directory (`cd test`). Note the change in the name prior to the $, it indicates the current directory. To obtain the full path, use `pwd`.

```bash
[dorian.rojas@accessnode test]$ pwd
/home/dorian.rojas/curso/test
```

> As mentioned above, the shell terminal is case sensitive and takes in consideration spaces as well. In a GUI, we name files and folders using spaces, periods, among other special characters. It is recommended that when working in Bash, names:
    1. Don't have spaces. This will complete change the meaning of the command.
    2. Don't start with - (dash).
    3. Stick to letters and numbers.
    4. If required, - (dash) or _ (underscore) can be used instead of spaces (e.g. `dorian_rojas` rather than `dorian rojas`).

Moreover, most commands have a help menu for flags that can be used (e.g. `ls -h`). In order to see all flags, the option `--help` must be indicated after the command.

```bash
[dorian.rojas@accessnode test]$ ls --help
Usage: ls [OPTION]... [FILE]...
List information about the FILEs (the current directory by default).
Sort entries alphabetically if none of -cftuvSUX nor --sort is specified.

Mandatory arguments to long options are mandatory for short options too.
  -a, --all                  do not ignore entries starting with .
  -A, --almost-all           do not list implied . and ..
      --author               with -l, print the author of each file
  -b, --escape               print C-style escapes for nongraphic characters
      --block-size=SIZE      with -l, scale sizes by SIZE when printing them;
                               e.g., '--block-size=M'; see SIZE format below
  -B, --ignore-backups       do not list implied entries ending with ~
  -c                         with -lt: sort by, and show, ctime (time of last
                               modification of file status information);
                               with -l: show ctime and sort by name;
                               otherwise: sort by ctime, newest first
[...]
```

When working with new tools this command becomes relevant as it allows to have an overview of the software's usage and requirements. In bioinformatics, displaying the 'how to run' menu is essential.

#### Removing command: careful while deleting files

Another relevant command is the remove command (`rm`). In Unix shell, deletion of files is permanent and there is no way to recover them after removal. An interesting approach to have a double confirmation when deleting is running `rm` with the flag `-i`.

```bash
[dorian.rojas@accessnode curso]$ rm -i note.txt
rm: remove regular file 'note.txt'? n
[dorian.rojas@accessnode curso]$ rm -ir test/
rm: remove directory 'test/'? n
[dorian.rojas@accessnode curso]$ ls -F
note.txt  test/
[dorian.rojas@accessnode curso]$
```

The `-i` flags prompts a removal confirmation that has to be answer using `Y` or `N` to delete or cancel the command, respectively. This gives the user a chance to check the files and directories prior to deletion.

#### Basic Bash commands tasks

Explore other basic commands individually by completing the following tasks. The answers or code samples to solve each task are at the final section of current page.

1. Go back to the parent directory
2. Print the complete path where you are located
3. Change the name of the directory 'test' to 'metagenomics'
4. Delete the 'metagenomics' directory
5. List the current running jobs

#### Text editing commands

Editing text files is another relevant function to understand in Bash as you have to repeatibly type and edit code. There are a lot of text editors (e.g. vi, vim, nano, sublime). We'll use 'nano', which is newer, simpler, and easier than other editor like vim. The basic command is `nano <file_name.txt>`.

```bash
nano note.txt
```

This command opens a screen where you can type. The lower section of the screen shows the menu options of the nano editor. It is key to learn how to save and close the document. For this, you have to press `ctrl + X`. If the document is empty, the screen will close and nothing will be created in the current directory (check with `ls`).

When something is typed in the screen, the `ctrl + X` command will prompt a question to 'Save modified buffer'. To save the changes make to the document, the user must type  `Y`. This will pop-up another inquiry to modify the name of the file. However, the basic commands contains the argument for file name and this can be skipped by typing `enter`.

In summary, to save a nano document with modification the command `ctrl + X, Y, enter` is required. Check the saved document in your current directory with `ls`.

> Files open in nano can be saved while typing using `ctrl + O` and then exiting.

```bash
[dorian.rojas@accessnode curso]$ nano note.txt
[dorian.rojas@accessnode curso]$ ls
note.txt
[dorian.rojas@accessnode curso]$
```

There are other text file commands.

Command|Function|Command|Function
-------|--------|-------|--------
`cat <file>`|Prints complete file|`head <file>`|Prints the first lines of file
`tail <file>`|Prints last lines of file|`more <file>`|Prints little by little, exiting with `q`
`wc <file>`|Counts words in file|`wc -l`|Counts lines in file
`cmd > <file>`|Saves STDOUT in file|`cat <file1> <file2> > <file3>`|Concatenates files into file 3
`grep "pattern" <file>`|Prints lines matching pattern|`awk`|Very diverse command for text manipulation

#### Text commands tasks

Explore the text commands through the following tasks.

1. Reopen the file and write 10 separate lines with different cities' name
2. Print the complete document
3. Print the file's last lines
4. Print the file's first lines
5. Try printing the lines little by little
6. Create a new document with 10 different people names in separate lines
7. Concatenate the cities and names documents
8. Print only one line corresponding to one of the names in the concatenated document
9. Count the number of lines in the concatenated document

#### Other useful Bash commands

Bash has other several commands.

Command|Function|Command|Function
-------|--------|-------|--------
`gzip`|Compresses into .gzip|`zip`|Compresses into .zip
`gunzip`|Decompresses a .gzip|`unzip`|Decompresses a .zip
`tar -cf <dir>`|Compresses dir into .tar|`tar -x <.tar>`|Decompresses a .tar
`wget <URL>`|Downloads file|`top`|Shows activity
`scp <server>:<path> <path>`|Copies files between computers|`exit`|Exits remote connection

Feel free to try and run these commands independently. Be careful and ask questions to your instructor if required.

#### Wildcards

Bash allows the use of wildcards to select multiple files simultaneously. The two most important wildcards are `*` and `?`. The asterik `*` represent zero or more characters, whereas `?` represents exactly one character. For a better understading, set a directory with the files `methane.txt` and `ethane.txt`. For instance, `*ethane.txt` would indicate both files, while `?ethane.txt` only `methane.txt`.

```bash
[dorian.rojas@accessnode curso]$ ls *ethane.txt
ethane.txt  methane.txt
[dorian.rojas@accessnode curso]$ ls ?ethane.txt
methane.txt
```

Wildcards are relevant in bioinformatics as they help to manage large datasets. For example, if you own a directory full of .fastq files and want to concatenate them into a single document, the command `cat *.fastq` would be efficient. If not using wildcards, each files name would have to be typed individually, which might become exhaustive if there are more than 5 files with large and complex names.

## Workshop methodology and slurm files

### Slurm files configuration and template

During this workshop, we will be working with slurm files. Slurm is a workload manager that enables the correct organization of jobs running in the computational cluster. Slurm files are defined with `.slurm` extentions. Briefly, they are normal text files with a characteristic header.

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="fastq"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL
```

The first line `#!/bin/bash` is the command that allows the computer to recognize the document as a .slurm file. The following lines `#SBATCH --partition=parallel` and `#SBATCH --account=parallel-24h` represent cluster/specific configurations. It is important to remark that all computational cluster are organized differently. Therefore, the introduction to the cluster is essential for its correct usage.

> The `#!/bin/bash` header is not a slurm-specific feature. Other commands to run code (e.g. `sh`) also require this header in the code file.

The next lines of the file represent the job's maximum running time, number of computation nodes to use, and the number of threads.

```vim
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
```

Subsequently, the job name and output files (commonly known as log files) are specified. Due to personal preferences, logs files are divided in error file (`.e`) and output file (`.o`). This allows an easy access to errors in case of failure.

```vim
#SBATCH --job-name="fastq"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
```

These last configs allow the slurm software to send emails indicating if the job has began, ended, or failed. This param is set to only `END` and `FAIL`. However, to receive a "began" email, the option `BEGIN` can be added. Also, you can simplify this by changing the options to `ALL`, which specifies to all of the above email types.

```vim
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL
```

> Note that most of these are personal configurations. The slurm software has a lot of `#SBATCH` flags that can be added and modified as preferred. If you have already worked with slurm and are incline to other configurations, feel free to modify the header accordingly.

Next, the .slurm file required to set the working directory using a `cd` command. Additionally, for most tools used in this workshop we are require to set the path to the singularity containers (`CTN_PATH`). Although this later is optional, it will simplify the code by avoiding the complete path.

```vim
cd /home/dorian.rojas/
CTN_PATH=/opt/ohpc/pub/containers/BIO/
```

This completes the main configuration to the template of the `.slurm` file. Following lines are for type the code of the specific tools running in that job. For instance, this is the complete `.slurm` file to the tool `fastq-dump` tool. This software allows to download data directly from the NCBI using SRA accessions.

```vim
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --account=parallel-24h
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="fastq"
#SBATCH -o zz-%x-%j.o
#SBATCH -e zz-%x-%j.e
#SBATCH --mail-user=dorian.rojas@ucr.ac.cr
#SBATCH --mail-type=END,FAIL

cd ./1-data/

CTN_PATH=/opt/ohpc/pub/containers/BIO/

$CTN_PATH/sra-tools-3.1.1.sif fastq-dump --version

for id in $@; do

echo "Downloading " $id
$CTN_PATH/sra-tools-3.1.1.sif fastq-dump --gzip \
        --skip-technical --readids --dumpbase --split-3 \
        --clip -O . $id

echo $id " downloaded"

done

date
time
```

### Workship methodology

For the workshop, two sample will be used to ensure the optimization of computational resources and avoid saturating the cluster. In order to have a better automatization of the pipeline, the slurm files are written under a simple interconection of different files. Automatization of workflows can be performed using other more efficient softwares (e.g. Nextflow, Snakemake). However, these require prior knowlegde on how to code the automatized pipeline. Therefore, we will use simple Bash commands.

Three files are required: a `accessions.txt`, a `batch.sh`, and a `.slurm`. The first document contains the samples' accessions in separate lines. The second file extract the accession from the first file using a `sed` command and prompts the `.slurm` to run using that accession as the input.

```bash
[dorian.rojas@accessnode test]$ cat accessions.txt
SRR9988205
SRR8555091

[dorian.rojas@accessnode test]$ cat batch.sh
#!/bin/bash

sample1="$(sed -n '1p' accessions.txt)"
sample2="$(sed -n '2p' accessions.txt)"

sbatch read_qc.slurm $sample1
sbatch read_qc.slurm $sample2
```

The slurm holds the job-specific code. In order to use the accession sent from the `batch.sh` command, the `.slurm` file is coded on a _for loop_. Loops are programming structures that allow the iteration of commands a specific number of times. In our case, the loop runs for the number of accession indicated in `batch.sh`.

Note `sed -n '1p' accessions.txt` only selects the first line of the `.txt`, meaning one single accession. However, this could be modified for any number of accessions. For instance, `sed -n '1,8p' accessions.txt` would run eight samples. The `batch.sh` file is coded for only one accession due to the number of sample we are analyzing.

With in `batch.sh`, it is indicated that the `.slurm` must run the _for loop_ using the sample accessions as sample name (`sbatch read_qc.slurm $sample1`). Therefore, in the `.slurm` file, the command `for id in $@; do` specifies to run the job for the number of ids (`$id`) present in `$sample1` (`$@`).

> In Bash, to call an object, the prompt must be typed before the name (e.g. `$id`, `$CTN_PATH`).

```vim
for id in $@; do

echo "Downloading " $id
$CTN_PATH/sra-tools-3.1.1.sif fastq-dump --gzip \
        --skip-technical --readids --dumpbase --split-3 \
        --clip -O . $id

echo $id " downloaded"

done
```

Enclosing the object's name in braces adds the following strings to the name in the code lecture. For example:

```bash
id = sample1
${id}.fastq = sample1.fastq
```

This comes handy in following codes when extensions and output names are required. For instance, during the quality control, input arguments requires the extension of each paired reads file (`${sample}_1.fastq`,`${sample}_2.fastq`).

```vim
for id in $@; do

echo "Working on " $id

mkdir 2-read_qc/$sample

metawrap read_qc -1 1-data/${sample}_1.fastq \
        -2 1-data/${sample}_2.fastq \
        -t 64 -o 2-read_qc/$sample

echo $id " done"

done
```

Finally, the sample selection based on `accession.txt` and `batch.sh` files allows the `.slurm` command to run throughout multiple samples, regardless of the name. To run the same job with completely differente sample, user have to create the `accession.txt` and `batch.sh` files with their respective new information. This permits to analyze samples more efficiently as the main `.slurm` command will only have to be created from scratch once. Some small modifications might still be needed (e.g. changing the working directory path, files paths).

Feel free to ask any questions regarding this scheme. It is essential to understand how this interplay works to complete the workshop. However, the practice will get you used to the commands.

## Task solutions

**Basic commands tasks:**

1. Go back to the parent directory

    ```bash
    cd ..
    ```

2. Print the complete path where you are located

    ```bash
    pwd
    ```

3. Change the name of the directory 'test' to 'metagenomics'

    ```bash
    mv test metagenomics
    ```

4. Delete the 'metagenomics' directory

    ```bash
    rm -r metagenomics
    ```

5. List the current running jobs

    ```bash
    squeue
    ```

**Text file commands task:**

1. Reopen the file and write 10 separate lines with different cities' name

    ```bash
    nano note.txt
    ```

2. Print the complete document

    ```bash
    cat note.txt
    ```

3. Print the file's last lines

    ```bash
    tail note.txt
    ```

4. Print the file's first lines

    ```bash
    head note.txt
    ```

5. Try printing the lines little by little

    ```bash
    more note.txt
    ```

6. Create a new document with 10 different people names in separate lines

    ```bash
    nano names.txt
    ```

7. Concatenate the cities and names documents

    ```bash
    cat note.txt names.txt > concat.txt
    ```

8. Print only one line corresponding to one of the names in the concatenated document

    ```bash
    grep "dorian" concat.txt
    ```

9. Count the number of lines in the concatenated document

    ```bash
    wc -l concat.txt
    ```
