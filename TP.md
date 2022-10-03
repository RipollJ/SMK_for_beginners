# TP snakemake

- Author: Julie Ripoll
- Contact: julie.ripoll87@gmail.com
- Date: 2022-09-14

Aim: First steps with snakemake

--------------------
## Installation

In this example we will use miniconda a package manager, to facilitate installation of programs and other tools.

Download miniconda at: https://docs.conda.io/en/latest/miniconda.html

(depend on your computer system)

Launch installation in a new terminal, here an example for Linux:
```bash
bash ./Miniconda3-py39_4.12.0-Linux-x86_64.sh
```

Update conda:
```bash
conda update conda
```

Install mamba (C++ implementation of conda, more efficient):
```bash
conda install -n base -c conda-forge mamba
```

Install snakemake (for Windows install snakemake-minimal instead):
```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Activate the snakemake environment:
```bash
conda activate snakemake
```

Try your installation:
```bash
snakemake --help
```

For more informations on conda command look at this [cheat sheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)


Create a directory for this exercise
```bash
mkdir exercise
# and move to this directory
cd exercise
```

## Download data

Download a sample with wget:
```bash
wget -P data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR026/SRR026762/SRR026762.fastq.gz
```
--------------
## Workflow

We will build a step-by-step snakemake workflow which takes two tools and home-made script.

The first tool is Cutadapt which cuts fastq files to get better quality reads.
Documentation: [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)

Create a snakefile and paste this rule into it:
```python
rule cutadapt:
    # Aim: removes adapter sequences from high-throughput sequencing reads
    # Use: cutadapt -a ADAPTER [options] [-o output.forward] [-p output.reverse]
     # <input.forward> <input.reverse>
    message:
        "cutadapt ---remove poly-A and adaptors--- on SRR026762"
    conda:
        "envs/quality.yml"
    input:
        "data/SRR026762.fastq.gz"
    output:
        "results/cutadapt/SRR026762.fastq.gz"
    params:
        quality = 35,
        length = 50,
        adapter = "'file:adapters.fa'"
    log:
        "results/cutadapt/SRR026762.log"
    shell:
        "cutadapt "
        "-a {params.adapter}' "
        "-q {params.quality} " # filter on quality thresold
        "-m {params.length} " # keep only read with minimal length defined and >
        "-o {output} "
        "{input} "
        "> {log}"
```
Note:
- this rule contains a **conda** flag with an environment called quality
- this rule uses one sample in input and return one sample in output
- here, the log file will contain the results print in the [terminal stdout](https://www.tutorialspoint.com/understanding-stdin-stderr-and-stdout-in-linux#)


-------------
We will build the conda file for the rule.

Create a conda environment configuration file "quality.yml" in a directory called envs.
```yaml
name: quality
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - multiqc=1.0
  - cutadapt=2.1
  - python=3.6
  - dataclasses
```
Note:
- the conda environment YAML file is composed of a name for the environment, channels for download from the online directory and the dependencies which contain all the required softwares.
- here, cutadapt is specified because as it will be used in the Snakefile rule

---------
How to execute a snakefile

as a first step, it is useful to test if the workflow is correctly defined and to estimate the amount of computation required.
The dry-run mode allows you to see the sequence of rules without running them.

Test running with the dry-run mode:
```bash
snakemake -s Snakefile --use-conda -j 1 -p -n
```

For direct running :
```bash
snakemake -s Snakefile --use-conda -j 1 -p
```

To remove outputs:
```bash
snakemake -s Snakefile -j 1 --delete-all-output
```
Notes:
- The -s tells Snakemake which snakefile to run.

- The --use-conda option tells Snakemake that a conda environment is required.

- The -j option tells Snakemake how many cores it can use (here only one job is executed), this a mandatory flag.

- The -p option tells Snakemake to also print the resulting shell command.

- The -n (or --dry-run) option, Snakemake will only show the running plan instead of actually running the steps.


------------
If we have more than one sample, snakemake can parallelize the rules so that they are executed per sample.
For this, download more samples in the data directory:
```bash
# install sra-tools
mamba install -c bioconda sra-tools
# download the NCBI archive
prefetch -O data/ SRR2931034
prefetch -O data/ SRR2931035
# extract fastq files
fastq-dump --split-3 -O SRR2931034/ --gzip data/SRR2931034/*.sra
fastq-dump --split-3 -O SRR2931035/ --gzip data/SRR2931035/*.sra

# reduce files for this TP
## just to gain time
gzip -dc data/SRR026762.fastq.gz | head -1000000 > data/SRR026762_1M.fastq
gzip -dc SRR2931034/SRR2931034.fastq.gz | head -1000000 > data/SRR2931034_1M.fastq
gzip -dc SRR2931035/SRR2931035.fastq.gz | head -1000000 > data/SRR2931035_1M.fastq
# zip them
gzip data/SRR026762_1M.fastq
gzip data/SRR2931034_1M.fastq
gzip data/SRR2931035_1M.fastq
```

-----------
Now, we can generalize the cutadapt rule to be used on multiple input files.
To do this, we will use wildcards.

Add a wildcard to the cutadapt rule:
```python
rule cutadapt:
    # Aim: removes adapter sequences from high-throughput sequencing reads
    # Use: cutadapt -a ADAPTER [options] [-o output.forward] [-p output.reverse]
     # <input.forward> <input.reverse>
    message:
        "cutadapt ---remove poly-A and adaptors--- on {wildcards.sample}"
    conda:
        "envs/quality.yml"
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/cutadapt/{sample}.fastq.gz"
    params:
        quality = 35,
        length = 50,
        adapter = "'file:adapters.fa'"
    log:
        "results/cutadapt/{sample}.log"
    shell:
        "cutadapt "
        "-a {params.adapter} "
        "-q {params.quality} " # filter on quality thresold
        "-m {params.length} " # keep only read with minimal length defined and >
        "-o {output} "
        "{input} "
        "> {log}"
```
Here, {sample} is the wildcard to be used.

This wildcard is a global wildcard and must be declared as such in the snakefile:
```python
SAMPLE, = glob_wildcards("data/{sample}.fastq.gz")
```

Now, Snakemake needs a rule called **all** or **defaults** to build the DAG of the jobs.

Add the rule all before the cutadapt rule:
```python
rule all:
    input:
        expand("results/cutadapt/{sample}.fastq.gz", sample = SAMPLE)

```
This rule takes as input the outputs of the cutadapt rule.
The expand function resolves the wildcards and allows snakemake to decompile the output tree.

---------------
Externalization of parameters is interesting for reusing the workflow without modifying the code.
For this, a configuration file can be used.

Create a configuration file "Config.yaml".
```yaml
## Params for CUTADAPT
cutadapt:
  quality: 35
  length: 50
  adapters: "'file:adapters.fa'"
```
Note: 
- a configuration file helps to reduce code errors by externalizing parameters that may differ between experiments.
- this step requires some modifications to the rule
- the configuration file can be declared in the Snakefile using the **configfile** tag: **"Config.yaml"** at the beginning of the file or at runtime in the command line.
```bash
snakemake -s Snakefile --use-conda --configfile Config.yaml -j 1 -p -n
```

Update the rule cutadapt:
```python
rule cutadapt:
    # Aim: removes adapter sequences from high-throughput sequencing reads
    # Use: cutadapt -a ADAPTER [options] [-o output.forward] [-p output.reverse]
     # <input.forward> <input.reverse>
    message:
        "cutadapt ---remove poly-A and adaptors--- on {wildcards.sample}"
    conda:
        "envs/quality.yml"
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/cutadapt/{sample}.fastq.gz"
    params:
        quality = config["cutadapt"]["quality"],
        length = config["cutadapt"]["length"],
        adapter = config["cutadapt"]["adapters"]
    log:
        "results/cutadapt/{sample}.log"
    shell:
        "cutadapt "
        "-a {params.adapter} "
        "-q {params.quality} " # filter on quality thresold
        "-m {params.length} " # keep only read with minimal length defined and >
        "-o {output} "
        "{input} "
        "> {log}"
```


Cutadapt can parallelize its process with an option for cores.
This parallelization can be added in the rule with the **threads** tag which takes the number of cores provided by the user.
You can declare the number of threads in the configuration file, to be externalized, which makes it easier to define a number of cores adapted to your computer.
It's also possible to use a function in the snakefile that automatically takes the number of available cores.

```python
rule cutadapt:
    # Aim: removes adapter sequences from high-throughput sequencing reads
    # Use: cutadapt -a ADAPTER [options] [-o output.forward] [-p output.reverse]
     # <input.forward> <input.reverse>
    message:
        "cutadapt ---remove poly-A and adaptors--- on {wildcards.sample}"
    conda:
        "envs/quality.yml"
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/cutadapt/{sample}.fastq.gz"
    params:
        quality = config["cutadapt"]["quality"],
        length = config["cutadapt"]["length"],
        adapter = config["cutadapt"]["adapters"]
    threads:
        3
    log:
        "results/cutadapt/{sample}.log"
    shell:
        "cutadapt "
        "-a {params.adapter} "
        "-q {params.quality} " # filter on quality thresold
        "-m {params.length} " # keep only read with minimal length defined and >
        "-o {output} "
        "-j {threads} "
        "{input} "
        "> {log}"
```

-----------
Now, we will create a script to test the aggregation rules.
In a terminal, create a new bash script and open it:

```bash
touch aggregation_test.py
# open with viw or other code editor
vim aggregation_test.py
```

Our goal is to see how the aggregation works in snakemake, using the **expand()** function or with a for loop.

Paste this code inside:
```python
#!/bin/python

from contextlib import redirect_stdout

with open(snakemake.output[0], 'w') as f:
    with redirect_stdout(f):
        print(snakemake.input)
        print(type(snakemake.input))
```

And a new rule in your snakefile
```python
rule test_aggregation_expand:
    # Aim: How works the aggregation
    input:
        expand("results/cutadapt/{sample}.fastq.gz", sample = SAMPLE)
    output:
        "results/aggregation_test_expand.txt"
    script:
        'aggregation_test.py'
```

Add the second test with the for loop:
```python
rule test_aggregation_loop:
    # Aim: How works the aggregation
    input:
        ["results/cutadapt/{sample}.fastq.gz".format(sample=sample) for sample in SAMPLE]
    output:
        "results/aggregation_test_loop.txt"
    script:
        'aggregation_test.py'
```

Now, you need to change the outputs in the rule **all**.
The outputs of cutadapt are not needed because the wildcard is resolved in the aggregation rules.

```python
rule all:
    input:
        test1 = "results/aggregation_test_expand.txt",
        test2 = "results/aggregation_test_loop.txt"

```
Note: this rule takes as input the outputs of the two independant aggregation rules.


----------
Create the DAG:

Snakemake specifies the DAG in "dot" language, using dot from Graphviz.

Install Graphviz in your snakemake environment:
```bash
conda install -c anaconda graphviz
```

Export the DAG of your workflow:
```bash
snakemake -s Snakefile --configfile Config.yaml --dag | dot -Tsvg > dag.svg
```

Execute the snakefile
```bash
snakemake -s Snakefile --use-conda --configfile Config.yaml -j 2 -p
```

-----------
## Extraction of reports

You can create automatic report from your snakemake directory, which is usefull for benchmarking, with this command:

```bash
snakemake -s Snakefile --configfile Config.yaml --report report.html
```

All information contained in the report (e.g. runtime statistics) is automatically collected after the snakefile is executed.

You can also define a specific output to add to the report file by declaring **report()** in the output.

This is an example:
```python
rule plot_commits_and_releases:
    input:
        "tables/git-log.csv"
    output:
        report("plots/commits+releases.svg", category="Plots", caption="report/commits+releases.rst")
    conda:
        "envs/stats.yaml"
    notebook:
        "notebooks/commits+releases.py.ipynb"
```
We do not apply it in this TP but I invite you to test it in one of your pipeline.

---------------

Some tools are useful to exploit the log files or outputs containing statistics. MultiQC summarizes the statistics of the steps performed by a tool.
Documentation: [MultiQC](https://multiqc.info/)

Add another rule to the snakefile:

```python
rule MultiQC:
    # Aim: Check quality of the results from another tool.
    message:
        """--- MultiQC reports on cutadapt outputs ---"""
    conda:
        "envs/quality.yml"
    input:
        expand("results/cutadapt/{sample}.fastq.gz", sample = SAMPLE)
    params:
        filename = "report_cutadapt.html", # report name
        module = 'cutadapt', # tool previously used
        outfolder = 'results/cutadapt/*.log' # log file
    output:
        directory('results/Quality/MQC_cutadapt/')
    log:
        'results/Quality/MQC_cutadapt.done'
    shell:
        "multiqc "
        "-n {params.filename} "
        "-m {params.module} "
        "-d {params.outfolder} "
        "--outdir {output} "
        "2>{log}"
```
Note:
- This new rule is dependant from the cutadapt rule because it uses its outputs files.
- This rule takes an argument in the output: **directory()**, this tag defines the type of output. Other tags can be used like **temp()** for temporary file that are deleted after the rule is executed.
- There is no wildcard in the output because of the use of **expand()** in the input. Expand() gives all the samples as input to the rule and not one by one. It's an aggregation.


Remember to add the outputs to the rule **all**.

-----------
## Checkpoints

One might have to filter out samples that do not pass the quality check (QC).

Because QC is an intermediate result of the same data analysis, it may be necessary to determine the part of the DAG that is downstream of QC only after QC has been finalized.

One option is to separate the QC and the actual analysis into two workflows, or to define a separate target rule for the QC, so that it can be performed manually upstream, before the actual analysis is started.

Alternatively, if QC is to occur automatically as part of the overall workflow, one can use Snakemake's conditional execution capabilities, using the checkpoint.

It is also possible to use checkpoints for cases where the output files are unknown before execution.

Here, we will check the outputs of an assembler. An assemblera allows to create a new reference genome or transcriptome based on your samples. They globally used iterative cycle of DeBruyn graph with different k-mer lengths.

Create a new environment called assembly.yml in the envs directory:
```yaml
name: assembly
channels:
  - bioconda
  - defaults
dependencies:
  - megahit
```

Now, add these rules to your snakefile and add the output to the rule all.
```python
rule commalist:
    # Aim: list input files in a comma separated list
    message:
        """--- Convert dict input files in a comma separated list ---"""
    input:
        expand("results/cutadapt/{sample}.fastq.gz", sample = SAMPLE)
    output:
        "results/cutadapt/list_samples.txt"
    shell:
        "echo {input} | sed ':a;N;$!ba;s/\\n/,/g' > {output}"

checkpoint assembly:
    # Aim: assembly of reads into contigs.
    message:
        """--- Megahit assembly in progress ---"""
    conda:
        "envs/assembly.yml"
    input:
        "results/cutadapt/list_samples.txt"
    threads:
        3
    params:
        klist = "29,39,59,79,99", # tool range 15-255
        kmin = "21",
        kmax = "141"
    output:
        directory("megahit_out")
    shell:
        "megahit "
        "-r $(paste {input}) "
        "-o {output} "
        "-t {threads} "
        "--k-list {params.klist} "
        "--k-min {params.kmin} "
        "--k-max {params.kmax} "

def get_file_names(wildcards):
    ck_output = checkpoints.assembly.get(**wildcards).output[0]
    global SMP
    SMP, = glob_wildcards(os.path.join(ck_output, "intermediate_contigs/{name}.contigs.fa"))
    return expand(os.path.join(ck_output, "intermediate_contigs/{name}.contigs.fa"), name=SMP)

rule check_contigs:
    input:
        get_file_names
    output:
        "results/contigs/summary.txt"
    shell:
        "wc -l {input} > {output}"

```
Note: 
- the function creates the new wildcards from the output in the directory of the assembly rule.
- you must specify the name of the rule in the function from which the new wildcards will be created.
- here, the next rule entry is a call to this control function. But it can also be in the rule all, if the wildcards are used in other rules.


-----------
Here is another example (from the snakemake documentation) combined with pipe() tag:
```python
rule all:
    input:
        "c.txt",

checkpoint a:
    output:
        "a.txt"
    shell:
        "touch {output}"

rule b1:
    output:
        pipe("b.pipe")
    shell:
        ""

rule b2:
    input:
        pipe("b.pipe")
    output:
        "b.txt"
    shell:
        "touch {output}"

rule c:
    input:
        "a.txt",
        "b.txt",
    output:
        "c.txt"
    shell:
        "touch {output}"
```

The **pipe()** tag allows you to pass the output from one rule to another in the stream without writing it to disk. This is different from **temp()** which writes the output and then deletes it.