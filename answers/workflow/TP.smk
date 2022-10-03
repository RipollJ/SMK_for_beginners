
SAMPLE, = glob_wildcards("data/{sample}.fastq.gz")


rule all:
    input:
        test1 = "results/aggregation_test_expand.txt",
        test2 = "results/aggregation_test_loop.txt",
        test3 = 'results/Quality/MQC_cutadapt.done',
        test4 = "results/contigs/summary.txt"


rule cutadapt:
    # Aim: removes adapter sequences from high-throughput sequencing reads
    # Use: cutadapt -a ADAPTER [options] [-o output.forward] [-p output.reverse]
     # <input.forward> <input.reverse>
    message:
        "cutadapt ---remove poly-A and adaptors--- on {wildcards.sample}"
    conda:
        "../envs/quality.yml"
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/cutadapt/{sample}.fastq.gz"
    params:
        quality = config["cutadapt"]["quality"],
        length = config["cutadapt"]["length"],
        adapter = config["cutadapt"]["adapters"]
    threads:
        8
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

rule test_aggregation_expand:
    # Aim: How works the aggregation
    input:
        expand("results/cutadapt/{sample}.fastq.gz", sample = SAMPLE)
    output:
        "results/aggregation_test_expand.txt"
    script:
        '../scripts/aggregation_test.py'


rule test_aggregation_loop:
    # Aim: How works the aggregation
    input:
        ["results/cutadapt/{sample}.fastq.gz".format(sample=sample) for sample in SAMPLE]
    output:
        "results/aggregation_test_loop.txt"
    script:
        '../scripts/aggregation_test.py'


rule MultiQC:
    # Aim: Check quality of the results from another tool.
    message:
        """--- MultiQC reports on cutadapt outputs ---"""
    conda:
        "../envs/quality.yml"
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

### second part

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
        "../envs/assembly.yml"
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
    SMP, = glob_wildcards(os.path.join(ck_output, "intermediate_contigs/k{name}.contigs.fa"))
    return expand(os.path.join(ck_output, "intermediate_contigs/k{name}.contigs.fa"), name=SMP)

rule check_contigs:
    input:
        get_file_names
    output:
        "results/contigs/summary.txt"
    shell:
        "wc -l {input} > {output}"
