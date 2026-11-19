# ATAC branch: align, filter, deduplicate, remove chrM and blacklist, call peaks
# (PRD section 7; ADR-001, ADR-009). ChIP shares these rules once its cards land.
BAM = RES + "/{study}/atac/bam/{sample}"
IDX = REF["bowtie2_index"]
_dna = config["dna_alignment"]


rule bowtie2_index:
    input:
        fasta=REF["fasta"],
    output:
        multiext(IDX, ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    params:
        prefix=IDX,
    threads: config["compute"]["threads_default"]
    log:
        RES + "/logs/bowtie2_index.log",
    conda:
        "../envs/align_dna.yaml"
    shell:
        "bowtie2-build --threads {threads} {input.fasta} {params.prefix} > {log} 2>&1"


rule bowtie2_align:
    input:
        reads=align_reads,
        index=multiext(IDX, ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    output:
        bam=temp(BAM + ".raw.bam"),
    params:
        prefix=IDX,
        preset=_dna["bowtie2"]["preset"],
        maxfrag=_dna["bowtie2"]["max_fragment_length"],
        extra=" ".join(_dna["bowtie2"]["extra_flags"]),
        reads=lambda wc, input: (
            f"-1 {input.reads[0]} -2 {input.reads[1]}" if wc.sample in PE_SAMPLES
            else f"-U {input.reads[0]}"
        ),
    threads: config["compute"]["threads_default"]
    log:
        bt2=RES + "/logs/{study}/{sample}.bowtie2.log",
    conda:
        "../envs/align_dna.yaml"
    shell:
        "bowtie2 -x {params.prefix} {params.reads} --{params.preset} -X {params.maxfrag} "
        "{params.extra} -p {threads} 2> {log.bt2} | samtools view -b -o {output.bam} -"


rule filter_pairs:
    input:
        BAM + ".raw.bam",
    output:
        temp(BAM + ".q30.bam"),
    params:
        mapq=_dna["filters"]["mapq_min"],
        flags=lambda wc: "-f 2 -F 1804" if wc.sample in PE_SAMPLES else "-F 1804",
    threads: config["compute"]["threads_default"]
    conda:
        "../envs/align_dna.yaml"
    shell:
        "samtools view -b -q {params.mapq} {params.flags} {input} "
        "| samtools sort -@ {threads} -o {output} -"


if _dna["duplicates"]["tool"] == "picard":

    rule dedup:
        input:
            BAM + ".q30.bam",
        output:
            bam=temp(BAM + ".dedup.bam"),
            metrics=RES + "/{study}/atac/bam/{sample}.dup_metrics.txt",
        log:
            RES + "/logs/{study}/{sample}.picard.log",
        conda:
            "../envs/align_dna.yaml"
        shell:
            "picard MarkDuplicates I={input} O={output.bam} M={output.metrics} "
            "REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT "
            "> {log} 2>&1"

else:

    rule dedup:
        input:
            BAM + ".q30.bam",
        output:
            bam=temp(BAM + ".dedup.bam"),
        threads: config["compute"]["threads_default"]
        conda:
            "../envs/align_dna.yaml"
        shell:
            "samtools sort -n -@ {threads} -o - {input} | samtools fixmate -m - - "
            "| samtools sort -@ {threads} -o - - | samtools markdup -r - {output.bam}"


rule index_dedup:
    input:
        BAM + ".dedup.bam",
    output:
        temp(BAM + ".dedup.bam.bai"),
    conda:
        "../envs/align_dna.yaml"
    shell:
        "samtools index {input}"


rule filter_regions:
    input:
        bam=BAM + ".dedup.bam",
        bai=BAM + ".dedup.bam.bai",
        blacklist=REF["blacklist"],
    output:
        bam=BAM + ".bam",
        bai=BAM + ".bam.bai",
    params:
        contigs=" ".join(config["reference"]["main_chromosomes"]),
    conda:
        "../envs/align_dna.yaml"
    shell:
        "samtools view -b {input.bam} {params.contigs} "
        "| bedtools intersect -v -abam stdin -b {input.blacklist} > {output.bam} "
        "&& samtools index {output.bam}"


rule macs3_callpeak:
    input:
        bam=BAM + ".bam",
        bai=BAM + ".bam.bai",
        chrom_sizes=REF["chrom_sizes"],
    output:
        narrowpeak=RES + "/{study}/atac/peaks/{sample}_peaks.narrowPeak",
        xls=RES + "/{study}/atac/peaks/{sample}_peaks.xls",
        summits=RES + "/{study}/atac/peaks/{sample}_summits.bed",
    params:
        flags=lambda wc: macs_flags(wc.sample),
        qvalue=config["peaks"]["qvalue"],
        gsize=lambda wc, input: genome_size(input.chrom_sizes),
        outdir=RES + "/{study}/atac/peaks",
    log:
        RES + "/logs/{study}/{sample}.macs3.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        "macs3 callpeak -t {input.bam} {params.flags} -g {params.gsize} -q {params.qvalue} "
        "--keep-dup all -n {wildcards.sample} --outdir {params.outdir} > {log} 2>&1"
