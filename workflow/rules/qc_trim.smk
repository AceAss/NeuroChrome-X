# Read QC and adapter trimming (PRD section 7: adapter trimming only).
if trimming_enabled():
    _trim_flags = "--detect_adapter_for_pe" if config["trimming"]["detect_adapters"] else ""
    _quality = "" if config["trimming"]["quality_trimming"] else "--disable_quality_filtering"

    rule trim_fastp_pe:
        input:
            r1=lambda wc: fastq_paths(wc.sample)[0],
            r2=lambda wc: fastq_paths(wc.sample)[1],
        output:
            r1=temp(RES + "/{study}/atac/trim/{sample}_R1.fastq.gz"),
            r2=temp(RES + "/{study}/atac/trim/{sample}_R2.fastq.gz"),
            json=RES + "/{study}/atac/trim/{sample}.fastp.json",
            html=RES + "/{study}/atac/trim/{sample}.fastp.html",
        wildcard_constraints:
            sample=wildcard_alternatives(PE_SAMPLES),
        params:
            adapter=_trim_flags,
            quality=_quality,
            min_len=config["trimming"]["min_length"],
        threads: config["compute"]["threads_default"]
        log:
            RES + "/logs/{study}/{sample}.fastp.log",
        conda:
            "../envs/qc.yaml"
        shell:
            "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
            "{params.adapter} {params.quality} --length_required {params.min_len} "
            "--thread {threads} -j {output.json} -h {output.html} > {log} 2>&1"

    rule trim_fastp_se:
        input:
            r1=lambda wc: fastq_paths(wc.sample)[0],
        output:
            r1=temp(RES + "/{study}/atac/trim/{sample}_R1.fastq.gz"),
            json=RES + "/{study}/atac/trim/{sample}.fastp.json",
            html=RES + "/{study}/atac/trim/{sample}.fastp.html",
        wildcard_constraints:
            sample=wildcard_alternatives(SE_SAMPLES),
        params:
            quality=_quality,
            min_len=config["trimming"]["min_length"],
        threads: config["compute"]["threads_default"]
        log:
            RES + "/logs/{study}/{sample}.fastp.log",
        conda:
            "../envs/qc.yaml"
        shell:
            "fastp -i {input.r1} -o {output.r1} {params.quality} "
            "--length_required {params.min_len} --thread {threads} "
            "-j {output.json} -h {output.html} > {log} 2>&1"


if config["qc"].get("fastqc", False):

    rule fastqc:
        input:
            lambda wc: fastq_paths(wc.sample),
        output:
            touch(RES + "/{study}/atac/fastqc/{sample}.done"),
        params:
            outdir=RES + "/{study}/atac/fastqc",
        threads: config["compute"]["threads_default"]
        log:
            RES + "/logs/{study}/{sample}.fastqc.log",
        conda:
            "../envs/qc.yaml"
        shell:
            "fastqc -t {threads} -o {params.outdir} {input} > {log} 2>&1"
