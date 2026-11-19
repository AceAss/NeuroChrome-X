# QC statistics, the QC table and MultiQC (PRD section 7; ADR-009).
_qc = config["qc"]["flags"]


rule stage_stats:
    input:
        raw=BAM + ".raw.bam",
        q30=BAM + ".q30.bam",
        dedup=BAM + ".dedup.bam",
        dedup_bai=BAM + ".dedup.bam.bai",
        final=BAM + ".bam",
        final_bai=BAM + ".bam.bai",
        peaks=RES + "/{study}/atac/peaks/{sample}_peaks.narrowPeak",
        gtf=REF["gtf"],
    output:
        RES + "/{study}/atac/stats/{sample}.tsv",
    params:
        layout=lambda wc: "PE" if wc.sample in PE_SAMPLES else "SE",
        contigs=" ".join(config["reference"]["main_chromosomes"]),
        mito=config["reference"]["mito_contig"],
        bt2_log=RES + "/logs/{study}/{sample}.bowtie2.log",
    conda:
        "../envs/qc.yaml"
    shell:
        "python {SCRIPTS}/py/make_qc_table.py sample --sample {wildcards.sample} "
        "--layout {params.layout} --raw-bam {input.raw} --q30-bam {input.q30} "
        "--dedup-bam {input.dedup} --final-bam {input.final} --bowtie2-log {params.bt2_log} "
        "--peaks {input.peaks} --gtf {input.gtf} --keep-contigs {params.contigs} "
        "--mito {params.mito} --out {output}"


rule qc_table:
    input:
        lambda wc: expand(
            RES + "/{study}/atac/stats/{sample}.tsv", study=wc.study, sample=ATAC_BY_STUDY[wc.study]
        ),
    output:
        RES + "/{study}/atac/qc/qc_table.tsv",
    params:
        frip=_qc["frip_min"],
        tss=_qc["tss_enrichment_min"],
        nrf=_qc["nrf_min"],
        pbc1=_qc["pbc1_min"],
    conda:
        "../envs/qc.yaml"
    shell:
        "python {SCRIPTS}/py/make_qc_table.py table --stats {input} --frip-min {params.frip} "
        "--tss-min {params.tss} --nrf-min {params.nrf} --pbc1-min {params.pbc1} --out {output}"


if config["qc"].get("multiqc", False):

    rule multiqc:
        input:
            fastp=lambda wc: expand(
                RES + "/{study}/atac/trim/{sample}.fastp.json",
                study=wc.study,
                sample=ATAC_BY_STUDY[wc.study],
            ),
            bt2=lambda wc: expand(
                RES + "/logs/{study}/{sample}.bowtie2.log",
                study=wc.study,
                sample=ATAC_BY_STUDY[wc.study],
            ),
        output:
            RES + "/{study}/atac/qc/multiqc_report.html",
        params:
            outdir=RES + "/{study}/atac/qc",
        log:
            RES + "/logs/{study}/multiqc.log",
        conda:
            "../envs/qc.yaml"
        shell:
            "multiqc -f -o {params.outdir} -n multiqc_report.html {input.fastp} {input.bt2} "
            "> {log} 2>&1"
