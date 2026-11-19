# Consensus peaks, counts and DESeq2 contrasts (PRD sections 7 and 8; ADR-008).
rule group_sheet:
    output:
        RES + "/{study}/{assay}/group_sheet.tsv",
    params:
        samples=config["paths"]["samples"],
    conda:
        "../envs/py_core.yaml"
    shell:
        "python {SCRIPTS}/py/write_group_sheet.py --samples {params.samples} "
        "--study {wildcards.study} --assay {wildcards.assay} --out {output}"


rule consensus_peaks:
    input:
        sheet=RES + "/{study}/atac/group_sheet.tsv",
        peaks=lambda wc: expand(
            RES + "/{study}/atac/peaks/{sample}_peaks.narrowPeak",
            study=wc.study,
            sample=ATAC_BY_STUDY[wc.study],
        ),
    output:
        RES + "/{study}/atac/consensus/consensus.bed",
    params:
        pattern=lambda wc: f"{RES}/{wc.study}/atac/peaks/{{sample}}_peaks.narrowPeak",
        build=config["reference"]["genome_build"],
        min_rep=config["peaks"]["consensus"]["min_replicates"],
    conda:
        "../envs/peaks.yaml"
    shell:
        "python {SCRIPTS}/py/consensus_peaks.py --samples {input.sheet} "
        "--peaks-pattern '{params.pattern}' --build {params.build} "
        "--min-replicates {params.min_rep} --out {output}"


rule count_peaks:
    input:
        bed=RES + "/{study}/atac/consensus/consensus.bed",
        bams=lambda wc: expand(
            RES + "/{study}/atac/bam/{sample}.bam", study=wc.study, sample=ATAC_BY_STUDY[wc.study]
        ),
        bais=lambda wc: expand(
            RES + "/{study}/atac/bam/{sample}.bam.bai", study=wc.study, sample=ATAC_BY_STUDY[wc.study]
        ),
    output:
        RES + "/{study}/atac/consensus/counts.tsv",
    params:
        samples=lambda wc: ATAC_BY_STUDY[wc.study],
    conda:
        "../envs/peaks.yaml"
    shell:
        "python {SCRIPTS}/py/count_peaks.py --bed {input.bed} --bams {input.bams} "
        "--samples {params.samples} --out {output}"


_deseq = config["differential"]["deseq2"]


rule deseq2_dar:
    input:
        counts=RES + "/{study}/atac/consensus/counts.tsv",
        sheet=RES + "/{study}/atac/group_sheet.tsv",
    output:
        RES + "/{study}/atac/dar/{contrast}.tsv",
    params:
        treatment=lambda wc: CONTRASTS[wc.contrast]["treatment"],
        control=lambda wc: CONTRASTS[wc.contrast]["control"],
        shrink=_deseq["lfc_shrink"],
        fdr=_deseq["fdr"],
        lfc=_deseq["min_abs_log2fc"],
        seed=config["project"]["seed"],
    log:
        RES + "/logs/{study}/{contrast}.deseq2.log",
    conda:
        "../envs/r_diff.yaml"
    shell:
        "Rscript {SCRIPTS}/R/deseq2_dar.R --counts {input.counts} --samples {input.sheet} "
        "--treatment {params.treatment} --control {params.control} --shrink {params.shrink} "
        "--fdr {params.fdr} --min_abs_log2fc {params.lfc} --seed {params.seed} "
        "--out {output} > {log} 2>&1"
