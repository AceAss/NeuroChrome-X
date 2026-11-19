# Coverage tracks: CPM bigWig per sample and per-condition mean (PRD section 7).
_bc = config["tracks"]["bamcoverage"]

if config["tracks"].get("enabled", True):

    rule bigwig_sample:
        input:
            bam=BAM + ".bam",
            bai=BAM + ".bam.bai",
        output:
            RES + "/{study}/atac/tracks/{sample}.bw",
        params:
            norm=_bc["normalize_using"],
            binsize=_bc["bin_size"],
            fragments=lambda wc: "--extendReads" if wc.sample in PE_SAMPLES else "",
        threads: config["compute"]["threads_default"]
        log:
            RES + "/logs/{study}/{sample}.bamcoverage.log",
        conda:
            "../envs/align_dna.yaml"
        shell:
            "bamCoverage -b {input.bam} -o {output} --normalizeUsing {params.norm} "
            "--binSize {params.binsize} {params.fragments} -p {threads} > {log} 2>&1"

    rule bigwig_condition_mean:
        input:
            lambda wc: [
                f"{RES}/{wc.study}/atac/tracks/{s}.bw" for s in CONDITIONS[(wc.study, wc.condition)]
            ],
        output:
            RES + "/{study}/atac/tracks/conditions/{condition}.mean.bw",
        params:
            binsize=_bc["bin_size"],
        threads: config["compute"]["threads_default"]
        log:
            RES + "/logs/{study}/{condition}.bigwigaverage.log",
        conda:
            "../envs/align_dna.yaml"
        shell:
            "bigwigAverage -b {input} -o {output} -bs {params.binsize} -p {threads} > {log} 2>&1"
