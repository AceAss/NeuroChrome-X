# Element annotation, expression filter and motif enrichment (PRD section 8; ADR-007).
_gene_filter = config["rna"]["gene_filter"]


def _expressed_input(wc):
    info = rna_inputs(wc.study)
    return [info[0]] if info else []


rule expressed_genes:
    input:
        counts=_expressed_input,
    output:
        RES + "/{study}/expressed_genes.txt",
    params:
        samples=lambda wc: (rna_inputs(wc.study) or ("", []))[1],
        min_count=_gene_filter["min_count"],
        min_samples=_gene_filter["min_samples"],
    conda:
        "../envs/py_core.yaml"
    shell:
        "python {SCRIPTS}/py/expressed_genes.py --counts {input.counts} "
        "--samples {params.samples} --min-count {params.min_count} "
        "--min-samples {params.min_samples} --out {output}"


_ann = config["annotation"]


def _annotate_inputs(wc):
    inputs = {
        "elements": f"{RES}/{wc.study}/atac/consensus/consensus.bed",
        "gtf": REF["gtf"],
        "ccres": REF["ccres"],
    }
    if rna_inputs(wc.study):
        inputs["expressed"] = f"{RES}/{wc.study}/expressed_genes.txt"
    return inputs


rule annotate_elements:
    input:
        unpack(_annotate_inputs),
    output:
        RES + "/{study}/atac/elements/regulatory_elements.tsv",
    params:
        expressed=lambda wc, input: f"--expressed {input.expressed}" if "expressed" in input.keys() else "",
        up=_ann["promoter"]["upstream"],
        down=_ann["promoter"]["downstream"],
        classes=" ".join(_ann["enhancer"]["tier1_ccre_classes"]),
        tier3=_ann["enhancer"]["tier3_min_tss_distance"],
        window=_ann["gene_link"]["window"],
    conda:
        "../envs/py_core.yaml"
    shell:
        "python {SCRIPTS}/py/annotate_elements.py --elements {input.elements} --gtf {input.gtf} "
        "--ccres {input.ccres} {params.expressed} --promoter-upstream {params.up} "
        "--promoter-downstream {params.down} --ccre-classes {params.classes} "
        "--tier3-min-tss-distance {params.tier3} --link-window {params.window} --out {output}"


_mot = config["motifs"]
_py = _mot.get("python_engine", {})

if _mot.get("engine", "python") == "python":

    rule motif_enrichment:
        input:
            dar=RES + "/{study}/atac/dar/{contrast}.tsv",
            fasta=REF["fasta"],
        output:
            RES + "/{study}/atac/motifs/{contrast}.{direction}.tsv",
        wildcard_constraints:
            direction="up|down",
        params:
            motifs=motif_args(),
            min_rel=_py.get("min_rel_score", 0.9),
            ratio=_py.get("ratio", 2),
            min_padj=_mot["background"].get("min_padj", 0.2),
            seed=config["project"]["seed"],
        conda:
            "../envs/py_core.yaml"
        shell:
            "python {SCRIPTS}/py/motif_enrichment.py --dar {input.dar} --fasta {input.fasta} "
            "--direction {wildcards.direction} {params.motifs} --min-rel-score {params.min_rel} "
            "--ratio {params.ratio} --background-min-padj {params.min_padj} "
            "--seed {params.seed} --out {output}"
