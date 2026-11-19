# Data acquisition rules. Only the synthetic toy dataset is implemented (PRD section 6);
# SRA/ENA downloads with provenance arrive with the real-data cards.
if config["reference"]["genome_build"] == "toy":
    TOY_DIR = str(Path(REF["fasta"]).parent)
    TOY_FASTQ = [p for r in ATAC if r.study_id == "toy" for p in fastq_paths(r.sample_id)]

    rule toy_data:
        output:
            fasta=REF["fasta"],
            gtf=REF["gtf"],
            chrom_sizes=REF["chrom_sizes"],
            blacklist=REF["blacklist"],
            ccres=REF["ccres"],
            manifest=TOY_DIR + "/MANIFEST.json",
            truth_regions=TOY_DIR + "/truth_regions.tsv",
            truth_genes=TOY_DIR + "/truth_genes.tsv",
            rna_counts=TOY_DIR + "/rna_counts.tsv",
            fastq=TOY_FASTQ,
        params:
            out=TOY_DIR,
            seed=config["project"]["seed"],
        log:
            RES + "/logs/toy_data.log",
        conda:
            "../envs/py_core.yaml"
        shell:
            "python tests/data/toy/make_toy.py --out {params.out} --seed {params.seed} > {log} 2>&1"
