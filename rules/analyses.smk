rule unzip_geno_data:
    input:
        "data/genotypes/geno.geno.gz",
    output:
        "data/genotypes/geno.geno"
    shell:
        "gzip -d -k {input}"

rule prepare_data_for_structure:
    input:
        geno = "data/genotypes/geno.geno",
        meta = "config/meta.csv"
    output:
        "data/structure/structure_input.tab"
    shell:
        "python scripts/angsd_to_structure.py -g {input.geno} -m {input.meta} -o {output}"
        
rule run_structure:
    input:
        data = "data/structure/structure_input.tab",
        mainparams = "data/structure/mainparams.txt",
        extraparams = "data/structure/extraparams.txt"
    output:
        "run_structure_checkpoint.txt"
    conda:
        "../envs/analyses.yml"
    shell:
        "structure -m {input.mainparams} -e {input.extraparams} "
        
# rule tess3:
#     input:
#         "config/meta.csv",
#         "data/genotypes/geno.geno"
#     output:
#         "results/dummy_file",
#     conda:
#         "../envs/analyses.yml"
#     script:
#         "../scripts/tess3.R"

