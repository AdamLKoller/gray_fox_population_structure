rule prepare_fst_input:
    # Calls the script `prepare_fst_input.py` to prepare input for the R package hierfstat (used to calculate #pairwise Fst)
    input:
        geno = "data/genotypes/geno5.geno",
        q_matrix = "data/structure/q_matrix_{k}"
    output:
        "data/genotypes/fst_input_{k}"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/prepare_fst_input.py --g {input.geno} -q {input.q_matrix} -o {output}
        """
        
rule get_pairwise_fst:
    # Calls the script `get_pairwise_fst.R` to generate pairwise_fst table (calculated according to Nei 1973)
    input:
        "data/genotypes/fst_input_{k}"
    output:
        "tables/pairwise_fst_{k}"
    conda:
        '../envs/fstats.yml'
    shell:
        """
        Rscript scripts/get_pairwise_fst.R -i {input} -o {output}
        """
        
rule prepare_fst_input_apriori:
    # Calls the script `prepare_fst_input.py` to prepare input for the R package hierfstat (used to calculate pairwise Fst)
    input:
        geno = "data/genotypes/geno5.geno",
        meta = "config/meta_apriori.csv",
        bam_list = "data/merged_bams/bams_4_rcr.txt"
        
    output:
        "data/genotypes/fst_input_apriori"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/prepare_fst_input_apriori.py --g {input.geno} -m {input.meta} -o {output} -b {input.bam_list}
        """
        
rule get_pairwise_fst_apriori:
    # Calls the script `get_pairwise_fst.R` to generate pairwise_fst table (calculated according to Nei 1973)
    input:
        "data/genotypes/fst_input_apriori"
    output:
        "tables/pairwise_fst_apriori"
    conda:
        '../envs/fstats.yml'
    shell:
        """
        Rscript scripts/get_pairwise_fst.R -i {input} -o {output}
        """
        
rule get_allele_stats:
    # Calculates basic population genetics statistics for each of the clusters identified by STRUCTURE with hierfstat
    input:
        "data/genotypes/fst_input_{k}"
    output:
        "tables/allele_stats_{k}.tab"
    conda:
        '../envs/fstats.yml'
    shell:
        """
        Rscript scripts/get_pop_allele_stats.R -i {input} -o {output}
        """
        
rule get_allele_stats_apriori:
    # Calculates basic population genetics statistics for each of the apriori sampling groups with hierfstat
    input:
        "data/genotypes/fst_input_apriori"
    output:
        "tables/allele_stats_apriori.tab"
    conda:
        '../envs/fstats.yml'
    shell:
        """
        Rscript scripts/get_pop_allele_stats.R -i {input} -o {output}
        """