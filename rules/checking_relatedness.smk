rule get_allele_frequency_file:
    # Creates file with allele frequencies (input for NgsRelate)
    input:
        "data/genotypes/geno.mafs.gz"
    output:
        "data/genotypes/freq"
    shell:
        """
        zcat {input} | cut -f5 |sed 1d > {output}
        """

rule run_NGS_relate:
    # Runs NgsRelate to get relatedness values for all pairs of samples
    # We do this to validate our technical replicates (should have high relatedness)
    # And to mark close relatives (rab > 0.45) for exclusion prior to running STRUCTURE
    # Within a highly related pair, the relative with the lower number of genomic sites is excluded
    input:
        glf = "data/genotypes/geno.glf.gz",
        freq = "data/genotypes/freq"
    output:
        "data/genotypes/ngs_relate_out"
    conda:
        "../envs/call_variants.yml"
    shell: 
        """
        ngsRelate  -g {input.glf} -f {input.freq}  -O {output} -r 0 -p 30 -c -e -n 339
        """
        
rule exclude_member_from_dyads:
    # Calls the script `scripts/get_meta_w_dyads.py` to create a new meta file meta_subset.csv that has a new column "to_exclude"
    # which marks close relatives that are to be excluded in STRUCTURE 
    input:
        ngs = "data/genotypes/ngs_relate_out",
        meta = "config/meta.csv",
        geno = "data/genotypes/geno.geno"
    output:
        "config/meta_subset.csv"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/get_meta_w_dyads.py -g {input.geno} -m {input.meta} -r {input.ngs} -o {output}
        """