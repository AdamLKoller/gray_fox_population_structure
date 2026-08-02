rule angsd_to_lfmm_format:
    # Calls the script `angsd_to_lfmm.py` to convert the output of ANGSD into .lfmm format so that it is readable by TESS3r
    input:
        geno = "data/genotypes/geno5.geno",
    output:
        "data/genotypes/geno5.lfmm"
    shell:
        """
        python scripts/angsd_to_lfmm.py -g {input.geno} -o {output}
        """
        
rule run_tess3r:
    input:
        lfmm = "data/genotypes/geno5.lfmm",
        meta = "config/meta_apriori_jittered.csv",
        bam_list = "data/merged_bams/bams_4_rcr.txt"
    output:
        "figures/tess3_cross_validation_plot.png"
    conda:
        "../envs/tess3.yml"
    shell:
        """
        Rscript scripts/run_tess3r.R -l {input.lfmm} -m {input.meta} -o figures/tess3 -b {input.bam_list}
        """