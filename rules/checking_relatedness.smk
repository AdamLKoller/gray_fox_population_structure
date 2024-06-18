rule get_allele_frequency_file:
    input:
        "data/genotypes/geno.mafs.gz"
    output:
        "data/genotypes/freq"
    shell:
        "zcat {input} | cut -f5 |sed 1d > {output}"

rule run_NGS_relate:
    input:
        glf = "data/genotypes/geno.glf.gz",
        freq = "data/genotypes/freq"
    output:
        "data/genotypes/ngs_relate_out"
    conda:
        "../envs/call_variants.yml"
    shell: 
        "ngsRelate  -g {input.glf} -f {input.freq}  -O {output} -r 0 -p 30 -c -e -n 339"
        
rule exclude_member_from_dyads:
    input:
        ngs = "data/genotypes/ngs_relate_out",
        meta = "config/meta.csv",
        geno = "data/genotypes/geno.geno"
    output:
        "config/meta_subset.csv"
    conda:
        "../envs/analyses.yml"
    shell:
        '''
        python scripts/get_meta_w_dyads.py -g {input.geno} -m {input.meta} -r {input.ngs} -o {output}
        '''