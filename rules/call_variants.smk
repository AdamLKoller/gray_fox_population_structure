rule call_genotypes:
    input:
        bam_index_files = expand("data/merged_bams/{sample}.bam.bai",sample=samples),
        bam_files =expand("data/merged_bams/{sample}.bam", sample=samples),
        list_of_bams = "data/merged_bams/bams.txt"
        
    output:
        "data/genotypes/geno.arg",
        "data/genotypes/geno.pos.gz",
        "data/genotypes/geno.counts.gz",
        "data/genotypes/geno.glf.gz",
        "data/genotypes/geno.glf.pos.gz",
        "data/genotypes/geno.mafs.gz",
        "data/genotypes/geno.geno.gz",
        "data/genotypes/geno.depthSample",
        "data/genotypes/geno.depthGlobal",
        
    conda:
        "../envs/call_variants.yml"
        
    shell:
        '''
        angsd -b {input.list_of_bams} -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -minInd 254 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 2 -doPost 1 -doGeno 4 -postCutoff 0.6 -out data/genotypes/geno
    '''
    