rule call_genotypes_1:
    # Calls genotypes using ANGSD
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
    params:
        minInd = int(len(samples) * 0.5)
        
    conda:
        "../envs/call_variants.yml"
        
    shell:
        """
        echo "angsd:     $(angsd --version)" >> versions.txt
        angsd -b {input.list_of_bams} -gl 1 -domajorminor 1 -dosnpstat 1 -doHWE 1 -snp_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-5 -domaf 1 -doGlf 3 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -minInd {params.minInd} -setMinDepthInd 5 -setMaxDepthInd 1000 -minMaf 0.01 -doCounts 1 -doDepth 1 -dumpCounts 2 -doPost 1 -doGeno 4 -skipTriallelic 1 -out data/genotypes/geno -nThreads 35
        """
    
rule unzip_geno_1:
    # Unzip .geno.gz output from ANGSD (hard-call genotypes)
    input:
        "data/genotypes/geno.geno.gz"
    output:
        "data/genotypes/geno.geno"
    shell:
        "gzip -d -k {input}"   
    
rule get_missing_loci_column_1:
    input:
        geno = "data/genotypes/geno.geno",
        meta = "config/meta.csv"
    output:
        "tables/missing_loci_1.csv"
    shell:
        """
        python scripts/get_missing_loci_data.py -g {input.geno} -m {input.meta} -c missing_loci_1 -o {output}
        """
    
rule get_bams_wo_replicates_1:
    # Filter out low quality replicates and samples with > 50% missing genotypes
    input:
        geno = "data/genotypes/geno.geno",
        meta = "config/meta.csv"
    output:
        "data/merged_bams/bams_filtered_1.txt"
    shell:
        """
        python scripts/get_bams_wo_reps.py -g {input.geno} -m {input.meta} -o {output}
        """
        
        

rule call_genotypes_2:
    # Calls genotypes after the removal of lower quality replicates
    input:
        bam_index_files = expand("data/merged_bams/{sample}.bam.bai",sample=samples),
        bam_files =expand("data/merged_bams/{sample}.bam", sample=samples),
        list_of_bams = "data/merged_bams/bams_filtered_1.txt"
        
    output:
        "data/genotypes/geno2.arg",
        "data/genotypes/geno2.pos.gz",
        "data/genotypes/geno2.counts.gz",
        "data/genotypes/geno2.glf.gz",
        "data/genotypes/geno2.glf.pos.gz",
        "data/genotypes/geno2.mafs.gz",
        "data/genotypes/geno2.geno.gz",
        "data/genotypes/geno2.depthSample",
        "data/genotypes/geno2.depthGlobal",
        
    conda:
        "../envs/call_variants.yml"
        
    shell:
        """
        MININD=(`wc -l {input.list_of_bams}`)
        MININD=$((MININD*0.8))

        
        angsd -b {input.list_of_bams} -gl 1 -domajorminor 1 -dosnpstat 1 -doHWE 1 -snp_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-5 -domaf 1 -doGlf 3 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -minInd $MININD -setMinDepthInd 5 -setMaxDepthInd 1000 -minMaf 0.01 -doCounts 1 -doDepth 1 -doBcf 1 -dumpCounts 2 -doPost 1 -doGeno 4 -skipTriallelic 1 -out data/genotypes/geno2 -nThreads 35
        """
        
rule unzip_geno_2:
    # Unzip .geno2.gz output from ANGSD (hard-call genotypes)
    input:
        "data/genotypes/geno2.geno.gz"
    output:
        "data/genotypes/geno2.geno"
    shell:
        "gzip -d -k {input}"   
    
###
# REMOVE INDIVIDUALS WITH > 25% data
###

###
# CALL GENOTYPES AGAIN
###
    
    
rule get_missing_loci_column_again: # may be redundant
    input:
        geno = "data/genotypes/geno2.geno",
        meta = "config/meta_subset.csv"
    output:
        "tables/missing_loci_2.csv"
    shell:
        """
        python scripts/get_missing_loci_data.py -g {input.geno} -m {input.meta} -c missing_loci_2 -o {output}
        """
        
rule get_allele_frequency_file:
    # Creates file with allele frequencies (input for NgsRelate)
    input:
        "data/genotypes/geno2.mafs.gz"
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
        glf = "data/genotypes/geno2.glf.gz",
        freq = "data/genotypes/freq",
        bamlist = 'data/merged_bams/bams_filtered.txt'
    output:
        "data/genotypes/ngs_relate_out"
    conda:
        "../envs/call_variants.yml"
    shell: 
        """
        echo "ngsRelate:     $(ngsRelate --version)" >> versions.txt
        N=(`wc -l < {input.bamlist}`)
        ngsRelate  -g {input.glf} -f {input.freq}  -O {output} -r 0 -p 30 -c -e -n $N -z {input.bamlist}
        """
        

        
rule exclude_member_from_dyads:
    # Calls the script `scripts/get_meta_w_dyads.py` to create a new meta file meta_subset.csv that has a new column "to_exclude"
    # which marks close relatives that are to be excluded in STRUCTURE 
    input:
        ngs = "data/genotypes/ngs_relate_out",
        meta = "config/meta.csv",
        geno = "data/genotypes/geno2.geno",
        bamlist = 'data/merged_bams/bams_filtered.txt'
    output:
        "config/meta_subset.csv"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/get_meta_w_dyads.py -g {input.geno} -m {input.meta} -r {input.ngs} -b {input.bamlist} -o {output}
        """