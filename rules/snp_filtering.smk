# geno > 50%
# imiss < 90%
# geno > 60%
# imiss < 70%
# geno > 70%
# imiss < 50%


rule get_list_of_bams_0:
    # Makes a list of bam file paths that will be used as input for ANGSD (all samples)
    input:
        bam_index_files = expand("data/merged_bams/{sample}.bam.bai",sample=samples),
        bam_files =expand("data/merged_bams/{sample}.bam", sample=samples)
        
    output:
        "data/merged_bams/bams0.txt"
    conda:
        "../envs/preprocessing.yml"
    shell:
        """
        ls data/merged_bams/*.bam > data/merged_bams/bams0.txt
        """
        

#######################################################################################################
##############################                   ROUND 1                 ##############################
#######################################################################################################
        
rule call_genotypes_1:
    input:
        bam_index_files = expand("data/merged_bams/{sample}.bam.bai",sample=samples),
        bam_files =expand("data/merged_bams/{sample}.bam", sample=samples),
        list_of_bams = "data/merged_bams/bams0.txt"
        
    output:
        "data/genotypes/geno1.arg",
        "data/genotypes/geno1.pos.gz",
        "data/genotypes/geno1.counts.gz",
        "data/genotypes/geno1.glf.gz",
        "data/genotypes/geno1.glf.pos.gz",
        "data/genotypes/geno1.mafs.gz",
        "data/genotypes/geno1.geno.gz",
        "data/genotypes/geno1.depthSample",
        "data/genotypes/geno1.depthGlobal",
        
    conda:
        "../envs/call_variants.yml"
        
    shell:
        """
        echo "angsd:     $(angsd --version)" >> versions.txt
        
        MININD=(`wc -l {input.list_of_bams}`)
        MININD=$(echo "$MININD * 0.5" | bc)

        angsd -b {input.list_of_bams} -gl 1 -domajorminor 1 -dosnpstat 1 -doHWE 1 -snp_pval 1e-5 -domaf 1 -doGlf 3 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -minInd $MININD -setMinDepthInd 5 -minMaf 0.01 -doCounts 1 -doDepth 1 -doBcf 1 -dumpCounts 2 -doPost 1 -doGeno 4 -skipTriallelic 1 -out data/genotypes/geno1 -nThreads 35
        """
        
rule unzip_geno_1:
    input:
        "data/genotypes/geno1.geno.gz"
    output:
        "data/genotypes/geno1.geno"
    shell:
        """
        gzip -d -k {input}
        """
        
rule apply_imiss_filter1:
    input:
        geno = "data/genotypes/geno1.geno",
        prev_bam_list = "data/merged_bams/bams0.txt",

    output:
        filtered_bams = "data/merged_bams/bams1.txt",
        log = "data/genotypes/imiss_filter1_log.txt"
    params:
        thresh = 0.9
    shell:
        """
        python scripts/apply_imiss_filter.py -g {input.geno} -b {input.prev_bam_list} -f {params.thresh} -o {output.filtered_bams} &> {output.log}
        """
        
#######################################################################################################
##############################                   ROUND 2                 ##############################
#######################################################################################################        
        
rule call_genotypes_2:
    input:
        bam_index_files = expand("data/merged_bams/{sample}.bam.bai",sample=samples),
        bam_files =expand("data/merged_bams/{sample}.bam", sample=samples),
        list_of_bams = "data/merged_bams/bams1.txt"
        
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
        echo "angsd:     $(angsd --version)" >> versions.txt
        
        MININD=(`wc -l {input.list_of_bams}`)
        MININD=$(echo "$MININD * 0.6" | bc)

        angsd -b {input.list_of_bams} -gl 1 -domajorminor 1 -dosnpstat 1 -doHWE 1 -snp_pval 1e-5 -domaf 1 -doGlf 3 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -minInd $MININD -setMinDepthInd 5 -minMaf 0.01 -doCounts 1 -doDepth 1 -doBcf 1 -dumpCounts 2 -doPost 1 -doGeno 4 -skipTriallelic 1 -out data/genotypes/geno2 -nThreads 35
        """
        
rule unzip_geno_2:
    input:
        "data/genotypes/geno2.geno.gz"
    output:
        "data/genotypes/geno2.geno"
    shell:
        """
        gzip -d -k {input}
        """
        
rule apply_imiss_filter2:
    input:
        geno = "data/genotypes/geno2.geno",
        prev_bam_list = "data/merged_bams/bams1.txt",
    output:
        filtered_bams = "data/merged_bams/bams2.txt",
        log = "data/genotypes/imiss_filter2_log.txt"
    params:
        thresh = 0.7
    shell:
        """
       python scripts/apply_imiss_filter.py -g {input.geno} -b {input.prev_bam_list} -f {params.thresh} -o {output.filtered_bams} &> {output.log}
        """

#######################################################################################################
##############################                   ROUND 3                 ##############################
#######################################################################################################        
        
rule call_genotypes_3:
    input:
        bam_index_files = expand("data/merged_bams/{sample}.bam.bai",sample=samples),
        bam_files =expand("data/merged_bams/{sample}.bam", sample=samples),
        list_of_bams = "data/merged_bams/bams2.txt"
        
    output:
        "data/genotypes/geno3.arg",
        "data/genotypes/geno3.pos.gz",
        "data/genotypes/geno3.counts.gz",
        "data/genotypes/geno3.glf.gz",
        "data/genotypes/geno3.glf.pos.gz",
        "data/genotypes/geno3.mafs.gz",
        "data/genotypes/geno3.geno.gz",
        "data/genotypes/geno3.depthSample",
        "data/genotypes/geno3.depthGlobal",
        
    conda:
        "../envs/call_variants.yml"
        
    shell:
        """
        echo "angsd:     $(angsd --version)" >> versions.txt
        
        MININD=(`wc -l {input.list_of_bams}`)
        MININD=$(echo "$MININD * 0.7" | bc)

        angsd -b {input.list_of_bams} -gl 1 -domajorminor 1 -dosnpstat 1 -doHWE 1 -snp_pval 1e-5 -domaf 1 -doGlf 3 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -minInd $MININD -setMinDepthInd 5 -minMaf 0.01 -doCounts 1 -doDepth 1 -doBcf 1 -dumpCounts 2 -doPost 1 -doGeno 4 -skipTriallelic 1 -out data/genotypes/geno3 -nThreads 35
        """
        
rule unzip_geno_3:
    input:
        "data/genotypes/geno3.geno.gz"
    output:
        "data/genotypes/geno3.geno"
    shell:
        """
        gzip -d -k {input}
        """
        
rule apply_imiss_filter3:
    input:
        geno = "data/genotypes/geno3.geno",
        prev_bam_list = "data/merged_bams/bams2.txt",
    output:
        filtered_bams = "data/merged_bams/bams3.txt",
        log = "data/genotypes/imiss_filter3_log.txt"
    params:
        thresh = 0.5
    shell:
        """
        python scripts/apply_imiss_filter.py -g {input.geno} -b {input.prev_bam_list} -f {params.thresh} -o {output.filtered_bams} &> {output.log}
        """
        
        
rule call_genotypes_4:
    input:
        bam_index_files = expand("data/merged_bams/{sample}.bam.bai",sample=samples),
        bam_files =expand("data/merged_bams/{sample}.bam", sample=samples),
        list_of_bams = "data/merged_bams/bams3.txt"
        
    output:
        "data/genotypes/geno4.arg",
        "data/genotypes/geno4.pos.gz",
        "data/genotypes/geno4.counts.gz",
        "data/genotypes/geno4.glf.gz",
        "data/genotypes/geno4.glf.pos.gz",
        "data/genotypes/geno4.mafs.gz",
        "data/genotypes/geno4.geno.gz",
        "data/genotypes/geno4.depthSample",
        "data/genotypes/geno4.depthGlobal",
        
    conda:
        "../envs/call_variants.yml"
        
    shell:
        """
        echo "angsd:     $(angsd --version)" >> versions.txt
        
        MININD=(`wc -l {input.list_of_bams}`)
        MININD=$(echo "$MININD * 0.7" | bc)

        angsd -b {input.list_of_bams} -gl 1 -domajorminor 1 -dosnpstat 1 -doHWE 1 -snp_pval 1e-5 -domaf 1 -doGlf 3 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -minInd $MININD -setMinDepthInd 5 -minMaf 0.01 -doCounts 1 -doDepth 1 -doBcf 1 -dumpCounts 2 -doPost 1 -doGeno 4 -skipTriallelic 1 -out data/genotypes/geno4 -nThreads 35
        """
        
rule unzip_geno_4:
    input:
        "data/genotypes/geno4.geno.gz"
    output:
        "data/genotypes/geno4.geno"
    shell:
        """
        gzip -d -k {input}
        """
        
#######################################################################################################
##############################           REMOVE CLOSE RELATIVES          ##############################
#######################################################################################################   

rule get_allele_frequency_file:
    # Creates file with allele frequencies (input for NgsRelate)
    input:
        "data/genotypes/geno4.mafs.gz"
    output:
        "data/genotypes/freq4"
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
        glf = "data/genotypes/geno4.glf.gz",
        freq = "data/genotypes/freq4",
        bamlist = 'data/merged_bams/bams3.txt'
    output:
        "data/genotypes/ngs_relate_rcr"
    conda:
        "../envs/call_variants.yml"
    shell: 
        """
        N=(`wc -l < {input.bamlist}`)
        ngsRelate  -g {input.glf} -f {input.freq}  -O {output} -r 0 -p 30 -c -e -n $N -z {input.bamlist}
        """
        
rule filter_close_relatives_and_replicates:
    input:
        geno = "data/genotypes/geno4.geno",
        prev_bam_list = "data/merged_bams/bams3.txt",
        relate = "data/genotypes/ngs_relate_rcr"
    output:
        close_relatives = "tables/close_relatives.tab",
        bam_list = "data/merged_bams/bams_rcr.txt"
    params:
        rab_cutoff = 0.45
    shell:
        """
        python scripts/apply_rcr_filter.py -g {input.geno} -b {input.prev_bam_list} -r {input.relate} -c {params.rab_cutoff} -s {output.close_relatives} -o {output.bam_list}
        """
        
        
rule apply_imiss_filter4:
    input:
        geno = "data/genotypes/geno4.geno",
        prev_bam_list = "data/merged_bams/bams3.txt",
    output:
        filtered_bams = "data/merged_bams/bams4.txt",
        log = "data/genotypes/imiss_filter4_log.txt"
    params:
        thresh = 0.25
    shell:
        """
        python scripts/apply_imiss_filter.py -g {input.geno} -b {input.prev_bam_list} -f {params.thresh} -o {output.filtered_bams} &> {output.log}
        """
    
rule combine_4_rcr_filters:
    input:
        four = "data/merged_bams/bams4.txt",
        rcr = "data/merged_bams/bams_rcr.txt"
    output:
        "data/merged_bams/bams_4_rcr.txt"
    shell:
        """
        comm -1 -2 <(sort {input.four}) <(sort {input.rcr}) > {output}
        """
        
#######################################################################################################
##############################           FINAL SNP CALLING          ##############################
#######################################################################################################  
        
rule call_genotypes_5:
    input:
        bam_index_files = expand("data/merged_bams/{sample}.bam.bai",sample=samples),
        bam_files =expand("data/merged_bams/{sample}.bam", sample=samples),
        list_of_bams = "data/merged_bams/bams_4_rcr.txt"
        
    output:
        "data/genotypes/geno5.arg",
        "data/genotypes/geno5.pos.gz",
        "data/genotypes/geno5.counts.gz",
        "data/genotypes/geno5.glf.gz",
        "data/genotypes/geno5.glf.pos.gz",
        "data/genotypes/geno5.mafs.gz",
        "data/genotypes/geno5.geno.gz",
        "data/genotypes/geno5.depthSample",
        "data/genotypes/geno5.depthGlobal",
        
    conda:
        "../envs/call_variants.yml"
        
    shell:
        """
        MININD=(`wc -l {input.list_of_bams}`)
        MININD=$(echo "$MININD * 0.8" | bc)

        angsd -b {input.list_of_bams} -gl 1 -domajorminor 1 -dosnpstat 1 -doHWE 1 -snp_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-5 -domaf 1 -doGlf 3 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -minInd $MININD -setMinDepthInd 5 -minMaf 0.01 -doCounts 1 -doDepth 1 -doBcf 1 -dumpCounts 2 -doPost 1 -doGeno 4 -skipTriallelic 1 -out data/genotypes/geno5 -nThreads 35
        """
        
rule unzip_geno_5:
    input:
        "data/genotypes/geno5.geno.gz"
    output:
        "data/genotypes/geno5.geno"
    shell:
        """
        gzip -d -k {input}
        """        
        
rule unzip_geno_5_counts:
    input:
        "data/genotypes/geno5.counts.gz"
    output:
        "data/genotypes/geno5.counts"
    shell:
        """
        gzip -d -k {input}
        """        
