rule unzip_geno_data:
    # Unzip .geno.gz output from ANGSD (hard-call genotypes)
    input:
        "data/genotypes/geno.geno.gz",
    output:
        "data/genotypes/geno.geno"
    shell:
        """
        gzip -d -k {input}
        """

rule prepare_data_for_structure:
    # Runs python script `angsd_to_structure` to transform ANGSD output into STRUCTURE input
    input:
        geno = "data/genotypes/geno.geno",
        meta = "config/meta_subset.csv"
    output:
        "data/structure/structure_input.tab"
    shell:
        """
        python scripts/angsd_to_structure.py -g {input.geno} -m {input.meta} -o {output}
        """
        
rule run_structure:
    # Runs STRUCTURE. Each k is designed to run in parallel.
    input:
        data = "data/structure/structure_input.tab",
        mainparams = "data/structure/mainparams.txt",
        extraparams = "data/structure/extraparams.txt"
    output:
        "data/structure/structure_output_{k}_f"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        structure -m {input.mainparams} -e {input.extraparams} -K {wildcards.k} -o data/structure/structure_output_{wildcards.k}
        """
        
rule harvest_structure:
    # Uses structureHarvester to generate summary with the likelihood of the data given k
    input:
        expand("data/structure/structure_output_{k}_f", k=range(1,11))
    output:
        "data/structure/summary.txt"
    shell:
        """
        structureHarvester.py --dir data/structure --out data/structure
        """
        
        
rule get_q_matrices:
    # Extracts the q-matrices from STRUCTURE's default output and stores into a csv for easy use later.
    # q-matrices show the proportion of each individual's ancestry to each of the k ancestral populations
    input:
        "data/structure/structure_output_{k}_f"
    output:
        "data/structure/q_matrix_{k}"
    shell:
        """
        python scripts/get_q_matrices.py -i {input} -o {output} -k {wildcards.k}
        """
    

rule get_delta_k_plot:
    # Uses the output from structureHarvester to generate k and delta k plots
    input:
        "data/structure/summary.txt"
    output:
        "figures/delta_k.png"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/get_delta_k_plot.py -i {input} -o {output}
        """
        
rule get_structure_barplot:
    # Generates structure barplots for 1 through k populations
    input:
        expand("data/structure/q_matrix_{k}", k =range(1,11))
    output:
        "figures/structure_barplot.png"
    conda:
        "../envs/analyses.yml"
    params:
        basename = "data/structure/q_matrix_",
        k = 4  
    shell:
        """
        python scripts/get_structure_plots.py -i {params.basename} -o {output} -k {params.k}
        """
        
rule transform_angsd_for_PCA:
    # Calls the script `angsd_to_pca.py` which is responsible for transforming ANGSD output data into a format usable for PCA
    # Encodes the homozygous for major allele as 0, heterozygous as 1, and homozygous minor as 2 (i.e. number of minor alleles)
    # Missing values are imputed with the average for that loci
    input:
        geno = "data/genotypes/geno.geno",
        meta = "config/meta_subset.csv"
    output:
        "data/genotypes/genotype_matrix.csv"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/angsd_to_pca.py -g {input.geno} -m {input.meta} -o {output} 
        """
    
    
rule get_PCA_plot:
    # Calls the script `get_pca_plot.py` to generate a PCA plot using sci-kit learn (2 dimensions)
    input:
        geno_matrix = "data/genotypes/genotype_matrix.csv",
        q_matrix = "data/structure/q_matrix_{k}"
    output:
        "figures/pca_plot_{k}.png"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/get_pca_plot.py -g {input.geno_matrix} -q {input.q_matrix} -o {output} -k {wildcards.k}
        """
        
rule get_pie_map:
    # Calls the script `get_pie_map` which plots each individual as a pie chart with their proportion of ancestry
    # from each of the k populations. Individuals are plotted with their latitude-longitude coordinates.
    input:
        q_matrix = "data/structure/q_matrix_{k}",
        meta = "config/meta_subset.csv"
    output:
        "figures/pie_map_{k}.png"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/get_pie_map.py -q {input.q_matrix} -m {input.meta} -o {output} -k {wildcards.k}
        """
        
rule get_combined_plots:
    # Calls the script `get_combined_plots.py` 
    # which produces a single figure with the STRUCTURE bar plot, pie map, and PCA plot for each k.
    input:
        q_matrix = "data/structure/q_matrix_{k}",
        meta = "config/meta_subset.csv",
        geno_matrix = "data/genotypes/genotype_matrix.csv"
    output:
        "figures/structure_PCA_pie_{k}.png"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/get_combined_plots.py -q {input.q_matrix} -g {input.geno_matrix} -m {input.meta} -o {output} -k {wildcards.k}
        """
    
rule get_heterozygosity_map:
    # Calls the script `get_heterozygosity_map.py` which maps heterozygosity (proportion of heterozygous loci)
    # across the geographic landscape. Heterozygosity densities are interpolated using pykrig (ordinary kriging, linear)
    input:
        meta = "config/meta_subset.csv",
        geno = "data/genotypes/geno.geno"
    output:
        "figures/heterozygosity_map.png"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/get_heterozygosity_map.py --g {input.geno} -m {input.meta} -o {output}
        """
        
rule prepare_fst_input:
    # Calls the script `prepare_fst_input.py` to prepare input for the R package hierfstat (used to calculate pairwise Fst)
    input:
        geno = "data/genotypes/geno.geno",
        q_matrix = "data/structure/q_matrix_{k}",
        meta = "config/meta_subset.csv"
    output:
        "data/genotypes/fst_input_{k}"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/prepare_fst_input.py --g {input.geno} -q {input.q_matrix} -m {input.meta} -o {output}
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
        
        
    