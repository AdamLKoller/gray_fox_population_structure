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
        meta = "config/meta_subset.csv"
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
        "data/structure/structure_output_{k}_f"
    conda:
        "../envs/analyses.yml"
    shell:
        "structure -m {input.mainparams} -e {input.extraparams} -K {wildcards.k} -o data/structure/structure_output_{wildcards.k}"
        
rule harvest_structure:
    input:
        expand("data/structure/structure_output_{k}_f", k=range(1,11))
    output:
        "data/structure/summary.txt"
    shell:
        "structureHarvester.py --dir data/structure --out data/structure"
        
        
rule get_q_matrices:
    input:
        "data/structure/structure_output_{k}_f"
    output:
        "data/structure/q_matrix_{k}"
    shell:
        "python scripts/get_q_matrices.py -i {input} -o {output} -k {wildcards.k}"
    

rule get_delta_k_plot:
    input:
        "data/structure/summary.txt"
    output:
        "figures/delta_k.png"
    conda:
        "../envs/analyses.yml"
    shell:
        "python scripts/get_delta_k_plot.py -i {input} -o {output}"
        
rule get_structure_barplot:
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
        "python scripts/get_structure_plots.py -i {params.basename} -o {output} -k {params.k}"
        

rule transform_angsd_for_PCA:
    input:
        geno = "data/genotypes/geno.geno",
        meta = "config/meta_subset.csv"
    output:
        "data/genotypes/genotype_matrix.csv"
    conda:
        "../envs/analyses.yml"
    shell:
        "python scripts/angsd_to_pca.py -g {input.geno} -m {input.meta} -o {output} "
    
    
rule get_PCA_plot:
    input:
        geno_matrix = "data/genotypes/genotype_matrix.csv",
        q_matrix = "data/structure/q_matrix_{k}"
    output:
        "figures/pca_plot_{k}.png"
    conda:
        "../envs/analyses.yml"
    shell:
        "python scripts/get_pca_plot.py -g {input.geno_matrix} -q {input.q_matrix} -o {output} -k {wildcards.k}"
        
rule get_pie_map:
    input:
        q_matrix = "data/structure/q_matrix_{k}",
        meta = "config/meta_subset.csv"
    output:
        "figures/pie_map_{k}.png"
    conda:
        "../envs/analyses.yml"
    shell:
        "python scripts/get_pie_map.py -q {input.q_matrix} -m {input.meta} -o {output} -k {wildcards.k}"
        
rule get_combined_plots:
    input:
        q_matrix = "data/structure/q_matrix_{k}",
        meta = "config/meta_subset.csv",
        geno_matrix = "data/genotypes/genotype_matrix.csv"
    output:
        "figures/structure_PCA_pie_{k}.png"
    conda:
        "../envs/analyses.yml"
    shell:
        "python scripts/get_combined_plots.py -q {input.q_matrix} -g {input.geno_matrix} -m {input.meta} -o {output} -k {wildcards.k}"
    
rule get_heterozygosity_map:
    input:
        meta = "config/meta_subset.csv",
        geno = "data/genotypes/geno.geno"
    output:
        "figures/heterozygosity_map.png"
    conda:
        "../envs/analyses.yml"
    shell:
        "python scripts/get_heterozygosity_map.py --g {input.geno} -m {input.meta} -o {output}"
        
    