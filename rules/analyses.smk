import shutil

rule prepare_data_for_structure:
    # Runs python script `angsd_to_structure` to transform ANGSD output into STRUCTURE input
    input:
        geno = "data/genotypes/geno2.geno",
        meta = "config/meta_subset.csv"
    output:
        "data/structure/structure_input.tab"
    shell:
        """
        python scripts/angsd_to_structure.py -g {input.geno} -m {input.meta} -o {output}
        """


rule create_structure_parameter_files:
    # Creates the parameter files used by STRUCTURE
    input:
        freq = "data/genotypes/freq",
        meta = 'config/meta_subset.csv'
    output:
        mainparams = "data/structure/mainparams.txt",
        extraparams = "data/structure/extraparams.txt"
        
    run:
        # Get the number of loci using shell command
        numloci = int(shell("wc -l < {input.freq}", read=True).strip())
        numinds = len(pd.read_csv('config/meta_subset.csv').loc[pd.read_csv('config/meta_subset.csv').to_exclude ==False])

        # Write mainparams file
        with open(output.mainparams, 'w') as mainparams:
            mainparams.write("#define MAXPOPS 10\n")
            mainparams.write("#define BURNIN 10000\n")
            mainparams.write("#define NUMREPS 25000\n")
            mainparams.write("#define INFILE data/structure/structure_input.tab\n")
            mainparams.write("#define OUTFILE data/structure/structure_output\n")
            mainparams.write(f"#define NUMINDS {numinds}\n")
            mainparams.write(f"#define NUMLOCI {numloci}\n")
            mainparams.write("#define PLOIDY 2\n")
            mainparams.write("#define MISSING -1\n")
            mainparams.write("#define ONEROWPERIND 0\n")
            mainparams.write("#define LABEL 0\n")
            mainparams.write("#define POPDATA 0\n")
            mainparams.write("#define POPFLAG 0\n")
            mainparams.write("#define LOCDATA 0\n")
            mainparams.write("#define PHENOTYPE 0\n")
            mainparams.write("#define EXTRACOLS 0\n")
            mainparams.write("#define MARKERNAMES 1\n")
            mainparams.write("#define RECESSIVEALLELES 0\n")
            mainparams.write("#define MAPDISTANCES 0\n")

        # Write extraparams file
        with open(output.extraparams, 'w') as extraparams:
            extraparams.write("#define RANDOMIZE 0\n")
            extraparams.write("#define NOADMIX 0\n")
        

        
        
rule run_structure:
    # Runs STRUCTURE. Each k is designed to run in parallel.
    input:
        data = "data/structure/structure_input.tab",
        mainparams = "data/structure/mainparams.txt",
        extraparams = "data/structure/extraparams.txt"
    output:
        "data/structure/structure_output_k{k}_rep{rep}_f"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        structure -m {input.mainparams} -e {input.extraparams} -K {wildcards.k} -o data/structure/structure_output_k{wildcards.k}_rep{wildcards.rep} -D {wildcards.rep}
        """
        
rule harvest_structure:
    # Uses structureHarvester to generate summary with the likelihood of the data given k
    input:
        expand("data/structure/structure_output_k{k}_rep{rep}_f", k=range(1,11), rep=range(1,4))
    output:
        "data/structure/summary.txt",
        "data/structure/evanno.txt"
    shell:
        """
        structureHarvester.py --dir data/structure --out data/structure --evanno
        """
        
rule pick_most_likely_reps:
    # Copies the most likely rep for each value of K and puts into a new file `structure_output_{k}_best`
    input:
        expand("data/structure/structure_output_k{k}_rep{rep}_f", k=range(1,11), rep=range(1,4))
    output:
        expand("data/structure/structure_output_{k}_best", k=range(1,11))
    run:
        k_to_best_rep = {str(k):(None,float('-inf')) for k in range(1,11)}
        for filename in os.listdir('data/structure'):
            if filename.endswith('_f'):
                k, rep = filename.split('/')[-1].split('_')[2][1:], filename.split('/')[-1].split('_')[-1]
                file = open('data/structure/'+filename, 'r')
                line = file.readline()
                while 'Estimated Ln Prob of Data' not in line:
                    line = file.readline()
                likelihood = float(line.split(' = ')[-1])
                if likelihood > k_to_best_rep[k][1]:
                    k_to_best_rep[k] = ('data/structure/'+filename, likelihood)
                    
        for k, best_rep in k_to_best_rep.items():
            print(f'The most likely rep for k={k} is {best_rep[0]} \t likelihood: {best_rep[1]}')
            shutil.copyfile(best_rep[0], f"data/structure/structure_output_{k}_best".replace(' ',''))
                
       
    
        
        
rule get_q_matrices:
    # Extracts the q-matrices from STRUCTURE's default output and stores into a csv for easy use later.
    # q-matrices show the proportion of each individual's ancestry to each of the k ancestral populations
    input:
        "data/structure/structure_output_{k}_best"
    output:
        "data/structure/q_matrix_{k}"
    shell:
        """
        python scripts/get_q_matrices.py -i {input} -o {output} -k {wildcards.k}
        """

rule get_delta_k_plot:
    # Uses the output from structureHarvester to generate k and delta k plots
    input:
        "data/structure/evanno.txt"
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
        geno = "data/genotypes/geno2.geno",
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
        

rule jitter_sample_locations:
    # Calls the script `jitter_locations.py` to jitter samples with the same location.
    # For some samples, just county or state data was collected. These points were given
    # longitude and latitude values equal to the centroid of their county or state.
    # We want to apply a jitter to these points so they do not overlap in visualizations.
    # The jittered locations are also used in TESS3
    input:
        meta = 'config/meta_subset.csv'
    output:
        'config/meta_subset_jittered.csv'
    shell:
        """
        python scripts/jitter_locations.py -m {input} -o {output}
        """
        
rule add_apriori_pops:
    # Calls the script `add_apriori_pops.py` to add another column 'POP'
    # containing apriori sample group to the meta file as defined in 
    # shapefiles in data/shapefiles/XYRed2017POP.*
    input:
        meta = 'config/meta_subset.csv'
    output:
        'config/meta_subset_apriori.csv'
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/add_apriori_pops.py -m {input} -o {output}
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
    
rule get_combined_plots_w_apriori:
    # Calls the script `get_combined_plots.py` 
    # which produces a single figure with the STRUCTURE bar plot, pie map, and PCA plot for each k.
    input:
        q_matrix = "data/structure/q_matrix_{k}",
        meta = "config/meta_subset_apriori.csv",
        geno_matrix = "data/genotypes/genotype_matrix.csv"
    output:
        "figures/structure_PCA_pie_apriori_{k}.png"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/get_combined_plots_apriori.py -q {input.q_matrix} -g {input.geno_matrix} -m {input.meta} -o {output} -k {wildcards.k}
        """
    
rule get_heterozygosity_map:
    # Calls the script `get_heterozygosity_map.py` which maps heterozygosity (proportion of heterozygous loci)
    # across the geographic landscape. Heterozygosity densities are interpolated using pykrig (ordinary kriging, linear)
    input:
        meta = "config/meta_subset.csv",
        geno = "data/genotypes/geno2.geno"
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
        geno = "data/genotypes/geno2.geno",
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
        
rule angsd_to_lfmm_format:
    # Calls the script `angsd_to_lfmm.py` to convert the output of ANGSD into .lfmm format so that it is readable by pcadapt.
    input:
        geno = "data/genotypes/geno2.geno",
        meta = "config/meta_subset.csv"
    output:
        "data/genotypes/geno2.lfmm"
    shell:
        """
        python scripts/angsd_to_lfmm.py -g {input.geno} -m {input.meta} -o {output}
        """
        
rule SNP_significance:
    # Calls the R script `outlier_SNPs.R` which uses pcadapt to calculate the significance of each SNP in explaining genetic variation, captured through a PCA projection.
    input:
        "data/genotypes/geno2.lfmm"
    output:
        "figures/snp_significance.png"
    conda:
        "../envs/pcadapt.yml"
    shell:
        """
        Rscript scripts/outlier_SNPs.R -i {input} -o {output}
        """
        
rule run_tess3r:
    input:
        lfmm = "data/genotypes/geno2.lfmm",
        meta = "config/meta_subset.csv"
    output:
        "figures/tess3_barplot_map_1.png",
        "figures/tess3_cross_validation_plot.png"
    conda:
        "../envs/tess3.yml"
    shell:
        """
        Rscript scripts/run_tess3r.R -l {input.lfmm} -m {input.meta} -o figures/tess3
        """

rule get_allele_stats:
    # Calls the script `get_pairwise_fst.R` to generate pairwise_fst table (calculated according to Nei 1973)
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
        
# rule get_iqtree_input:
#     input:
#         geno = "data/genotypes/geno2.geno",
#         meta = "config/meta_subset.csv"
#     output:
#         "data/genotypes/geno2.phy"
#     shell:
#         """
#         python scripts/angsd_to_phyla.py -g {input.geno} -m {input.meta} -o {output}
#         """
    
# rule run_iqtree:
#     input:
#         "data/genotypes/geno2.phy"
#     output:
#         "iqtree_checkpoint.txt"
#     conda:
#         '../envs/phyla.yml'
#     shell:
#         """
#         iqtree -s {input} -m GTF+ASC -nt 30
#         """
        
        
    