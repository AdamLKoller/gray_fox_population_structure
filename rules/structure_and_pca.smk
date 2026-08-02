import shutil


# STRUCTURE
################################################################################################################


rule prepare_data_for_structure:
    # Runs python script `angsd_to_structure` to transform ANGSD output into STRUCTURE input
    input:
        geno = "data/genotypes/geno5.geno",
    output:
        "data/structure/structure_input.tab"
    shell:
        """
        python scripts/angsd_to_structure.py -g {input.geno} -o {output}
        """


rule create_structure_parameter_files:
    # Creates the parameter files used by STRUCTURE
    input:
        freq = "data/genotypes/geno5.geno",
        bam_list = "data/merged_bams/bams_4_rcr.txt"
    output:
        mainparams = "data/structure/mainparams.txt",
        extraparams = "data/structure/extraparams.txt"
        
    run:
        # Get the number of loci using shell command
        numloci = int(shell("wc -l < {input.freq}", read=True).strip())
        numinds = int(shell("wc -l < {input.bam_list}", read=True).strip())

        # Write mainparams file
        with open(output.mainparams, 'w') as mainparams:
            mainparams.write("#define MAXPOPS 10\n")
            mainparams.write("#define BURNIN 25000\n")
            mainparams.write("#define NUMREPS 100000\n")
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
        echo "structure:     $(structure --version)" >> versions.txt
        structure -m {input.mainparams} -e {input.extraparams} -K {wildcards.k} -o data/structure/structure_output_k{wildcards.k}_rep{wildcards.rep} -D {wildcards.rep} &> data/structure/structure_log_k{wildcards.k}_rep{wildcards.rep}
        """
        
rule harvest_structure:
    # Uses structureHarvester to generate summary with the likelihood of the data given k
    input:
        expand("data/structure/structure_output_k{k}_rep{rep}_f", k=range(1,11), rep=range(1,11))
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
        expand("data/structure/structure_output_k{k}_rep{rep}_f", k=range(1,11), rep=range(1,11))
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
        
        
        
# PCA
################################################################################################################

rule transform_angsd_for_PCA:
    # Calls the script `angsd_to_pca.py` which is responsible for transforming ANGSD output data into a format usable for PCA
    # Encodes the homozygous for major allele as 0, heterozygous as 1, and homozygous minor as 2 (i.e. number of minor alleles)
    # Missing values are imputed with the average for that loci. This has the effect of pulling points towards the origin of
    # the ordinate space, which is undesirable. However, since missing data is quite low, this is less of an issue.
    # Additionally, if clusters persist despite being pulled to the center, they must be robust, meaning averaging only decreases
    # our power to detect clustering, but should not lead to false positives, only false negatives.
    input:
        geno = "data/genotypes/geno5.geno",
    output:
        "data/genotypes/genotype_matrix.csv"
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/angsd_to_pca.py -g {input.geno} -o {output} 
        """
        
        
# Generating figures 2 and 3
##################################################################################################################

rule get_combined_plots:
    # Calls the script `get_combined_plots.py` 
    # which produces a single figure with the STRUCTURE bar plot, pie map, and PCA plot for each k.
    # Among the 10 plots for each K produced, are figures 2 and 3 of our paper.
    # The coloring can be quite janky, and will, without a doubt, have an issue depending on the K as a result of my poor programming.
    # The list of colors within `get_combined_plots.py` should be adjusted for the particular K you are interested in--my apologies.
    input:
        q_matrix = "data/structure/q_matrix_{k}",
        meta = "config/meta_apriori_jittered.csv",
        geno_matrix = "data/genotypes/genotype_matrix.csv",
        bam_list = "data/merged_bams/bams_4_rcr.txt"
    output:
        "figures/structure_PCA_pie_{k}.png"
    conda:
        "../envs/plotting.yml"
    shell:
        """
        python scripts/get_combined_plots.py -q {input.q_matrix} -g {input.geno_matrix} -m {input.meta} -o {output} -k {wildcards.k} -b {input.bam_list}
        """





