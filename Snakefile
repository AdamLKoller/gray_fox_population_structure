import pandas as pd

meta = pd.read_csv('config/meta.csv')

barcodes = meta['barcode'].unique().tolist()
samples = meta['Sample_ID'].tolist()

        
rule all:
    input:
        "figures/delta_k.png",
        expand("figures/structure_PCA_pie_{k}.png", k = range(1,11)),
        #expand("figures/structure_PCA_pie_apriori_{k}.png", k = range(1,11)),
        "figures/heterozygosity_map.png",
        "tables/deduped_reads.csv",
        "tables/trimmed_reads.csv",
        "tables/aligned_reads.csv",
        #expand("tables/pairwise_fst_{k}", k = range(2,11)),
        "figures/snp_significance.png",
        expand("tables/allele_stats_{k}.tab", k = range(2,5)),
        "figures/tess3_cross_validation_plot.png",
        "tables/pairwise_fst_apriori",
        "tables/aligned_reads.csv",
        "tables/allele_stats_apriori.tab",
        "data/genotypes/geno5.counts",
        "tables/raw_reads.csv"
        
        


include:'rules/get_data.smk'
include:'rules/preprocessing.smk'
include:'rules/snp_filtering.smk'
include:'rules/analyses.smk'

    

            

        

        
        
        

        
    
    
        