import pandas as pd

meta = pd.read_csv('config/meta.csv')

barcodes = meta['barcode'].unique().tolist()
samples = meta['Sample_ID'].tolist()

        
rule all:
    input:
        "figures/structure_barplot.png",
        "figures/delta_k.png",
        expand("figures/pie_map_{k}.png", k = range(1,11)),
        expand("figures/structure_PCA_pie_{k}.png", k = range(1,11)),
        expand("figures/structure_PCA_pie_apriori_{k}.png", k = range(1,11)),
        "figures/heterozygosity_map.png",
        #"tables/deduped_reads.csv",
        #"tables/trimmed_reads.csv",
        expand("tables/pairwise_fst_{k}", k = range(2,11)),
        "figures/snp_significance.png",
        "figures/tess3_barplot_map_1.png",
        expand("tables/allele_stats_{k}.tab", k = range(2,5)),
        "figures/tess3_cross_validation_plot.png",
        #"iqtree_checkpoint.txt"
        #"tables/aligned_reads.csv"
        
        


include:'rules/get_data.smk'
include:'rules/preprocessing.smk'
include:'rules/call_variants.smk'
include:'rules/analyses.smk'
    

            

        

        
        
        

        
    
    
        