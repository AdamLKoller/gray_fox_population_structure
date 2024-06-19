import pandas as pd

meta = pd.read_csv('config/meta.csv')

barcodes = meta['barcode'].unique().tolist()
samples = meta['Sample_ID'].tolist()

        
rule all:
    input:
        "figures/structure_barplot.png",
        "figures/delta_k.png",
        expand("figures/pca_plot_{k}.png", k = range(1,11)),
        expand("figures/pie_map_{k}.png", k = range(1,11)),
        expand("figures/structure_PCA_pie_{k}.png", k = range(1,11)),
        "figures/heterozygosity_map.png",
        "data/genotypes/ngs_relate_out",
        "tables/deduped_reads.csv",
        "tables/trimmed_reads.csv",
        expand("tables/pairwise_fst_{k}", k = range(2,11))
        


include:'rules/get_data.smk'
include:'rules/preprocessing.smk'
include:'rules/call_variants.smk'
include:'rules/analyses.smk'
include: 'rules/checking_relatedness.smk'
    

            

        

        
        
        

        
    
    
        