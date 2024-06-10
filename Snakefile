import pandas as pd

meta = pd.read_csv('config/meta.csv')

barcodes = meta['barcode'].unique().tolist()
samples = meta['Sample_ID'].tolist()

        
rule all:
    input:
        "run_structure_checkpoint.txt"


include:'rules/get_data.smk'
include:'rules/preprocessing.smk'
include:'rules/call_variants.smk'
include:'rules/analyses.smk'
    

            

        

        
        
        

        
    
    
        