rule add_apriori_pops:
    # Calls the script `add_apriori_pops.py` to add another column 'POP'
    # containing apriori sample group to the meta file as defined in 
    # shapefiles in data/shapefiles/XYRed2017POP.*
    input:
        meta = 'config/meta.csv'
    output:
        'config/meta_apriori.csv'
    conda:
        "../envs/analyses.yml"
    shell:
        """
        python scripts/add_apriori_pops.py -m {input} -o {output}
        """ 

rule jitter_sample_locations:
    # Calls the script `jitter_locations.py` to jitter samples with the same location.
    # For some samples, just county or state data was collected. These points were given
    # longitude and latitude values equal to the centroid of their county or state.
    # We want to apply a jitter to these points so they do not overlap in visualizations.
    # The jittered locations are also used in TESS3
    input:
        meta = 'config/meta_apriori.csv'
    output:
        'config/meta_apriori_jittered.csv'
    conda:
        "../envs/jitter.yml"
    shell:
        """
        python scripts/jitter_locations.py -m {input} -o {output}
        """