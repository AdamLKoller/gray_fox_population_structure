rule get_refgenome:
    # Downloads the gray fox reference genome from NCBI
    output:
        'data/ref_genome/GCA_032313775.1_UCinereo1.0_genomic.fna'
    params:
        ftp = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/032/313/775/GCA_032313775.1_UCinereo1.0/GCA_032313775.1_UCinereo1.0_genomic.fna.gz',
        outdir = 'data/ref_genome'
    shell:
        """
        wget -P {params.outdir} {params.ftp}
        gunzip {output}.gz
        """
        
rule build_genome:
    # Builds the reference genome index files neccesary for Bowtie2 alignment
    input:
        reference = 'data/ref_genome/GCA_032313775.1_UCinereo1.0_genomic.fna'
        
    output:
        'data/ref_genome/UCinero_ref.1.bt2',
        'data/ref_genome/UCinero_ref.2.bt2',
        'data/ref_genome/UCinero_ref.3.bt2',
        'data/ref_genome/UCinero_ref.4.bt2',
        'data/ref_genome/UCinero_ref.rev.1.bt2',
        'data/ref_genome/UCinero_ref.rev.2.bt2'
    params:
        basename = 'data/ref_genome/UCinero_ref'
    shell: 
        """
        bowtie2-build {input.reference} {params.basename} -p 35 
        """