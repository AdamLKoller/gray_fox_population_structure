rule decompress:
    input:
        "data/lane{lane}/raw/{library}_S{library}_L00{lane}_R1_001.fastq.gz"
    output:
        temp("data/lane{lane}/raw/{library}_S{library}_L00{lane}_R1_001.fastq")
    conda:
        "../envs/preprocessing.yml"
    shell:
        "gzip -d -k {input}"
        
rule demultiplex:
    # export PATH="/home/LC/kollad01/gray_fox_east_structure/programs/2bRAD_denovo/:$PATH"
    input:
        "data/lane{lane}/raw/{library}_S{library}_L00{lane}_R1_001.fastq"
    output:
        temp("data/lane{lane}/raw/{library}_log.txt")
    conda:
        "../envs/preprocessing.yml"
        
    shell:
        '''
        trim2bRAD_2barcodes_dedup.pl input={input} site=".{{12}}CGA.{{6}}TGC.{{12}}|.{{12}}GCA.{{6}}TCG.{{12}}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{{4}} &> {output}
        '''

        
def get_sample_name(library,barcode):
    try:
        row = meta.loc[(meta.library==library) & (meta.barcode==barcode)]
        return row.id.values[0]
    except:
        return None
    
rule clean_and_rename:
    input:
        expand("data/lane{lane}/raw/{library}_log.txt", lane=[1,2], library = range(1,39))
        
    output:
        expand("data/lane{lane}/deduped/{sample}.dedup", lane=[1,2], sample=samples)
       
    run:
        for lane in ['1','2']:
            directory = f'data/lane{lane}/raw/'
            for filename in os.listdir(directory):
                if filename.endswith('.tr0'):
                    library, barcode = int(filename.split('_')[0]), filename.split('_')[-1][0:4]
                    sample_name = get_sample_name(library, barcode)
                    if sample_name:
                        os.rename(directory+filename, f'data/lane{lane}/deduped/{sample_name}.dedup')
                    else:
                        os.remove(directory+filename)
        
        
rule align_reads_to_ref:
    input:
        ref_idx1 = 'data/ref_genome/UCinero_ref.1.bt2',
        ref_idx2 = 'data/ref_genome/UCinero_ref.2.bt2',
        ref_idx3 = 'data/ref_genome/UCinero_ref.3.bt2',
        ref_idx4 = 'data/ref_genome/UCinero_ref.4.bt2',
        ref_idx1_rev = 'data/ref_genome/UCinero_ref.rev.1.bt2',
        ref_idx2_rev = 'data/ref_genome/UCinero_ref.rev.2.bt2',
        read_file = "data/lane{lane}/deduped/{sample}.dedup",
    output:
        temp("data/lane{lane}/sams/{sample}.sam")
    conda:
        "../envs/align_reads_to_ref.yml"
    params:
        genome_basename = 'data/ref_genome/UCinero_ref'
    shell:
        '''
        bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x {params.genome_basename} -U {input.read_file} -S {output}
        '''
rule merge_sams:
    input:
        lane1_file = "data/lane1/sams/{sample}.sam",
        lane2_file = "data/lane2/sams/{sample}.sam",
    output:
        "data/merged_bams/{sample}.sam"
    conda:
        "../envs/preprocessing.yml"
    shell:
        "samtools merge -o {output} {input.lane1_file} {input.lane2_file}"
        
        
rule convert_sams_to_bams:
    input:
        "data/merged_bams/{sample}.sam"
    output:
        "data/merged_bams/{sample}.bam"
    conda:
        "../envs/preprocessing.yml"
    shell:
        "samtools sort -O bam -o {output} {input}"

rule build_bam_index_files:
    input:
        "data/merged_bams/{sample}.bam"
    output:
        "data/merged_bams/{sample}.bam.bai"
    conda:
        "../envs/preprocessing.yml"
    shell:
        "samtools index {input}"
          
rule make_list_of_bams:
    input:
        bam_index_files = expand("data/merged_bams/{sample}.bam.bai",sample=samples),
        bam_files =expand("data/merged_bams/{sample}.bam", sample=samples)
        
    output:
        "data/merged_bams/bams.txt"
    conda:
        "../envs/preprocessing.yml"
    shell:
        "ls data/merged_bams/*.bam > data/merged_bams/bams.txt"
        
