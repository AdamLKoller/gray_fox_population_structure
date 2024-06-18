from subprocess import check_output

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
        return row.Sample_ID.values[0]
    except:
        return None
    
rule clean_and_rename:
    input:
        expand("data/lane{lane}/raw/{library}_log.txt", lane=[1,2], library = range(1,39))
        
    output:
        expand("data/lane{lane}/deduped/{sample}.dedup", lane=[1,2], sample=samples)
       
    run:
        for lane in ['1','2']:
            directory = f'data/lane{lane}/raw/'.replace(' ','')
            for filename in os.listdir(directory):
                if filename.endswith('.tr0'):
                    library, barcode = int(filename.split('_')[0]), filename.split('_')[-1][0:4]
                    sample_name = get_sample_name(library, barcode)
                    if sample_name:                        
                        os.rename(directory+filename, f'data/lane{lane}/deduped/{sample_name}.dedup'.replace(' ',''))
                        print(directory+filename, 'has been renamed to ', f'data/lane{lane}/deduped/{sample_name}.dedup'.replace(' ',''))
                    else:        
                        os.remove(directory+filename)
                        print(directory+filename, 'has been removed')
        
rule count_deduped_reads:
    input:
        expand("data/lane{lane}/deduped/{sample}.dedup", lane=[1,2], sample=samples)
    output:
        "tables/deduped_reads.csv"
    run:
        shell("touch tables/deduped_reads.csv")
        deduped_reads_file = open("tables/deduped_reads.csv",'w')
        deduped_reads_file.write("Sample_ID,deduped_count_lane1,deduped_count_lane2\n")
        directory = 'data/lane1/deduped/'
        for sample_id in [x.split('/')[-1].split('.')[0] for x in os.listdir(directory) if x.endswith('.dedup')]:
            file_l1, file_l2 = f'data/lane1/deduped/{sample_id}.dedup', f'data/lane2/deduped/{sample_id}.dedup' 
            read_count_l1 = wc(file_l1.replace(' ',''))//4
            read_count_l2 = wc(file_l2.replace(' ',''))//4
            deduped_reads_file.write(f'{sample_id},{read_count_l1},{read_count_l2}\n'.replace(' ',''))

rule trim_reads:
    input:
        "data/lane{lane}/deduped/{sample}.dedup"
    output:
        "data/lane{lane}/trimmed/{sample}.trim"
    conda:
        "../envs/preprocessing.yml"
    shell:
        """
        cutadapt -q 15,15 -m 25 -o {output} {input}
        """

def wc(filename):
    return int(check_output(["wc", "-l", filename]).split()[0])
        
rule count_trimmed_reads:
    input:
        expand("data/lane{lane}/trimmed/{sample}.trim", lane=[1,2], sample=samples)
    output:
        "tables/trimmed_reads.csv"
    run:
        shell("touch tables/trimmed_reads.csv")
        trimmed_reads_file = open("tables/trimmed_reads.csv",'w')
        trimmed_reads_file.write("Sample_ID,trimmed_count_lane1,trimmed_count_lane2\n")
        directory = 'data/lane1/trimmed/'
        for sample_id in [x.split('/')[-1].split('.')[0] for x in os.listdir(directory) if x.endswith('.trim')]:
            file_l1, file_l2 = f'data/lane1/trimmed/{sample_id}.trim', f'data/lane2/trimmed/{sample_id}.trim' 
            read_count_l1 = wc(file_l1.replace(' ',''))//4
            read_count_l2 = wc(file_l2.replace(' ',''))//4
            trimmed_reads_file.write(f'{sample_id},{read_count_l1},{read_count_l2}\n'.replace(' ',''))
                    
                    
        
rule align_reads_to_ref:
    input:
        ref_idx1 = 'data/ref_genome/UCinero_ref.1.bt2',
        ref_idx2 = 'data/ref_genome/UCinero_ref.2.bt2',
        ref_idx3 = 'data/ref_genome/UCinero_ref.3.bt2',
        ref_idx4 = 'data/ref_genome/UCinero_ref.4.bt2',
        ref_idx1_rev = 'data/ref_genome/UCinero_ref.rev.1.bt2',
        ref_idx2_rev = 'data/ref_genome/UCinero_ref.rev.2.bt2',
        read_file = "data/lane{lane}/trimmed/{sample}.trim",
    output:
        sam = temp("data/lane{lane}/sams/{sample}.sam"),
        log = "data/lane{lane}/logs/{sample}_align.log"
    conda:
        "../envs/align_reads_to_ref.yml"
    params:
        genome_basename = 'data/ref_genome/UCinero_ref'
    shell:
        '''
        bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x {params.genome_basename} -U {input.read_file} -S {output.sam} > {output.log}
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
        
