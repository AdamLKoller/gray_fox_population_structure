from subprocess import check_output

rule merge_libraries:
    # Merges libraries across lanes prior to demultiplexing
    input:
        lane1_file = "data/lane1/{library}_S{library}_L001_R1_001.fastq.gz",
        lane2_file = "data/lane2/{library}_S{library}_L002_R1_001.fastq.gz"
    output:
        temp("data/merged_reads/library_{library}.fastq.gz")
    threads:
        2
    shell:
        """
        cat {input.lane1_file} {input.lane2_file} > {output} 
        """
        

rule decompress_fastq:
    # Decompresses raw .fastq files
    input:
        "data/merged_reads/library_{library}.fastq.gz"
    output:
        temp("data/merged_reads/library_{library}.fastq")
    threads:
        2
    shell:
        """
        gzip -d -k {input}
        """        
        
rule demultiplex:
    # Demultiplexes reads according to their 4bp barcode
    # This step creates extra files not associated with a sample (will be removed later)
    input:
        "data/merged_reads/library_{library}.fastq"
    output:
        "data/merged_reads/library_{library}_demulti_log.txt"
    conda:
        "../envs/preprocessing.yml"
    threads:
        2
    shell:
        """
        trim2bRAD_2barcodes_dedup.pl input={input} site=".{{12}}CGA.{{6}}TGC.{{12}}|.{{12}}GCA.{{6}}TCG.{{12}}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{{4}} &> {output}
        """

        
def get_sample_name(library,barcode):
    """
    Helper function for the `clean_and_rename` rule
    """
    try:
        row = meta.loc[(meta.library==library) & (meta.barcode==barcode)]
        return row.Sample_ID.values[0]
    except:
        return None
    
rule clean_and_rename:
    # Removes extra files from demultiplexing that were not associated with a sample
    # Those that are associated with a sample are renamed and moved to the deduped/ directory
    input:
        expand("data/merged_reads/library_{library}_demulti_log.txt", library = range(1,39))
        
    output:
        expand("data/merged_reads/{sample}.dedup", sample=samples)
       
    run:
        directory = 'data/merged_reads/'
        for filename in os.listdir(directory):
            if filename.endswith('.tr0'):
                library, barcode = int(filename.split('_')[1]), filename.split('_')[-1].rstrip('.tr0')
                sample_name = get_sample_name(library, barcode)
                if sample_name:                        
                    os.rename(directory+filename, directory+sample_name+'.dedup')
                    print(directory+filename, 'has been renamed to', directory+sample_name+'.dedup')
                    
                else:        
                    os.remove(directory+filename)
                    print(directory+filename, 'has been removed')
                    
def wc(filename):
    """
    helper function for `count_trimmed_reads`
    gets the word count for a file
    """
    return int(check_output(["wc", "-l", filename]).split()[0])
        
        
rule count_deduped_reads:
    # Counts the number of reads for each deduped sample and stores in tables/deduped_reads.csv
    input:
        expand("data/merged_reads/{sample}.dedup", sample=samples)
    output:
        "tables/deduped_reads.csv"
    run:
        shell("touch tables/deduped_reads.csv")
        deduped_reads_file = open("tables/deduped_reads.csv",'w')
        deduped_reads_file.write("Sample_ID,deduped_count\n")
        for file in [x for x in os.listdir('data/merged_reads') if x.endswith('.dedup')]:
            read_count = wc('data/merged_reads/'+file) // 4
            sample_id = file.split('/')[-1].split('.')[0]
            deduped_reads_file.write(sample_id + ',' + str(read_count) + '\n')
            
rule trim_reads:
    # Trims low-quality bp at the end of reads using cutadapt. 
    # Reads with < 25bp of quality bases are removed entirely
    input:
        "data/merged_reads/{sample}.dedup"
    output:
        temp("data/merged_reads/{sample}.trim")
    conda:
        "../envs/preprocessing.yml"
    shell:
        """
        echo "cutadapt:     $(cutadapt --version)" >> versions.txt
        cutadapt -q 15,15 -m 25 -o {output} {input}
        """


rule count_trimmed_reads:
    # Counts the number of reads for each trimmed sample and stores in tables/trimmed_reads.csv
    input:
        expand("data/merged_reads/{sample}.trim", sample=samples)
    output:
        "tables/trimmed_reads.csv"
    run:
        shell("touch tables/trimmed_reads.csv")
        trimmed_reads_file = open("tables/trimmed_reads.csv",'w')
        trimmed_reads_file.write("Sample_ID,trimmed_count\n")
        for file in [x for x in os.listdir('data/merged_reads') if x.endswith('.trim')]:
            read_count = wc('data/merged_reads/' + file) // 4
            sample_id = file.split('/')[-1].split('.')[0]
            trimmed_reads_file.write(sample_id + ',' + str(read_count) + '\n')
                    
rule align_reads_to_ref:
    # Aligns reads to the gray fox reference genome using bowtie2
    input:
        ref_idx1 = 'data/ref_genome/UCinero_ref.1.bt2',
        ref_idx2 = 'data/ref_genome/UCinero_ref.2.bt2',
        ref_idx3 = 'data/ref_genome/UCinero_ref.3.bt2',
        ref_idx4 = 'data/ref_genome/UCinero_ref.4.bt2',
        ref_idx1_rev = 'data/ref_genome/UCinero_ref.rev.1.bt2',
        ref_idx2_rev = 'data/ref_genome/UCinero_ref.rev.2.bt2',
        read_file = "data/merged_reads/{sample}.trim",
    output:
        sam = temp("data/merged_bams/{sample}.sam"),
        log = temp("data/merged_bams/{sample}_align.log")
    conda:
        "../envs/align_reads_to_ref.yml"
    params:
        genome_basename = 'data/ref_genome/UCinero_ref'
    shell:
        """
        echo "bowtie2:     $(bowtie2 --version)" >> versions.txt
        bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x {params.genome_basename} -U {input.read_file} -S {output.sam} &> {output.log}
        """
        
rule get_alignment_rate_data:
    # Gathers alignment rate statistics for each sample in stores in tables/aligned_reads.csv
    input:
        expand("data/merged_bams/{sample}_align.log", sample=samples)
    output:
        "tables/aligned_reads.csv",
        
    run:
        
        shell("touch tables/aligned_reads.csv")
        table_file = open("tables/aligned_reads.csv", 'w')
        table_file.write("Sample_ID,aligned_reads,aligned_0_times,aligned_exactly_once,aligned_more_than_once,overall_alignment_rate\n")
        
        for file in [x for x in os.listdir('data/merged_reads/') if x.endswith('align.log')]:
            sample_id = file.split('/')[-1].split('_')[0]
            log_file = open(log_file, 'r')
            to_write += log_file.readline().split(' ')[0] + ','
            log_file.readline()
            to_write += log_file.readline().lstrip().split(' ')[0] + ','
            to_write += log_file.readline().lstrip().split(' ')[0] + ','
            to_write += log_file.readline().lstrip().split(' ')[0] + ','
            to_write += str(float(log_file.readline().lstrip().split(' ')[0].rstrip('%'))/100) + ','
            to_write = to_write.rstrip(',')
            to_write += '\n'
            table_file.write(to_write)
        
                
rule convert_sams_to_bams:
    # Converts sam files to bam files so that they can be used as input for ANGSD
    input:
        "data/merged_bams/{sample}.sam"
    output:
        "data/merged_bams/{sample}.bam"
    conda:
        "../envs/preprocessing.yml"
    shell:
        """
        samtools sort -O bam -o {output} {input}
        """

rule build_bam_index_files:
    # Creates bam index files
    input:
        "data/merged_bams/{sample}.bam"
    output:
        "data/merged_bams/{sample}.bam.bai"
    conda:
        "../envs/preprocessing.yml"
    shell:
        """
        samtools index {input}
        """
          
