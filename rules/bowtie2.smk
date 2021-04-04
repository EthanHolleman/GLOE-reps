
rule download_hg_bt_index:
    output:
        index_dir=directory('rawdata/bowtie2/hg19_index'),
        downloaded_zip='rawdata/bowtie2/hg19_index/hg19.zip',
    shell:'''
    mkdir --parents {output.index_dir}
    curl https://genome-idx.s3.amazonaws.com/bt/hg19.zip -o {output.downloaded_zip}
    unzip {output.downloaded_zip} -d {output.index_dir}
    '''


rule map_reads:
    conda: 
        'envs/bowtie2.yml'
    input:
        sample_reads='output/{sample}/trimmomatic/{sample}.trimmed.fastq.gz',
        bt_index='rawdata/bowtie2/hg19_index',
    output:
        temp('output/{sample}/mapped/{sample}.sam')
    threads: 12
    shell:'''
    bowtie2 -q -x {input.bt_index}/hg19 -U {input.sample_reads} \
    -p {threads} -S {output}
    '''
    

