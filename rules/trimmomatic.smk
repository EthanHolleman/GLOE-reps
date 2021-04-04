
rule trimmomatic:
    conda: 
        '../envs/trimmomatic.yml'
    input:
        reads='rawdata/GLOE-seq/{sample}.fastq.gz',
        primers='rawdata/primers/TruSeq3-SE.fa'
    output:
        'output/{sample}/trimmomatic/{sample}.trimmed.fastq.gz'
    threads: 16
    shell:'''
    trimmomatic SE -threads {threads} -phred33 {input.reads} {output} \
    ILLUMINACLIP:{input.primers}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
    '''