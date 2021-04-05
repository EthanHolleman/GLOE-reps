# process bowtie2 output files into desired formats and trim low quality
# alignments

rule sam_to_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/{sample}/bowtie2/{sample}.sam'
    output:
        'output/{sample}/process_alignment/{sample}.bam'
    shell:'''
    samtools view -bhSu -o {output} {input}
    '''


rule sort_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/{sample}/process_alignment/{sample}.bam'
    output:
        temp('output/{sample}/process_alignment/{sample}.sorted.bam')
    params:
        sort='output/{sample}_temp'
    threads: 16
    shell:'''
    mkdir --parents {params.sort}
    samtools sort -O bam -T {params.sort} --threads {threads} -o {output} {input}
    rm -r {params.sort}
    '''


rule trim_bam_low_qual_alignments_all:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/{sample}/process_alignment/{sample}.sorted.bam'
    output:
        temp('output/{sample}/process_alignment/{sample}.trim.all.bam')
    shell:'''
    samtools view -q 30 -bhu -o {output} {input}
    '''


rule trim_bam_low_qual_alignments_footloop_only:
    input:
        footloop_amplicons='resources/footprinted_sites.bed',
        sorted_bam='output/{sample}/process_alignment/{sample}.sorted.bam'
    output:
        temp('output/{sample}/process_alignment/{sample}.trim.footloop.bam')
    shell:'''
    samtools view -q 30 -L {input.footloop_amplicons} -bhu -o {output} \
    {input.sorted_bam}
    '''


rule sort_trimmed_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/{sample}/process_alignment/{sample}.trim.{region}.bam'
    output:
        temp('output/{sample}/process_alignment/{sample}.sorted.trim.{region}.bam')
    threads: 16
    params:
        sort='output/alignment/{sample}.{region}.trim.temp'
    shell:'''
    mkdir --parents {params.sort}
    samtools sort -O bam -T {params.sort} --threads {threads} -o {output} {input}
    rm -r {params.sort}
    '''


rule bam_to_bed:
    conda: 
        '../envs/bedtools.yml'
    input:
        'output/{sample}/process_alignment/{sample}.sorted.trim.{region}.bam'
    output:
        temp('output/{sample}/process_alignment/{sample}.sorted.trim.{region}.bed')
    shell:'''
    bedtools bamtobed -i {input} > {output}
    '''

# rule remove_low_read_count_breaks:
#     conda:
#         '../envs/bedtools.yml'









