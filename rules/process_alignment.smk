# process bowtie2 output files into desired formats and trim low quality
# alignments

rule sam_to_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/{sample}/mapped/{sample}.sam'
    output:
        'output/{sample}/mapped/{sample}.bam'
    shell:'''
    samtools view -bhSu -o {output} {input}
    '''


rule sort_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/{sample}/mapped/{sample}.bam'
    output:
        temp('output/{sample}/mapped/{sample}.sorted.bam')
    params:
        sort='output/{sample}_temp'
    threads: 16
    shell:'''
    mkdir --parents {params.sort}
    samtools sort -O bam -T {params.sort} --threads {threads} -o {output} {input}
    rm -r {params.sort}
    '''


rule trim_bam_low_qual_alignments:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/{sample}/mapped/{sample}.sorted.bam'
    output:
        temp('output/{sample}/mapped/{sample}.trim.all.bam')
    shell:'''
    samtools view -q 30 -bhu -o {output} {input}
    '''


rule trim_bam_low_qual_alignments_footloop_only:
    input:
        footloop_amplicons='resources/footprinted_sites.bed'
        sorted_bam='output/{sample}/mapped/{sample}.sorted.bam'
    output:
        temp('output/{sample}/mapped/{sample}.trim.footloop.bam')
    shell:'''
    samtools view -q 30 -L {input.footloop_amplicons} -bhu -o {output} \
    {input.sorted_bam}
    '''


rule sort_trimmed_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/{sample}/mapped/{sample}.trim.{region}.bam'
    output:
        'output/{sample}/mapped/{sample}.sorted.trim.{region}.bam'
    threads: 16
    params:
        sort='output/alignment/{sample}.{region}.trim.temp'
    shell:'''
    mkdir --parents {params.sort}
    samtools sort -O bam -T {params.sort} --threads {threads} -o {output} {input}
    rm -r {params.sort}
    '''


# rule bam_to_bed:
#     conda: 
#         '../envs/bedtools.yml'
#     input:
#         'output/{sample}/mapped/{sample}.sorted.trim.bam'
#     output:
#         temp('output/{sample}/mapped/{sample}.sorted.trim.bed')
#     shell:'''
#     bedtools bamtobed -i {input} > {output}
#     '''

# rule remove_low_read_count_breaks:
#     conda:
#         '../envs/bedtools.yml'


rule read_depth:
    input:
        'output/{sample}/mapped/{sample}.sorted.trim.{region}.bam'
    output:
        'output/{sample}/mapped/{sample}.sorted.trim.depth.{region}.bed'
    shell:'''
    samtools depth {input} > {output}
    '''


rule remove_shallow_breaks:
    conda:
        '../envs/R.yml'
    input:
        'output/{sample}/mapped/{sample}.sorted.trim.depth.{region}.bed'
    output:
        'output/{sample}/mapped/{sample}.sorted.trim.depth.deep.{region}.bed'
    shell:'''
    Rscript scripts/call_breaks.R {input} {output}
    '''








