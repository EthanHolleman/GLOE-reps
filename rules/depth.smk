
rule bed_to_bam:
    conda:
        '../envs/bedtools.yml'
    input:
        bed='output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.bed',
        genome='rawdata/bowtie2/hg19_index'
    params:
        index_name='hg19_index'
    output:
        temp('output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.bam')
    shell:'''
    bedtools bedtobam -i {input.bed} -g {input.genome}/{params.index_name} > \
    {output}
    '''

rule read_depth:
    input:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.bam'
    output:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.depth.bed'
    shell:'''
    samtools depth {input} > {output}
    '''


rule remove_shallow_breaks:
    conda:
        '../envs/R.yml'
    input:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.depth.bed'
    output:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.depth.deep.bed'
    shell:'''
    Rscript scripts/call_breaks.R {input} {output}
    '''