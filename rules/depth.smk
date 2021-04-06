
rule bed_to_bam:
    conda:
        '../envs/bedtools.yml'
    input:
        bed='output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.{strand}.bed',
        genome='rawdata/hg19/hg19.chrom.sizes'
    params:
        out_dir = 'output/{sample}/depth/{mode}/'
    output:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.{strand}.bam'
    shell:'''
    mkdir -p {params.out_dir}
    bedtools bedtobam -i {input.bed} -g {input.genome} > \
    {output}
    '''

rule sort_bam_for_read_depth:
    conda:
        '../envs/samtools.yml'
    input:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.{strand}.bam'
    output:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.sorted.{strand}.bam'
    params:
        sort='output/{sample}/depth/{mode}/temp'
    threads: 16
    shell:'''
    mkdir -p {params.sort}
    samtools sort -O bam -T {params.sort} --threads {threads} {input} > {output}
    rm -r {params.sort}
    '''


rule read_depth:
    input:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.sorted.{strand}.bam'
    output:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.depth.{strand}.bed'
    shell:'''
    samtools depth {input} > {output}
    '''


rule remove_shallow_breaks_fwd:
    conda:
        '../envs/R.yml'
    input:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.depth.fwd.bed'
    output:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.depth.deep.fwd.bed'
    params:
        strand='+'
    shell:'''
    Rscript scripts/call_breaks.R {input} {output} {params.strand}
    '''

rule remove_shallow_breaks_rev:
    conda:
        '../envs/R.yml'
    input:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.depth.rev.bed'
    output:
        'output/{sample}/depth/{mode}/{sample}.{mode}.sorted.trim.{region}.depth.deep.rev.bed'
    params:
        strand='-'
    shell:'''
    Rscript scripts/call_breaks.R {input} {output} {params.strand}
    '''