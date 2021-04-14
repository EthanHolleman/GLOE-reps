
rule bed_to_coverage_bedgraph:
# be careful does not consider strand
    conda:
        '../envs/bedtools.yml'
    input:
        bed='{filepath}.bed',
        chrom_sizes='rawdata/hg19/hg19.chrom.sizes'
    output:
        bedgraph='{filepath}.coverage.bedgraph',
        sort=temp('{filepath}.sorted.bed')
    shell:'''
    bedtools sort -i {input} > {output.sort}
    bedtools genomecov -bg -i {input.bed} -g {input.chrom_sizes} \
    > {output.bedgraph}
    '''


rule bed_to_bedgraph:
    input:
        '{file}.bed'
    output:
        '{file}.bedgraph'
    shell:'''
    cut -f 1,2,3,5 {input} > {output} && [[ -s {output} ]]
    '''


rule bedgraph_to_bigwig:
    conda:
        '../envs/ucsc.yml'
    input:
        '{filepath}.bedgraph',
        chrom_sizes='rawdata/hg19/hg19.chrom.sizes'
    output:
        '{filepath}.bw'
    shell:'''
    bedGraphToBigWig {input} {input.chrom_sizes} {output}
    '''
