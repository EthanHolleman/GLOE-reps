

rule compare_replicate_closeness:
    conda:
        '../envs/bedtools.yml'
    input:
        rep1='output/{sample_1}/mapped/{sample_1}.sorted.trim.depth.deep.{region}.bed',
        rep2='output/{sample_2}/mapped/{sample_2}.sorted.trim.depth.deep.{region}.bed'
    output:
        'output/compare_replicates/{sample_1}-{sample_2}/closest.{region}.bed'
    shell:'''
    bedtools closest -d -t first -s -a {input.rep1} -b {input.rep2} > {output}
    '''

rule compare_replicate_coverage:
    conda:
        '../envs/bedtools.yml'
    input:
        rep1='output/{sample_1}/mapped/{sample_1}.sorted.trim.depth.deep.{region}.bed',
        rep2='output/{sample_2}/mapped/{sample_2}.sorted.trim.depth.deep.{region}.bed'
    output:
         'output/compare_replicates/{sample_1}-{sample_2}/coverage.{region}.bed'
    shell:'''
    bedtools coverage -a -s -sorted {input.rep1} -b {input.rep2} > {output}
    '''


