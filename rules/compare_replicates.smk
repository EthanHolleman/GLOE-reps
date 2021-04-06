

rule compare_replicate_closeness:
    conda:
        '../envs/bedtools.yml'
    input:
        rep1='output/{sample_1}/depth/{mode}/{sample_1}.{mode}.sorted.trim.{region}.depth.deep.{strand}.bed',
        rep2='output/{sample_2}/depth/{mode}/{sample_2}.{mode}.sorted.trim.{region}.depth.deep.{strand}.bed',
        genome='rawdata/hg19/hg19.chrom.sizes'
    output:
        'output/compare_replicates/{sample_1}-{sample_2}/{mode}/closest.{region}.{strand}.bed'
    params:
        out_dir='output/compare_replicates/{sample_1}-{sample_2}'
    shell:'''
    mkdir -p {params}
    bedtools closest -d -t first -s -g {input.genome} -a {input.rep1} -b {input.rep2} > {output}
    '''


rule compare_replicate_coverage:
    conda:
        '../envs/bedtools.yml'
    input:
        rep1='output/{sample_1}/depth/direct/{sample_1}.direct.sorted.trim.{region}.depth.deep.{strand}.bed',
        rep2='output/{sample_2}/depth/direct/{sample_2}.direct.sorted.trim.{region}.depth.deep.{strand}.bed'
    output:
        'output/compare_replicates/{sample_1}-{sample_2}/coverage.{region}.bed'
    params:
        out_dir='output/compare_replicates/{sample_1}-{sample_2}'
    shell:'''
    mkdir -p {params.out_dir}
    bedtools coverage -a -s -sorted {input.rep1} -b {input.rep2} > {output}
    '''

# think we will just combine boxplots later on for now 

rule plot_replicate_closeness:
    conda:
        '../envs/R.yml'
    input:
        closness='output/compare_replicates/{sample_1}-{sample_2}/direct/closest.{region}.{strand}.bed',
        samples='samples/GLOE_samples.csv'
    output:
        'output/compare_replicates/{sample_1}-{sample_2}/plots/closest/{strand}/boxplot_{region}.png'
    params:
        out_dir = 'output/compare_replicates/{sample_1}-{sample_2}/plots/closest/{strand}'
    shell:'''
    mkdir -p {params.out_dir}
    Rscript scripts/plot_replicate_dists.R {input.closness} {input.samples} {output}
    '''
        
