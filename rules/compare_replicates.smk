

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
        out_dir='output/compare_replicates/{sample_1}-{sample_2}/{strand}'
    shell:'''
    mkdir -p {params.out_dir}
    bedtools coverage -a -s -sorted {input.rep1} -b {input.rep2} > {output}
    '''


rule intersect_replicates:
    # intersect replicates to calcuate proportion of overlapping breaks
    # between replicates, this could be what? Barplot of somekind
    # or venn diagram might actually be kind of nice.
    conda:
        '../envs/bedtools.yml'
    input:
        rep1='output/{sample_1}/depth/direct/{sample_1}.direct.sorted.trim.{region}.depth.deep.{strand}.bed',
        rep2='output/{sample_2}/depth/direct/{sample_2}.direct.sorted.trim.{region}.depth.deep.{strand}.bed'
    output:
        'output/compare_replicates/{sample_1}-{sample_2}/intersect.{region}.bed'
    shell:'''
    bedtools intersect -wa -wb -sorted -a {input.rep1} -b {input.rep2} > {output}
    '''

rule unique_breaks_rep1:
    conda:
        '../envs/bedtools.yml'
    input:
        rep1='output/{sample_1}/depth/direct/{sample_1}.direct.sorted.trim.{region}.depth.deep.{strand}.bed',
        rep2='output/{sample_2}/depth/direct/{sample_2}.direct.sorted.trim.{region}.depth.deep.{strand}.bed'
    output:
        'output/compare_replicates/{sample_1}-{sample_2}/intersect.unique.{sample_1}.{region}.bed'
    shell:'''
    bedtools intersect -v -a {input.rep1} -b {input.rep2} > {output}
    '''


rule unique_breaks_rep2:
    conda:
        '../envs/bedtools.yml'
    input:
        rep1='output/{sample_1}/depth/direct/{sample_1}.direct.sorted.trim.{region}.depth.deep.{strand}.bed',
        rep2='output/{sample_2}/depth/direct/{sample_2}.direct.sorted.trim.{region}.depth.deep.{strand}.bed'
    output:
        'output/compare_replicates/{sample_1}-{sample_2}/intersect.unique.{sample_2}.{region}.bed'
    shell:'''
    bedtools intersect -v -a {input.rep2} -b {input.rep1} > {output}
    '''


rule breaks_venn_diagram:
    conda:
        '../envs/bedtools.yml'
    input:
        unique_rep1='output/compare_replicates/{sample_1}-{sample_2}/intersect.unique.{sample_1}.{region}.bed',
        unique_rep2='output/compare_replicates/{sample_1}-{sample_2}/intersect.unique.{sample_2}.{region}.bed',
        intersect='output/compare_replicates/{sample_1}-{sample_2}/intersect.{region}.bed'
    output:
        'output/compare_replicates/{sample_1}-{sample_2}/plots/venn/{strand}/venn_{region}.png'
    params:
        out_dir='output/compare_replicates/{sample_1}-{sample_2}/plots/venn/{strand}'
    shell:'''
    Rscript scripts/nick_venn_diagram.R {input.unique_rep1} {input.unique_rep2} {intersect}
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
        
