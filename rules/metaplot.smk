

rule seperate_fwd_strand_hg19_genes:
    input:
        'rawdata/hg19/hg19_apprisplus.bed'
    output:
        'rawdata/hg19/hg19_apprisplus.fwd.bed'
    shell:'''
    grep "+" {input} > {output}
    '''


rule seperate_rev_strand_hg19_genes:
    input:
        'rawdata/hg19/hg19_apprisplus.bed'
    output:
        'rawdata/hg19/hg19_apprisplus.rev.bed'
    shell:'''
    grep "-" {input} > {output}
    '''
    


rule intersect_hg19_genes:
    conda:
        '../envs/bedtools.yml'
    input:
        genes='rawdata/hg19/hg19_apprisplus.{strand}.bed',
        breaks='output/call_peaks/swapped/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.bed'
    params:
        out_dir='output/metaplot/intersection/{peak_call_method}'
    output:
        'output/metaplot/intersection/{peak_call_method}/{control}.vs._{treatment_a}_{treatment_b}.{mode}.{region}.{strand}.intersect.bed'
    shell:'''
    mkdir -p {params.out_dir}
    bedtools intersect -wa -wb -a {input.genes} -b {input.breaks} > {output}
    '''


rule convert_intersection_to_metabed:
    conda:
        '../envs/python.yml'
    input:
        intersection='output/metaplot/intersection/{peak_call_method}/{control}.vs._{treatment_a}_{treatment_b}.{mode}.{region}.{strand}.intersect.bed',
        genes='rawdata/hg19/hg19_apprisplus.{strand}.bed'
        # need the genes file again in case any genes did not intersect with
        # the data
    threads: 12
    output:
        'output/metaplot/metabed/{peak_call_method}/{control}.vs._{treatment_a}_{treatment_b}.{mode}.{region}.{strand}.meta.bed'
    params:
        n_bins=100,
        out_dir='output/metaplot/metabed/{peak_call_method}'
    shell:'''
    mkdir -p {params.out_dir}
    python scripts/intersection_to_metabed.py {input.intersection} {input.genes} \
    {output} --threads {threads} --n_bins {params.n_bins}
    '''


rule make_metaplot:
    conda:
        '../envs/R.yml'
    input:
        'output/metaplot/metabed/{peak_call_method}/{control}.vs._{treatment_a}_{treatment_b}.{mode}.{region}.{strand}.meta.bed'
    output:
        'output/metaplot/plots/{peak_call_method}/{control}.vs._{treatment_a}_{treatment_b}.{mode}.{region}.{strand}.meta.png'
    params:
        out_dir='output/metaplot/plots/{peak_call_method}/'
    shell:'''
    mkdir -p {params.out_dir}
    Rscript scripts/metaplot.R {input} {output}
    '''
