
rule bed_to_coverage_bed:
    conda:
        '../envs/bedtools.yml'
    input:
        bed='output/{sample}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
        chrom_sizes='rawdata/hg19/hg19.chrom.sizes'
    output:
        'output/rfd/bed_coverage/{sample}.{mode}.{strand}.{region}.coverage.bedgraph'
    params:
        out_dir='output/rfd/bed_coverage'
    shell:'''
    mkdir -p {params.out_dir}
    bedtools genomecov -bg -i {input.bed} -trackline -g {input.chrom_sizes} \
    > {output}
    '''


rule sort_bedgraph:
    conda:
        '../envs/bedtools.yml'
    input:
        'output/rfd/bed_coverage/{sample}.{mode}.{strand}.{region}.coverage.bedgraph'
    output:
        temp('output/rfd/bedgraph/{sample}.{mode}.{strand}.{region}.sorted.coverage.bedgraph')
    shell:'''
    bedtools sort -i {input} > {output}
    '''


rule bed_to_bigwig:
    conda:
        '../envs/ucsc.yml'
    input:
        bedgraph='output/rfd/bedgraph/{sample}.{mode}.{strand}.{region}.sorted.coverage.bedgraph',
        chrom_sizes='rawdata/hg19/hg19.chrom.sizes'
    output:
        temp('output/rfd/bigwig/{sample}.{mode}.{strand}.{region}.bw')
    params:
        out_dir='output/rfd/bigwig'
    shell:'''
    mkdir -p {params.out_dir}
    bedGraphToBigWig {input.bedgraph} {input.chrom_sizes} {output}
    '''


rule bigwig_compare_strands_subtract:
    conda:
        '../envs/deeptools.yml'
    input:
        fwd_bigwig='output/rfd/bigwig/{sample}.{mode}.fwd.{region}.bw',
        rev_bigwig='output/rfd/bigwig/{sample}.{mode}.rev.{region}.bw'
    threads: 12
    output:
        'output/rfd/subtract/{sample}.{mode}.{region}.rev-fwd.bw'
    params:
        out_dir='output/rfd/subtract/'
    shell:'''
    mkdir -p {params.out_dir}
    bigwigCompare --numberOfProcessors {threads} -b1 {input.rev_bigwig} \
    -b2 {input.fwd_bigwig} -bs 1000 --operation subtract -of bigwig -o {output}
    '''


rule bigwig_compare_strands_add:
    conda:
        '../envs/deeptools.yml'
    input:
        fwd_bigwig='output/rfd/bigwig/{sample}.{mode}.fwd.{region}.bw',
        rev_bigwig='output/rfd/bigwig/{sample}.{mode}.rev.{region}.bw'
    threads: 12
    output:
        'output/rfd/add/{sample}.{mode}.{region}.rev+fwd.bw'
    params:
        out_dir='output/rfd/add'
    shell:'''
    mkdir -p {params.out_dir}
    bigwigCompare --numberOfProcessors {threads} -b1 {input.rev_bigwig} \
    -b2 {input.fwd_bigwig} -bs 1000 --operation add -of bigwig -o {output}
    '''


rule rfd:
    input:
        subtract='output/rfd/subtract/{sample}.{mode}.{region}.rev-fwd.bw',
        add='output/rfd/add/{sample}.{mode}.{region}.rev+fwd.bw'
    output:
        'output/rfd/rfd_bw/{sample}.{mode}.{region}.1000.rfd.bw'
    threads: 12
    shell:'''
    bigwigCompare --numberOfProcessors {threads} -b1 {input.subtract} \
    -b2 {input.add} -bs 1000 --pseudocount 0 --operation ratio -of bigwig \
    -o {output}
    '''