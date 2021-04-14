import pandas as pd

GLOE_SAMPLES = pd.read_csv(
    'samples/GLOE_samples.csv', sep=','
).set_index('Sample Name', drop=False)



rule index_bam:
    conda:
        '../envs/samtools.yml'
    input:
        'output/{sample}/process_alignment/sorted/trim.all.bam'
    output:
        'output/{sample}/process_alignment/sorted/trim.all.bai'
    shell:'''
    samtools index {input} {output}
    '''


rule sort_stranded_summits:
    conda: 
        '../envs/bedtools.yml'
    input:
        'output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.strand.bed'
    output:
        'output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.{region}.macs2_peak_call_summits.{strand}.all.sorted.bed'
    shell:'''
    bedtools sort -i {input} > {output} && [[ -s {output} ]]
    '''


# rule peak_called_bam_to_bedgraph:
#     input:
#         'output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.{region}.macs2_peak_call_summits.{strand}.all.sorted.bed'
#     output:
#         'output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.{region}.macs2_peak_call_summits.{strand}.all.sorted.bedgraph'
#     shell:'''
#     cut -f 1,2,3,5 {input} > {output} && [[ -s {output} ]]
#     '''


# rule peak_called_bedgraph_to_bigwig:
#     conda:
#         '../envs/ucsc.yml'
#     input:
#         bedgraph='output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.{region}.macs2_peak_call_summits.{strand}.all.sorted.bedgraph',
#         chromSizes='rawdata/hg19/hg19.chrom.sizes'
#     output:
#         'output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.{region}.macs2_peak_call_summits.{strand}.all.sorted.bw'
#     shell:'''
#     bedGraphToBigWig {input.bedgraph} {input.chromSizes} {output} && [[ -s {output} ]]
#     '''


rule multibigwigcompare:
    conda:
        '../envs/deeptools.yml'
    input:
        normal='output/call_peaks/normal/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.{region}.macs2_peak_call_summits.{strand}.all.sorted.bw',
        swapped='output/call_peaks/swapped/{sample_c}_{sample_d}.vs.{sample_a}_{sample_b}_macs2.{mode}.{region}.macs2_peak_call_summits.{strand}.all.sorted.bw'
    output:
        'output/fingerprint/{sample_a}_{sample_b}.{sample_c}_{sample_d}.{mode}.{region}.macs2_peak_call_summits.{strand}.all.npz'
    shell:'''
    mkdir -p output/fingerprint
    multiBigwigSummary bins --smartLabels -b {input.normal} {input.swapped} -o {output}
    '''


rule compute_matrix:
    conda:
        '../envs/deeptools.yml'
    input:
        normal='output/call_peaks/normal/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.{region}.macs2_peak_call_summits.{strand}.all.sorted.bw',
        swapped='output/call_peaks/swapped/{sample_c}_{sample_d}.vs.{sample_a}_{sample_b}_macs2.{mode}.{region}.macs2_peak_call_summits.{strand}.all.sorted.bw',
        genes='rawdata/hg19/hg19_apprisplus.{strand}.bed'
    threads: 12
    output:
        'output/call_peaks/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.{region}.macs2_peak_call_summits.{strand}.all.sorted.matrix.gz'
    shell:'''
    computeMatrix scale-regions -p {threads} -S {input.normal} {input.swapped} \
    -R {input.genes} -o {output}
    '''


rule compute_matrix_aligned_reads:
    conda:
        '../envs/deeptools.yml'
    input:
        bigwig='output/{sample}/reorient_alignments/{mode}/bigawk.sorted.trim.{region}.coverage.bw',
        genes='/home/ethollem/projects/GLOE-reps/rawdata/hg19/hg19_apprisplus.bed'
    threads: 2
    output:
        'output/fingerprint/matrix/{sample}.{mode}.{region}.genes.matrix.gz'
    shell:'''
    computeMatrix scale-regions -p {threads} -S {input.bigwig} \
    -R {input.genes} -o {output}
    '''


rule plot_profile_single_sample:
    conda:
        '../envs/deeptools.yml'
    input:
        'output/fingerprint/matrix/{sample}.{mode}.{region}.{profile_type}.matrix.gz'
    output:
        'output/fingerprints/plots/{sample}.{mode}.{region}.profile.{profile_type}.png'
    params:
        sample=lambda wildcards: wildcards['sample']
    shell:'''
    plotProfile -m {input} -o {output} --averageType mean --sampleLabel {params.sample}
    '''


rule coverage_reorriented_reads_over_genes:
    # calculate histograms of gene coverage for reorriented reads
    # doing this as the metaplots of called reads show much greater
    # signal in genes for revese calling 
    conda:
        '../envs/bedtools.yml'
    input:
        bed='output/{sample}/reorient_alignments/{mode}/bigawk.sorted.trim.{region}.bed',
        genes='rawdata/hg19/hg19_apprisplus.bed'
    output:
        hist='output/fingerprint/bed_coverage_genes/{sample}.{mode}.{region}.genes.coverage.hist.bed',
        sort_bed=temp('output/{sample}/reorient_alignments/{mode}/bigawk.sorted.trim.{region}.sorted.bed')
    params:
        out_dir = 'output/fingerprint/bed_coverage_genes'
    shell:'''
    mkdir -p {params.out_dir}
    bedtools sort -i {input.bed} > {output.sort_bed}
    bedtools coverage -sorted -hist -s -a {input.genes} -b {output.sort_bed}
    '''


rule plot_profile:
    conda:
        '../envs/deeptools.yml'
    input:
        'output/call_peaks/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.{region}.macs2_peak_call_summits.{strand}.all.sorted.matrix.gz'
    output:
        'output/fingerprint/plots/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}.{mode}.{region}.{strand}.profile.png'
    shell:'''
    plotProfile -m {input} -o {output} --averageType mean --samplesLabel normal reversed
    '''


rule plotCorrelation:
    conda:
        '../envs/deeptools.yml'
    input:
        'output/fingerprint/{sample_a}_{sample_b}.{sample_c}_{sample_d}.{mode}.{region}.macs2_peak_call_summits.{strand}.all.npz'
    output:
        'output/fingerprint/plots/{sample_a}_{sample_b}.{sample_c}_{sample_d}.{mode}.{region}.{strand}.correlation.png'
    params:
        cor_method='pearson',
        what_to_plot='scatterplot'
    shell:'''
    plotCorrelation --corData {input} --corMethod {params.cor_method} \
    --whatToPlot {params.what_to_plot} -o {output}
    '''

# rule fingerprint_samples:
#     conda:
#         '../envs/deeptools.yml'
#     input:
#         indexs=expand('output/{sample}/process_alignment/sorted/trim.all.bai', sample=GLOE_SAMPLES['Sample Name']),
#         bam=expand('output/{sample}/process_alignment/sorted/trim.all.bam', sample=GLOE_SAMPLES['Sample Name'])
#     output:
#         'output/fingerprint/sorted.trim.all.fingerprint.png'
#     threads: 12
#     shell:'''
#     plotFingerprint -p {threads} --smartLabels -b {input.bam} -o {output}
#     '''
