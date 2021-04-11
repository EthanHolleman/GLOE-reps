

rule seperate_footloop_strands:
    input:
        'rawdata/footloop/footloop_all.bed'
    output:
        fwd='rawdata/footloop/footloop_all.fwd.bed',
        rev='rawdata/footloop/footloop_all.rev.bed'
    shell:'''
    grep "+" {input} > {output.fwd}
    grep "-" {input} > {output.rev}
    '''


rule truncate_footloops:
    input:
        'rawdata/footloop/footloop_all.{strand}.bed'
    output:
        temp('output/truncated_footloop/footloop_all.{strand}.trunc.bed')
    shell:'''
    mkdir -p output/truncated_footloop/
    /home/ethollem/anaconda3/bin/Rscript scripts/truncate_footloops.R {input} {output}
    '''


rule sort_truncated_footloop_strands:
    input:
        'output/truncated_footloop/footloop_all.{strand}.trunc.bed'
    output:
        'output/truncated_footloop/footloop_all.{strand}.trunc.sorted.bed'
    shell:'''
    bedtools sort -i {input} > {output}
    '''


rule closest_footloop_swapped_mac_calls:
    # Get distance from footloop reads to peaks called using macs2 with
    # treatment and control input swapped in order to remove okazaki
    # fragment signal
    conda:
        '../envs/bedtools.yml'
    input:
        footloop='output/truncated_footloop/footloop_all.{strand}.trunc.sorted.bed',
        breaks='output/call_peaks/swapped/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/fwd/macs2_peak_call_summits.bed'
    output:
        'output/compare_footloop/closest/{control}.vs.{treatment_a}_{treatment_b}_macs2/footloop_closest.{mode}.{region}.{strand}.bed'
    params:
        out_dir='output/compare_footloop/closest/{control}.vs.{treatment_a}_{treatment_b}_macs2'
    shell:'''
    mkdir -p {params.out_dir}
    bedtools closest -t first -D a -a {input.footloop} -b {input.breaks} > {output}
    '''


rule plot_closest_footloop_swapped_mac_calls:
    conda:
        '../envs/R.yml'
    input:
        fwd='output/compare_footloop/closest/{control}.vs.{treatment_a}_{treatment_b}_macs2/footloop_closest.{mode}.{region}.fwd.bed',
        rev='output/compare_footloop/closest/{control}.vs.{treatment_a}_{treatment_b}_macs2/footloop_closest.{mode}.{region}.rev.bed'
    output:
        'output/call_peaks/swapped/{control}.vs.{treatment_a}_{treatment_b}_macs2/plots/{region}.{control}.vs.{treatment_a}_{treatment_b}.{mode}.png'
    shell:'''
    Rscript scripts/plt_called_peaks_footloop_dists.R {input.fwd} {input.rev} {output}
    '''