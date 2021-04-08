rule call_peaks:
    conda:
        '../envs/macs2.yml'
    input:
    # output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.fwd.bed
        treatment='output/{treatment}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
        control='output/{control}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed'

    output:
        'output/call_peaks/normal/{treatment}.vs.{control}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.bed'
    params:
        experiment_name='macs2_peak_call',
        out_dir='output/call_peaks/normal/{treatment}.vs.{control}_macs2/{mode}/{region}/{strand}'
    shell:'''
    mkdir -p {params.out_dir}
    macs2 callpeak -t {input.treatment} -c {input.control} -n {params.experiment_name} \
    --outdir {params.out_dir} -m 5 50 -g 1.20E+07 --bw 200 --format BED --extsize 1 \
    --nomodel --shift 0 --keep-dup all
    '''

rule call_peaks_swapped:
    conda:
        '../envs/macs2.yml'
    input:
    # output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.fwd.bed
        treatment='output/{treatment}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
        control='output/{control}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed'

    output:
        'output/call_peaks/swapped/{control}.vs.{treatment}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.bed'
    params:
        experiment_name='macs2_peak_call',
        out_dir='output/call_peaks/swapped/{control}.vs.{treatment}_macs2/{mode}/{region}/{strand}'
    shell:'''
    mkdir -p {params.out_dir}
    macs2 callpeak -t {input.control} -c {input.treatment} -n {params.experiment_name} \
    --outdir {params.out_dir} -m 5 50 -g 1.20E+07 --bw 200 --format BED --extsize 1 \
    --nomodel --shift 0 --keep-dup all
    '''

rule remove_called_peaks:
    input:
        peaks='output/call_peaks/{treatment}.vs.{control}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.bed',
        treatment_bed='output/{treatment}/depth/{mode}/deep/{region}.{strand}.sorted.bed'
    output:
        'output/{treatment}/ok_free/{mode}/{strand}/ok_free.control_{control}.{region}.bed'
    params:
        out_dir='output/{treatment}/ok_free/{mode}/{strand}'
    shell:'''
    mkdir -p {params.out_dir}
    bedtools intersect -v -a {input.treatment_bed} -b {input.peaks} > {output}
    '''

    # Why fos and jun were targeted? These 