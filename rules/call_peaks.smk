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
        treatment_a='output/{treatment_a}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
        treatment_b='output/{treatment_b}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
        control='output/{control}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed'
    output:
        'output/call_peaks/swapped/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.bed'
    params:
        experiment_name='macs2_peak_call',
        out_dir='output/call_peaks/swapped/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/{strand}'
    shell:'''
    mkdir -p {params.out_dir}
    macs2 callpeak -t {input.control} -c {input.treatment_a} {input.treatment_b} -n {params.experiment_name} \
    --outdir {params.out_dir} -m 5 50 -g 1.20E+07 --bw 200 --format BED --extsize 1 \
    --nomodel --shift 0 --keep-dup all
    '''