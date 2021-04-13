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


rule call_peaks_swapped_average_sample:
    conda:
        '../envs/macs2.yml'
    input:
    # output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.fwd.bed
        treatment_a='output/{treatment_a}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
        treatment_b='output/{treatment_b}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
        control_a='output/{control_a}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
        control_b='output/{control_b}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed'
    output:
        'output/call_peaks/swapped/{control_a}_{control_b}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.bed'
    params:
        experiment_name='macs2_peak_call',
        out_dir='output/call_peaks/swapped/{control_a}_{control_b}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/{strand}'
    shell:'''
    mkdir -p {params.out_dir}
    macs2 callpeak -t {input.control_a} {input.control_b} -c {input.treatment_a} {input.treatment_b} -n {params.experiment_name} \
    --outdir {params.out_dir} -m 5 50 -g 1.20E+07 --bw 200 --format BED --extsize 1 \
    --nomodel --shift 0 --keep-dup all
    '''


rule call_peaks_normal_average_sample:
    conda:
        '../envs/macs2.yml'
    input:
    # output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.fwd.bed
        treatment_a='output/{treatment_a}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
        treatment_b='output/{treatment_b}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
        control_a='output/{control_a}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
        control_b='output/{control_b}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed'
    output:
        'output/call_peaks/normal/{treatment_a}_{treatment_b}.vs.{control_a}_{control_b}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.bed'
    params:
        experiment_name='macs2_peak_call',
        out_dir='output/call_peaks/normal/{treatment_a}_{treatment_b}.vs.{control_a}_{control_b}_macs2/{mode}/{region}/{strand}'
    shell:'''
    mkdir -p {params.out_dir}
    macs2 callpeak -t {input.treatment_a} {input.treatment_b} -c {input.control_a} {input.control_b} -n {params.experiment_name} \
    --outdir {params.out_dir} -m 5 50 -g 1.20E+07 --bw 200 --format BED --extsize 1 \
    --nomodel --shift 0 --keep-dup all
    '''


rule add_strand_to_summits_fwd_average_sample:
    input:
        'output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2/{mode}/{region}/fwd/macs2_peak_call_summits.bed'
    output:
        temp('output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2/{mode}/{region}/fwd/macs2_peak_call_summits.strand.bed')
    shell:'''
    awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""+"}}' {input} > {output}
    '''


rule add_strand_to_summits_rev_average_sample:
    input:
        'output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2/{mode}/{region}/rev/macs2_peak_call_summits.bed'
    output:
        temp('output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2/{mode}/{region}/rev/macs2_peak_call_summits.strand.bed')
    shell:'''
    awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""-"}}' {input} > {output}
    '''


rule peaks_stranded_all_average_sample:
    input:
        fwd='output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2/{mode}/{region}/fwd/macs2_peak_call_summits.strand.bed',
        rev='output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2/{mode}/{region}/rev/macs2_peak_call_summits.strand.bed'
    output:
        'output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.{region}.macs2_peak_call_summits.strand.all.bed'
    shell:'''
    cat {input.fwd} {input.rev} > {output}
    '''


rule course_annotate_peaks_average_sample:
    input:
        'output/call_peaks/{peak_call_method}/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.{region}.macs2_peak_call_summits.strand.all.bed'
    output:
        'output/call_peaks/{peak_call_method}/annotate/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}.{mode}.{region}.anno.all.tsv'
    params:
        out_dir='output/call_peaks/{peak_call_method}/annotate'
    shell:'''
    mkdir -p {params.out_dir}
    /home/ethollem/anaconda3/envs/test/bin/Rscript scripts/annotate_peaks.R {input} {output}
    '''


rule plot_course_annotation_average_sample:
    conda:
        '../envs/R.yml'
    input:
        normal='output/call_peaks/normal/annotate/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}.{mode}.{region}.anno.all.tsv',
        swapped='output/call_peaks/swapped/annotate/{sample_c}_{sample_d}.vs.{sample_a}_{sample_b}.{mode}.{region}.anno.all.tsv'
    output:
        'output/call_peaks/plots/{sample_a}_{sample_b}.{sample_c}_{sample_d}.{mode}.{region}.anno.all.png'
    params:
        out_dir='output/call_peaks/plots'
    shell:'''
    mkdir -p {params.out_dir}
    Rscript scripts/plot_annotations.R {input.normal} {input.swapped} {output}
    '''


rule add_strand_to_summits_fwd:
    input:
        'output/call_peaks/{peak_call_method}/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/fwd/macs2_peak_call_summits.bed'
    output:
        temp('output/call_peaks/{peak_call_method}/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/fwd/macs2_peak_call_summits.strand.bed')
    shell:'''
    awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""+"}}' {input} > {output}
    '''


rule add_strand_to_summits_rev:
    input:
        'output/call_peaks/{peak_call_method}/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/rev/macs2_peak_call_summits.bed'
    output:
        temp('output/call_peaks/{peak_call_method}/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/rev/macs2_peak_call_summits.strand.bed')
    shell:'''
    awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""-"}}' {input} > {output}
    '''


rule peaks_stranded_all:
    input:
        fwd='output/call_peaks/{peak_call_method}/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/fwd/macs2_peak_call_summits.strand.bed',
        rev='output/call_peaks/{peak_call_method}/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/rev/macs2_peak_call_summits.strand.bed'
    output:
        'output/call_peaks/{peak_call_method}/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/macs2_peak_call_summits.strand.all.bed'
    shell:'''
    cat {input.fwd} {input.rev} > {output}
    '''


rule course_annotate_peaks:
    input:
        'output/call_peaks/{peak_call_method}/{control}.vs.{treatment_a}_{treatment_b}_macs2/{mode}/{region}/macs2_peak_call_summits.strand.all.bed'
    output:
        'output/call_peaks/{peak_call_method}/annotate/{control}.vs.{treatment_a}_{treatment_b}_{mode}.{region}.anno.all.tsv'
    params:
        out_dir='output/call_peaks/{peak_call_method}/annotate'
    shell:'''
    mkdir -p {params.out_dir}
    /home/ethollem/anaconda3/envs/test/bin/Rscript scripts/annotate_peaks.R {input} {output}
    '''



# rm /home/ethollem/projects/GLOE-reps/output/call_peaks/swapped/GSM3939125.vs.GSM4305465_GSM4305466_macs2/direct/all/macs2_peak_call_summits.strand.all.bed /home/ethollem/projects/GLOE-reps/output/call_peaks/swapped/GSM3939125.vs.GSM4305466_GSM4305465_macs2/direct/all/macs2_peak_call_summits.strand.all.bed /home/ethollem/projects/GLOE-reps/output/call_peaks/swapped/GSM3939125.vs.GSM4305466_GSM4305466_macs2/direct/all/macs2_peak_call_summits.strand.all.bed /home/ethollem/projects/GLOE-reps/output/call_peaks/swapped/GSM3939126.vs.GSM4305465_GSM4305465_macs2/direct/all/macs2_peak_call_summits.strand.all.bed /home/ethollem/projects/GLOE-reps/output/call_peaks/swapped/GSM3939126.vs.GSM4305465_GSM4305466_macs2/direct/all/macs2_peak_call_summits.strand.all.bed /home/ethollem/projects/GLOE-reps/output/call_peaks/swapped/GSM3939126.vs.GSM4305466_GSM4305465_macs2/direct/all/macs2_peak_call_summits.strand.all.bed /home/ethollem/projects/GLOE-reps/output/call_peaks/swapped/GSM3939126.vs.GSM4305466_GSM4305466_macs2/direct/all/macs2_peak_call_summits.strand.all.bed



