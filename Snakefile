include: "rules/call_peaks.smk"
include: "rules/download_data.smk"
include: "rules/bowtie2.smk"
include: "rules/trimmomatic.smk"
include: "rules/process_alignment.smk"
include: "rules/compare_replicates.smk"
include: "rules/reorient_aligments.smk"
include: "rules/depth.smk"
include: "rules/remove_okazaki.smk"
include: "rules/compare_footloop.smk"
include: "rules/metaplot.smk"
include: "rules/rfd.smk"



wildcard_constraints:
    strand="fwd|rev",
    mode="indirect|direct",
    control="GSM[A-Za-z0-9]+"


import itertools

MODES = ['direct']
STRANDS = ['fwd', 'rev']
REGIONS = ['footloop', 'all']
PEAK_CALL_METHODS = ['normal', 'swapped']

# swapped peak calling is calling peaks with the control as the treatment

def group_replicates(samples):
    all_samples = []
    for _, sample in samples.iterrows():
        modification = sample['modification']
        replicates = samples.loc[samples['modification'] == modification]
        replicates = tuple([rep['Sample Name'] for _, rep in replicates.iterrows()])
        all_samples.append(replicates)
    return all_samples


replicate_groups = group_replicates(GLOE_SAMPLES)
paired_replicates_boxplot_output = expand(
    'output/compare_replicates/{rep[0]}-{rep[1]}/plots/closest/{mode}/{strand}/boxplot_{region}.png',
    rep=replicate_groups, allow_missing=True
)

treatment=list(GLOE_SAMPLES.loc[GLOE_SAMPLES['modification'] == 'siRNA LIG1']['Sample Name'])
control=list(GLOE_SAMPLES.loc[GLOE_SAMPLES['modification'] == 'siRNA control']['Sample Name'])
macs2_peaks = expand(
    'output/call_peaks/normal/{treatment}.vs.{control}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.bed',
    treatment=treatment, control=control,
    mode=MODES, strand=STRANDS, region=REGIONS,
    allow_missing=True
)


metaplot_metabed_files = expand(
    'output/metaplot/metabed/{peak_call_method}/{control}.vs._{treatment_a}_{treatment_b}.{mode}.all.{strand}.meta.bed',
    treatment_a=treatment, treatment_b=treatment, strand=STRANDS, mode=MODES,
    control=control, peak_call_method=PEAK_CALL_METHODS

)

metaplots = expand(
    'output/metaplot/plots/{peak_call_method}/{control}.vs._{treatment_a}_{treatment_b}.{mode}.all.{strand}.meta.png',
    treatment_a=treatment, treatment_b=treatment, strand=STRANDS, mode=MODES,
    control=control, peak_call_method=PEAK_CALL_METHODS, 
)

rfd = expand(
    'output/rfd/rfd_bw/{sample}.{mode}.all.1000.rfd.bw',
    sample=treatment, mode=MODES
)

all_peaks_stranded_average_normal = expand(
    'output/call_peaks/normal/{treatment_a}_{treatment_b}.vs.{control_a}_{control_b}_macs2/{mode}/all/macs2_peak_call_summits.strand.all.bed',
    zip, treatment_a=treatment, treatment_b=treatment, mode=MODES,
    control_a=control, control_b=control
)

all_peaks_stranded_average_normal = expand(
    'output/call_peaks/normal/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.all.macs2_peak_call_summits.strand.all.bed',
    zip, sample_a=treatment, sample_b=treatment[::-1], mode=MODES,
    sample_d=control[::-1], sample_c=control
)

all_peaks_stranded_average_swapped = expand(
    'output/call_peaks/swapped/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.all.macs2_peak_call_summits.strand.all.bed',
    zip, sample_a=control, sample_b=control[::-1], mode=MODES,
    sample_d=treatment[::-1], sample_c=treatment
)

all_peak_stranded_average_plots = expand(
    'output/call_peaks/plots/{sample_c}_{sample_d}.vs.{sample_a}_{sample_b}.{mode}.all.anno.all.png',
     zip, sample_a=control, sample_b=control[::-1], mode=MODES,
    sample_d=treatment[::-1], sample_c=treatment
)


rule all:
    input:
        #macs2_peaks_swap,
        #macs2_peaks,
        #closest_footloops_swapped_peaks,
        #closest_footloops_swapped_peaks_plots
        #metaplot_intersection_files
        #metaplot_metabed_files
        #metaplots
        #rfd
        #all_peaks_stranded
        #plot_annotations,
        #test_chipseeker,
        all_peaks_stranded_average_swapped,
        all_peaks_stranded_average_normal,
        all_peak_stranded_average_plots
        #closest_footloops_swapped_peaks_multi_treatment
        #ok_free_bed_files
        #"1"
        # expand(
        #     'output/{sample}/depth/{mode}/bam/{region}.{strand}.sorted.bam',
        #     sample=treatment, mode=MODES, region=REGIONS, strand=STRANDS
        # )



