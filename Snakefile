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



wildcard_constraints:
    strand="fwd|rev",
    mode="indirect|direct",
    control="GSM[A-Za-z0-9]+"


import itertools

MODES = ['direct']
STRANDS = ['fwd', 'rev']
REGIONS = ['footloop', 'all']


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

macs2_peaks_swap = expand(
    'output/call_peaks/swapped/{control}.vs.{treatment}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.bed',
    treatment=treatment, control=control,
    mode=MODES, strand=STRANDS, region=REGIONS,
    allow_missing=True
)

ok_free_bed_files = expand(
    'output/{treatment}/ok_free/{mode}/{strand}/ok_free.control_{control}.{region}.bed',
    treatment=treatment, strand=STRANDS, region=REGIONS, mode=MODES,
    control=control
)

closest_footloops_swapped_peaks = expand(
    'output/compare_footloop/closest/{control}.vs.{treatment}_macs2/footloop_closest.{mode}.{strand}.bed',
    treatment=treatment, strand=STRANDS, region=REGIONS, mode=MODES,
    control=control
)

closest_footloops_swapped_peaks_plots = expand(
    'output/compare_footloop/closest/{control}.vs.{treatment}_macs2/plots/{control}.vs.{treatment}.{mode}.png',
    treatment=treatment, strand=STRANDS, region=REGIONS, mode=MODES,
    control=control
)




rule all:
    input:
        #macs2_peaks_swap,
        #macs2_peaks,
        #closest_footloops_swapped_peaks,
        closest_footloops_swapped_peaks_plots
        #ok_free_bed_files
        #"1"
        # expand(
        #     'output/{sample}/depth/{mode}/bam/{region}.{strand}.sorted.bam',
        #     sample=treatment, mode=MODES, region=REGIONS, strand=STRANDS
        # )



