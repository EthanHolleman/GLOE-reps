# rule call_peaks:
#     conda:
#         '../envs/macs2.yml'
#     input:
#     # output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.fwd.bed
#         treatment='output/{treatment}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed',
#         control='output/{control}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed'

#     output:
#         'output/call_peaks/{sample_1}.vs.{sample_2}_macs2/{mode}/{region}/{strand}'
#     params:
#         experiment_name='{treatment}.vs.{control}.{mode}.{strand}.{region}_macs2'
#     shell:'''
#     macs2 callpeak -t {treatment} -c {control} -n {params.experiment_name} \
#     --outdir {output} -m 5 50 -g 1.20E+07 --bw 200 --format BED --extsize 1 \
#     --nomodel --shift 0 --keep-dup
#     '''

# //vars for task macs2 from catalog ChIPseq, version 1
# MACS2_TARGETS="targets.txt" // targets file describing the samples
# MACS2_MFOLD="-m 5 50"	// range of enrichment ratio (default: 5,50)
# MACS2_GSIZE="-g " + ESSENTIAL_MACS2_GSIZE // the mappable genome size
# MACS2_BWIDTH="--bw " + Integer.toString(ESSENTIAL_FRAGLEN)	  // bandwidth use for model building
# MACS2_INPUT=BED // where the bed files are stored
# MACS2_FORMAT="--format BED"
# MACS2_EXTRA="--extsize 1 --nomodel --shift 0 --keep-dup " + ESSENTIAL_DUP	