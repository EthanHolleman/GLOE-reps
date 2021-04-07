# Use macs to to call peaks between control and test samples



        #     BED=\$(basename $input) &&
        #     grep \$BED $MACS2_TARGETS | while read -r TARGET; do
        #         TREATMENT=\$(       echo \$TARGET | cut -f1 -d" ") &&
        #         Tname=\$(   echo \$TARGET | cut -f2 -d" ") &&
        #         CONTROL=\$(    echo \$TARGET | cut -f3 -d" ") &&
        #         Cname=\$(echo \$TARGET | cut -f4 -d" ") &&
        #         CompName=\$(echo \$TARGET | cut -f5 -d" ") &&
        #         TREATMENTFOR=\${TREATMENT%%.bed}".for.bed" && 
        #         TREATMENTREV=\${TREATMENT%%.bed}".rev.bed" && 
        #         TREATMENTnameFOR=\${Tname}"_for" && 
        #         TREATMENTnameREV=\${Tname}"_rev" && 
        #         CONTROLFOR=\${CONTROL%%.bed}".for.bed" && 
        #         CONTROLREV=\${CONTROL%%.bed}".rev.bed" && 
        #         CONTROLnameFOR=\${Cname}"_for" && 
        #         CONTROLnameREV=\${Cname}"_rev";   
                
        #         if [ "\$BED" != "\$INPUT" ]; then
        #             echo "\${Tname} vs \${Cname}" >> $output &&
        #             macs2 callpeak -t $MACS2_INPUT/\$TREATMENTFOR -c $MACS2_INPUT/\$CONTROLFOR -n \${TREATMENTnameFOR}.vs.\${CONTROLnameFOR}_macs2 $MACS2_FLAGS &&
        #             macs2 callpeak -t $MACS2_INPUT/\$TREATMENTREV -c $MACS2_INPUT/\$CONTROLREV -n \${TREATMENTnameREV}.vs.\${CONTROLnameREV}_macs2 $MACS2_FLAGS &&
        #             if [ \$? -ne 0 ]; then rm $output; fi &&
        #             awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" "+"}' \${TREATMENTnameFOR}.vs.\${CONTROLnameFOR}_macs2_summits.bed > \${CompName}_macs2.bed &&
        #             awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" "-"}' \${TREATMENTnameREV}.vs.\${CONTROLnameREV}_macs2_summits.bed >> \${CompName}_macs2.bed &&
        #             bedtools sort -i \${CompName}_macs2.bed > \${CompName}_macs2_summits.bed &&
        #             awk '{if(NR>20)print}' \${TREATMENTnameFOR}.vs.\${CONTROLnameFOR}_macs2_peaks.xls > \${CompName}_macs2_peaks.xls &&
        #             awk '{if(NR>21)print}' \${TREATMENTnameREV}.vs.\${CONTROLnameREV}_macs2_peaks.xls >> \${CompName}_macs2_peaks.xls &&
        #             cat \${TREATMENTnameFOR}.vs.\${CONTROLnameFOR}_macs2_peaks.narrowPeak \${TREATMENTnameREV}.vs.\${CONTROLnameREV}_macs2_peaks.narrowPeak > \${CompName}_macs2_peaks.narrowPeak && 
        #             rm \${TREATMENTnameFOR}.vs.\${CONTROLnameFOR}_macs2* &&
        #             rm \${TREATMENTnameREV}.vs.\${CONTROLnameREV}_macs2* &&
        #             rm \${CompName}_macs2.bed    &&
        #             mv \${CompName}_macs2* $output.dir ;
        #         fi;
        #     done
        # ""","macs2"

rule call_peaks:
    conda:
        '../envs/macs2.yml'
    input:
    # output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.fwd.bed
        treatment='output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.{strand}.bed',
        control='output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.{strand}.bed'
    output:
        dir('output/call_peaks/treatment_{treatment.sample}.control_{control.sample}_macs2/{mode}/{strand}/{region}')
    params:
        experiment_name='{treatment.sample}.vs.{control.sample}.{mode}.{strand}.{region}_macs2'
    shell:'''
    macs2 callpeak -t {input.treatment} -c {input.control} -n {params.experiment_name} --outdir {output}
    '''

# This would produce peak regions that we assume to be okazaki fragments so want
# want to go back and get non-spontatnous breaks
# Use the summit files to remove spontanous breaks and see where they occur in
# the genome that would be cools
