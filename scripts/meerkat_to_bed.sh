#!/bin/sh

input=$1
sample=$2
outDir=$3
#varType="SPAN=-1"
varType="'EVENT=(del_inss|del_inso|inss|inso|transl_inter)_'"

# escape all backslashes, quotes, and dollar signs that are not within square brackets

cmd_add_linenum="awk '{FS=OFS=\"\\t\"}{print \$0, NR}'"
cmd_remove_headers="grep -v \"^#\""
cmd_find_intraSV="grep -v -E $varType"
cmd_find_interSV="grep -E $varType"

# C]21:48119850]
cmd_get_intra_coords="awk '{FS=OFS=\"\\t\"} {match(\$5, \"[\\\\[\\\\]]([^:]+):([0-9]+)[\\\\[\\\\]]\", a); chr=a[1];pos1=a[2]} {print \$1, \$2, pos1, \$NF}'"

cmd_sort_intra_coords="awk '{FS=OFS=\"\\t\"} (\$2<\$3){print \$0} (\$2>\$3){print \$1,\$3,\$2,\$4} (\$2==\$3){print \$1,\$2,\$3+1,\$4}'"
cmd_get_end1_coords="awk '{FS=OFS=\"\\t\"} (\$2!=0){pos1=\$2-1; pos2=\$2} (\$2==0){pos1=\$2; pos2=\$2+1} {print \$1,pos1,pos2,\$NF}'"

# C]21:48119850]
cmd_get_end2_coords="awk '{FS=OFS=\"\\t\"} {match(\$5, \"[\\\\[\\\\]]([^:]+):([0-9]+)[\\\\[\\\\]]\", a); chr=a[1];pos1=a[2]} (pos1!=0){print chr,pos1-1,pos1,\$NF} (pos1==0){print chr,pos1,pos1+1,\$NF}'"

if [[ $input = *.gz ]]; then
    cmd_intra="zcat ${input} | ${cmd_add_linenum} | ${cmd_remove_headers} | ${cmd_find_intraSV} | ${cmd_get_intra_coords} | ${cmd_sort_intra_coords} > ${outDir}${sample}_intra.bed"
    cmd_end1="zcat ${input} | ${cmd_add_linenum} | ${cmd_remove_headers} | ${cmd_find_interSV} | ${cmd_get_end1_coords} > ${outDir}${sample}_end1.bed"
    cmd_end2="zcat ${input} | ${cmd_add_linenum} | ${cmd_remove_headers} | ${cmd_find_interSV} | ${cmd_get_end2_coords} > ${outDir}${sample}_end2.bed"
else
    cmd_intra="${cmd_add_linenum} ${input} | ${cmd_remove_headers} | ${cmd_find_intraSV} | ${cmd_get_intra_coords} | ${cmd_sort_intra_coords} > ${outDir}${sample}_intra.bed"
    cmd_end1="${cmd_add_linenum} ${input} | ${cmd_remove_headers} | ${cmd_find_interSV} | ${cmd_get_end1_coords} > ${outDir}${sample}_end1.bed"
    cmd_end2="${cmd_add_linenum} ${input} | ${cmd_remove_headers} | ${cmd_find_interSV} | ${cmd_get_end2_coords} > ${outDir}${sample}_end2.bed"
fi

echo $cmd_intra
eval $cmd_intra
echo $cmd_end1
eval $cmd_end1
echo $cmd_end2
eval $cmd_end2
