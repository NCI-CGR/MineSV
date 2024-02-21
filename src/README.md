
# Installation of Truvari 2.0
```
git clone https://github.com/spiralgenetics/truvari.git
cd truvari/
python setup.py sdist bdist_wheel 
pip install pandas==1.0.3
pip install dist/Truvari-2.0.0.dev0.tar.gz
```

# Test run

```
python mocca_bench.py -b ../data/dream_modified.vcf.gz  -c ../data/IS1_intra.bed  -o test -f ../../refGenomes/Homo_sapiens_assembly19.fasta
371 857 260

```

# Svaba 
```bash

cd ~/CL/SV/mocca-bench/src

python mocca_bench.py -b ../data/dream_modified.vcf.gz  -c ../../test/IS1_4_new/output/compare_and_annotate/svaba/bed_input/IS1_intra.bed  -o test -f ../../refGenomes/Homo_sapiens_assembly19.fasta
371 782 304

### 1-base shift
2020-06-16 14:15:01,691 [INFO] Match: X 90613041        90628050        90613043        90628050        0.999867
2020-06-16 14:15:01,691 [INFO] Match: X 95229151        95233275        95229152        95233276        0.999758

# base vcf
X       90613042        .       T       <INV>   100     PASS    SOMATIC;SVTYPE=INV;END=90628050;SVLEN=15008     GT      ./.


# svaba
910 X__90613043__1048576886:1__T__T]X:90628050]__99__PASS__DISC_MAPQ=37;EVDNC=ASDIS;HOMSEQ=A;MAPQ=60;MATEID=1048576886:2;MATENM=0;NM=0;NUMPARTS=2;SCTG=c_23_90625501_90650501_19C;SPAN=15007;SVTYPE=BND__GT:AD:DP:GQ:PL:SR:DR:LR:LO__0/0:0:18:4.8:0,4.8,52.8:0:0:4.882:0__1/1:63:9:16.8:184.8,16.8,0:22:44:-184.8:184.8
911 X__90628050__1048576886:2__A__A]X:90613043]__99__PASS__DISC_MAPQ=37;EVDNC=ASDIS;HOMSEQ=A;MAPQ=60;MATEID=1048576886:1;MATENM=0;NM=0;NUMPARTS=2;SCTG=c_23_90625501_90650501_19C;SPAN=15007;SVTYPE=BND__GT:AD:DP:GQ:PL:SR:DR:LR:LO__0/0:0:18:4.8:0,4.8,52.8:0:0:4.882:0__1/1:63:9:16.8:184.8,16.8,0:22:44:-184.8:184.8
912 X__95229152__645788328:1__A__A[X:95233276[__99__PASS__DISC_MAPQ=31;EVDNC=ASDIS;MAPQ=60;MATEID=645788328:2;MATENM=0;NM=0;NUMPARTS=2;SCTG=c_23_95207001_95232001_2C;SPAN=4124;SVTYPE=BND__GT:AD:DP:GQ:PL:SR:DR:LR:LO__0/0:0:20:5.4:0,5.4,59.4:0:0:5.425:0__1/1:36:4:9.6:105.6,9.6,0:17:20:-105.6:105.6
913 X__95233276__645788328:2__C__]X:95229152]C__99__PASS__DISC_MAPQ=34;EVDNC=ASDIS;MAPQ=60;MATEID=645788328:1;MATENM=0;NM=0;NUMPARTS=2;SCTG=c_23_95207001_95232001_2C;SPAN=4124;SVTYPE=BND__GT:AD:DP:GQ:PL:SR:DR:LR:LO__0/0:0:20:5.4:0,5.4,59.4:0:0:5.425:0__1/1:36:4:9.6:105.6,9.6,0:17:20:-105.6:105.6

# breakdancer
python mocca_bench.py -b ../data/dream_modified.vcf.gz  -c ../../test/IS1_4_new/output/compare_and_annotate/breakdancer/bed_input/IS1_intra.bed  -o test -f ../../refGenomes/Homo_sapiens_assembly19.fasta
371 93411 124

# delly
python mocca_bench.py -b ../data/dream_modified.vcf.gz  -c ../../test/IS1_4_new/output/compare_and_annotate/delly/bed_input/IS1_intra.bed  -o test -f ../../refGenomes/Homo_sapiens_assembly19.fasta
371 857 260

# manta ?!
python mocca_bench.py -b ../data/dream_modified.vcf.gz  -c ../../test/IS1_4_new/output/compare_and_annotate/manta/bed_input/IS1_intra.bed  -o test -f ../../refGenomes/Homo_sapiens_assembly19.fasta
```

---
## Test on GIAB data
GIAB uses the same reference genome but its VCF is in different format, where type of SV is specified in SVTYPE in INFO.
```bash

awk '{if ($3=="hs37d5") print}'  ../../test/GIAB_germline/output/compare_and_annotate/svaba/bed_input/HG002_NA24385_son_intra.bed 
hs37d5  34796304        hs37d5  34796304        454453178:2     S       S]hs37d5:23832218]      42      LOWICSUPPORT      EVDNC=ASSMB;INSERTION=AAATAGAAGCGGCATCA;MAPQ=60;MATEID=454453178:1;MATENM=0;NM=0;NUMPARTS=2;SCTG=c_86_34790001_34815001_11C;SPAN=10964086;SVTYPE=BND      GT:AD:DP:GQ:PL:SR:DR:LR:LO      1/1:14:7:3.9:42.9,3.9,0:14:0:-42.91:42.91 152087  152087



python  mocca_bench.py -b ../data/HG002_SVs_Tier1_v0.6.vcf.gz  -c ../../test/GIAB_germline/output/compare_and_annotate/breakdancer/bed_input/HG002_NA24385_son_intra.bed   -o test -f ../../refGenomes/Homo_sapiens_assembly19.fasta --includebed ../data/HG002_SVs_Tier1_v0.6.bed --passonly

9641 1603 641

################# BWA realign ###################
### manta germline call from bwa realign (run by Bari)
python  mocca_bench.py -b ../data/HG002_SVs_Tier1_v0.6.vcf.gz  -c ../../test/2x250_bwa_germline_bari/compare_and_annotate/manta/bed_input/AJ_2x250_son_intra.bed   -o test -f ../../refGenomes/Homo_sapiens_assembly19.fasta --includebed ../data/HG002_SVs_Tier1_v0.6.bed --passonly --onebased
9641 3858 3079

mv ~/dn/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed ../data
python  mocca_bench.py -b ../data/HG002_SVs_Tier1_v0.6.vcf.gz  -c ../../test/2x250_bwa_germline_bari/compare_and_annotate/manta/bed_input/AJ_2x250_son_intra.bed   -o test -f ../../refGenomes/Homo_sapiens_assembly19.fasta --includebed ../data/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed --passonly --onebased


### svaba germline call from bwa realign (run by Bari)
python  mocca_bench.py -b ../data/HG002_SVs_Tier1_v0.6.vcf.gz  -c ../../test/2x250_bwa_germline_bari/compare_and_annotate/svaba/bed_input/AJ_2x250_son_intra.bed   -o test -f ../../refGenomes/Homo_sapiens_assembly19.fasta --includebed ../data/HG002_SVs_Tier1_v0.6.bed --passonly

9641 8360 870

### breakdancer germline call from bwa realign (run by Bari)
python  mocca_bench.py -b ../data/HG002_SVs_Tier1_v0.6.vcf.gz  -c ../../test/2x250_bwa_germline_bari/compare_and_annotate/breakdancer/bed_input/AJ_2x250_son_intra.bed   -o test -f ../../refGenomes/Homo_sapiens_assembly19.fasta --includebed ../data/HG002_SVs_Tier1_v0.6.bed --passonly  

9641 8307 1200
```

## Check some cases using suviz

https://svviz.readthedocs.io/en/latest/

https://github.com/nspies/svviz2

[:boom:] ERROR: Failed building wheel for seqlib
```bash
cd .. 

python3 -m venv ~/venv/svviz2_venv
source ~/venv/svviz2_venv/bin/activate
pip3 install -U git+git://github.com/nspies/svviz2.git

### failed with svviz2 and try with svviz1
pip install -U svviz


### run
### germline
### 2020-06-24 21:47:40,717 [INFO] Match: DEL       13      32176764        32176871        32176743        32176891        0.722973
-p 35429 
cd /home/zhuw10/CL/SV 

### failed 
module load R
pip install rpy2
svviz -t del -b ../GIAB/data/bam/HG002_BWA_2x250.bam  refGenomes/Homo_sapiens_assembly19.fasta 13      32176764        32176871 --dotplots -e  my.pdf
OSError: cannot load library '/DCEG/Resources/Tools/R/R-3.6.2/lib64/R/lib/libR.so': libpcre.so.1: cannot open shared object file: No such file or directory

# using web 
svviz -t del -b ../GIAB/data/bam/HG002_BWA_2x250.bam  refGenomes/Homo_sapiens_assembly19.fasta 13      32176764        32176871 -p 35429 


### check using truvari
truvari --base /home/english/data/giab_calls/HG002_SVs_Tier1_v0.6.vcf.gz --comp /home/english/data/ajtrio/manual_builder/37/results/HG002.bg/variants.vcf.gz --giabreport --passonly --includebed ../data/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed -o results_large_hc


truvari bench -p 0 -b ../data/HG002_SVs_Tier1_v0.6.vcf.gz -c ../data/HG002_SVs_Tier1_v0.6.vcf.gz -o TestvcfvsGIABv0.6 --passonly --includebed ../data/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed -r 2000 --giabreport


```


### debug: no INS 
```
python  mocca_bench.py -b ../data/HG002_SVs_Tier1_v0.6.vcf.gz  -c ../../test/2x250_bwa_germline_bari/compare_and_annotate/manta/bed_input/AJ_2x250_son_intra.bed   -o test -f ../../refGenomes/Homo_sapiens_assembly19.fasta --includebed ../data/HG002_SVs_Tier1_v0.6.bed --passonly --onebased  2> foo
9641 5554 3247

wc -l ../../test/2x250_bwa_germline_bari/compare_and_annotate/manta/bed_input/AJ_2x250_son_intra.bed
10963 ../../test/2x250_bwa_germline_bari/compare_and_annotate/manta/bed_input/AJ_2x250_son_intra.bed

module load bedtools

bedtools intersect -a ../data/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed -b ../../test/2x250_bwa_germline_bari/compare_and_annotate/manta/bed_input/AJ_2x250_son_intra.bed -f 1e-10 -F 1 |wc -l 
5417

bedtools intersect -a ../data/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed -b ../../test/2x250_bwa_germline_bari/compare_and_annotate/manta/bed_input/AJ_2x250_son_intra.bed -f 1e-10 -F 1 > manta_germline.bed
```

### manually check whether there is any overlap between INS and the prediction
```
bedtools intersect -a manta_germline.bed -b ../data/HG002_SVs_Tier1_v0.6.vcf.gz -wa

conda install -c bioconda bedops
zcat ../data/HG002_SVs_Tier1_v0.6.vcf.gz | dos2unix > giab.vcf
vcf2bed   <giab.vcf  > giab.bed
Error: Could not find newline in intermediate buffer; check input [34634 | 6475 | 41109]
       Please check that your input contains Unix newlines (cat -A) or increase TOKENS_MAX_LENGTH in BEDOPS.Constants.hpp and recompile BEDOPS.
```

### convert giab.vcf to bed using SUVIVOR
```
grep PASS giab.vcf > giab_pass.vcf
SURVIVOR vcftobed  giab_pass.vcf 50 50000  giab_flt.bed

```

### has PASS selected? NO and the coordinates are confusing: I need drop the first 6 columns and print 12,13,15 as the first 
```
vcf:
1       421058  HG2_Ill_SpiralSDKrefine_2       GAGCATCCTTGCAGCTGCAGGCTTCAGTCTACCAGAGAATGTGAGGTGTTATTCTTCTAGGGCAGTGGTTAGAAAAGAAAATGAAAGTAGCAGTACTCTTTTCCTAATGCAACCATAGATGGATGATCAGAATTTGTAATCCATAAGGTAGAAG      A
       20      LongReadHomRef  ClusterIDs=HG2_Ill_SpiralSDKrefine_2:HG2_PB_PB10Xdip_2158;NumClusterSVs=2;ExactMatchIDs=HG2_Ill_SpiralSDKrefine_2;NumExactMatchSVs=1;ClusterMaxShiftDist=0;ClusterMaxSizeDiff=0;ClusterMaxEditDist=0.00649350649350649;PBcalls=1;Illcalls=1;TenXcalls=0;CGcalls=0;PBexactcalls=0;Illexactcalls=1;TenXexactcalls=0;CGexactcalls=0;HG2count=2;HG3count=0;HG4count=0;NumTechs=2;NumTechsExact=1;SVLEN=-153;DistBack=34378;DistForward=15514;DistMin=15514;DistMinlt1000=FALSE;MultiTech=TRUE;MultiTechExact=FALSE;SVTYPE=DEL;END=421211;sizecat=100to299;DistPASSHG2gt49Minlt1000=.;DistPASSMinlt1000=.;MendelianError=.;HG003_GT=./.;HG004_GT=./.;TRall=FALSE;TRgt100=FALSE;TRgt10k=FALSE;segdup=TRUE;REPTYPE=SUBSDEL;BREAKSIMLENGTH=0;REFWIDENED=1:421058-421211 GT:GTcons1:PB_GT:PB_REF:PB_ALT:PBHP_GT:PB_REF_HP1:PB_ALT_HP1:PB_REF_HP2:PB_ALT_HP2:TenX_GT:TenX_REF_HP1:TenX_ALT_HP1:TenX_REF_HP2:TenX_ALT_HP2:ILL250bp_GT:ILL250bp_REF:ILL250bp_ALT:ILLMP_GT:ILLMP_REF:ILLMP_ALT:BNG_LEN_DEL:BNG_LEN_INS:nabsys_svm    ./.:./.:./.:0:0:./.:0:0:0:0:./.:0:0:0:0:./.:0:0:./.:.:.:.:.:.   

converted bed: 
1       421056  421056  1       421210  421210  HG2_Ill_SpiralSDKrefine_2       ,       +       -       DEL     1       421057  1       421211


awk '{if($7=="PASS") print}' giab.vcf > giab_pass.vcf
wc -l giab_pass.vcf
12745 giab_pass.vcf

SURVIVOR vcftobed  giab_pass.vcf 50 50000  giab_flt.bed
wc -l giab.bed
29421 giab.bed

wc -l giab_flt.bed
12503 giab_flt.bed

awk -v OFS='\t' '{print $12,$13,$15,$7,$8,$9,$10,$11,$12,$13,$14,$15}' giab_flt.bed > giab_flt2.bed

bedtools intersect -a giab_flt2.bed -b ../data/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed -f 1 -wa  > giab_v6.bed

wc -l giab_v6.bed 
9185 giab_v6.bed

```
### Now find out the overlap between giab_v6.bed and prediction
```
bedtools intersect -a giab_flt2.bed -b manta_germline.bed -wa -wb > foo.bed
wc -l foo.bed
3302

cut -f 8 foo.bed | sort | uniq -c
   3069 DEL
    233 INS

For example:
1       3748983 3748984 HG2_PB_HySA_172 ,       +       -       INS     1       3748983 1       3748984 VS 1       3748936 3749008
1       3748984 HG2_PB_HySA_172 G       GAACCAGGAGTGTAGTGTAGGTCGTGACCAGGAGTGTAGTGTAGGTCATAACCAGGAGTGTAGTGTAGCACATTAGGTCGTA      20      PASS    ClusterIDs=HG2_PB_HySA_172:HG3_PB_pbsv_108:HG2_Ill_250bpfermikitraw_78:HG2_PB_pbsv_106:HG2_PB_PB10Xdip_2140:HG2_PB_PB10Xdip_2139:HG4_PB_HySA_144;NumClusterSVs=7;ExactMatchIDs=HG2_PB_HySA_172:HG4_PB_HySA_144;NumExactMatchSVs=2;ClusterMaxShiftDist=0.0952380952380952;ClusterMaxSizeDiff=0.0952380952380952;ClusterMaxEditDist=0.300578034682081;PBcalls=6;Illcalls=1;TenXcalls=0;CGcalls=0;PBexactcalls=2;Illexactcalls=0;TenXexactcalls=0;CGexactcalls=0;HG2count=5;HG3count=1;HG4count=1;NumTechs=2;NumTechsExact=1;SVLEN=81;DistBack=-25;DistForward=37;DistMin=-25;DistMinlt1000=TRUE;MultiTech=TRUE;MultiTechExact=FALSE;SVTYPE=INS;END=3748984;sizecat=50to99;DistPASSHG2gt49Minlt1000=FALSE;DistPASSMinlt1000=TRUE;MendelianError=FALSE;HG003_GT=0/1;HG004_GT=0/1;TRall=TRUE;TRgt100=TRUE;TRgt10k=FALSE;segdup=FALSE;REPTYPE=DUP;BREAKSIMLENGTH=1;REFWIDENED=1:3748983-3748986       GT:GTcons1:PB_GT:PB_REF:PB_ALT:PBHP_GT:PB_REF_HP1:PB_ALT_HP1:PB_REF_HP2:PB_ALT_HP2:TenX_GT:TenX_REF_HP1:TenX_ALT_HP1:TenX_REF_HP2:TenX_ALT_HP2:ILL250bp_GT:ILL250bp_REF:ILL250bp_ALT:ILLMP_GT:ILLMP_REF:ILLMP_ALT:BNG_LEN_DEL:BNG_LEN_INS:nabsys_svm    1/1:1/1:1/1:4:45:1/1:4:30:1:9:./.:5:2:4:7:0/1:4:12:./.:.:.:.:.:.
```

### Re-test all 


```
python my_vcf2bed.py -i ../../bench_ws/giab/tp-base.vcf > my-tp-base.bed

python my_vcf2bed.py -i ../../bench_ws/synthetic/is1.vcf  | grep -P "DUP|INV|DEL" > my_is1.bed

wc -l *.bed
    3302 foo.bed
   29421 giab.bed
   12503 giab_flt2.bed
   12503 giab_flt.bed
   29421 giab_good.bed
    9185 giab_v6.bed
    5417 manta_germline.bed
     371 my_is1.bed
    9361 my-tp-base.bed
  111484 total
```