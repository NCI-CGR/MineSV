#!/usr/bin/env perl


@ARGV == 3 or die "$0 <INPUT VCF file> <Sample ID> <output Dir> \n";  

my $input=shift;
my $sample = shift;
my $output_dir = shift;

##################################################
# To convert GRIDSS V4.2 VCF file to bed format 
##################################################

# Three files to be outputed here:
# ${outDir}${sample}_intra.bed: chr, start, end and line ID; start should be less than end. if start==end, end=start+1.
# ${outDir}${sample}_end1.bed
# ${outDir}${sample}_end2.bed 

### Sample input
# zcat ../../call_germline/output/HG002_NA24385_son.vcf.gz | less
# 1       54712   gridss0fb_21o   T       TTTTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTC[1:54713[ 706.37
#   SINGLE_ASSEMBLY AS=0;ASC=1X;ASQ=0.00;ASRP=1;ASSR=22;BA=0;BANRP=0;BANRPQ=0.00;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BEALN=22:36615720|-|53M|0;BEID=asm0-80445,asm0-80446;BEIDH=0,0;BEIDL=0,0;BMQ=70.00;BMQN=70.00;BMQX=70.00;BQ=16.53;BSC=1;BSCQ=16.53;BUM=0;BUMQ=0.00;BVF=0;CAS=0;CASQ=0.00;CQ=706.37;EVENT=gridss0fb_21;IC=0;IHOMPOS=0,0;IQ=0.00;MATEID=gridss0fb_21h;MQ=60.00;MQN=60.00;MQX=60.00;RAS=2;RASQ=706.37;REF=21;REFPAIR=0;RP=0;RPQ=0.00;SB=0.4375;SC=1X;SR=0;SRQ=0.00;SVTYPE=BND;VF=14      GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF   .:0.00:1:22:0:0.00:0:0.00:0.00:0:0:16.53:1:16.53:0:0.00:0:0.00:0:0.00:706.37:706.37:21:0:0:0.00:0:0.00:14
# 1       54713   gridss0fb_21h   T       ]1:54712]TTTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCT 706.37
#   SINGLE_ASSEMBLY AS=2;ASC=1X;ASQ=706.37;ASRP=1;ASSR=22;BA=0;BANRP=0;BANRPQ=0.00;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BEALN=GL000220.1:119381|+|53M|0;BEID=asm0-80445,asm0-80446;BEIDH=0,0;BEIDL=0,0;BMQ=70.00;BMQN=70.00;BMQX=70.00;BQ=163.31;BSC=6;BSCQ=107.62;BUM=3;BUMQ=55.69;BVF=6;CAS=0;CASQ=0.00;CQ=706.37;EVENT=gridss0fb_21;IC=0;IHOMPOS=0,0;IQ=0.00;MATEID=gridss0fb_21o;MQ=60.00;MQN=60.00;MQX=60.00;RAS=0;RASQ=0.00;REF=21;REFPAIR=0;RP=0;RPQ=0.00;SB=0.43243244;SC=1X;SR=0;SRQ=0.00;SVTYPE=BND;VF=14 GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF   .:706.37:1:22:0:0.00:0:0.00:0.00:0:0:163.31:6:107.62:3:55.69:6:0.00:0:0.00:706.37:0.00:21:0:0:0.00:0:0.00:14

##################################################
# Essential information to capture:
# - line number: line number of the current file handle is stored in the special variable $.
# - chr1, pos1, chr2, pos2
#  
#  chr2:pos2 could be obtained from the column 
# We may need ignore short indels:
# 1       104160  gridss0fb_55o   A       AACACACACAC[1:104161[ 
 

# open file handles
my $basefn = $output_dir.'/'.$sample;

$input = "zcat $input |" if ($input =~ 'gz$'); 

open(FIN, $input) or die $!; 
open(INTRA, '>', $basefn."_intra.bed") or die $!; 
open(INTER1, '>', $basefn."_end1.bed") or die $!;
open(INTER2, '>', $basefn."_end2.bed") or die $!;

my $i=0;
while(<FIN>){
    next if /^#/;
    chomp;
    my ($chr1, $pos1, $id, $ref, $alt, $qual, $filter, $info, $format) = split(/\t/,$_);
    
    # next if $id =~ /b$/; # ignore single BND
    next if $id !~ /o$/ | $filter ne "PASS"; # ignore the parter BND and focus "o"

    my $NR = $.; 
    
    my %info_hash = map { /([^=]+)=?(\S*)/ } split m{;}, $info; 

    my ($left, $orient, $chr2, $pos2,$right) = ($alt =~ /^(.*)([\[\]])(.+):(.+)\2(.*)$/);

    if($chr1 eq $chr2){
        # intrachromosomal
        if($pos1 <= $pos2){
            print INTRA join("\t", $chr1, $pos1, $pos2, $NR)."\n";
        }else{
            print INTRA join("\t", $chr1, $pos2, $pos1, $NR)."\n";
        }
    }else{
        # interchromosomal
        $pos1 = ($pos1==0)?0:$pos1-1; 
        $pos2 = ($pos2==0)?0:$pos2-1; 
        print INTER1 join("\t", $chr1, $pos1, $pos1+1, $NR)."\n";
        print INTER2 join("\t", $chr2, $pos2, $pos2+1, $NR)."\n";
    }


    # print join("\t", $chr1, $pos1, $id, $ref, $alt, $info_hash{'EVENT'}, $NR, $left, $orient, $chr2, $pos2,$right, $filter)."\n";
 
}

close(FIN);
close(INTRA);
close(INTER1);
close(INTER2);


# identify the SV types
# ref: https://raw.githubusercontent.com/stat-lab/EvalSVcallers/master/scripts/vcf_convert/convert_GRIDSS_vcf.pl
sub get_SVtype_and_SVlen{

}