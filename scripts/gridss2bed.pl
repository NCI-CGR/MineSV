#!/usr/bin/env perl -w


@ARGV == 3 or die "$0 <INPUT VCF file> <Sample ID> <output Dir> \n";  

my $input=shift;
my $sample = shift;
my $output_dir = shift;
my $MIN_SV_LEN = 50;

##################################################
# To convert GRIDSS V4.2 VCF file to bed format 
# Wei Zhu
# 2020-10-27 13:50:09
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

    
    my ($svType, $svLen, $chr2, $pos2, $orient1) = &get_SVtype_and_SVlen($chr1, $pos1,  $ref, $alt);

    next if $svType eq "INV" & $orient1 eq "-"; # keep one pair is sufficient

    if($chr1 eq $chr2){
        # intrachromosomal
        next if $svLen < $MIN_SV_LEN; # ignore short SVs

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


    # print join("\t", $svType, $chr1, $pos1, $id, $ref, $alt, $info_hash{'EVENT'}, $NR, $left, $orient, $chr2, $pos2,$right, $filter)."\n";
 
}

close(FIN);
close(INTRA);
close(INTER1);
close(INTER2);


# identify the SV types
# ref: https://raw.githubusercontent.com/stat-lab/EvalSVcallers/master/scripts/vcf_convert/convert_GRIDSS_vcf.pl (hard to follow)
# https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9

# 1       25204182        25210853        UID_34  DUP     6671
# 1       32176186        32193833        UID_45  INV     17647
# 1       36419037        36421453        UID_53  DEL     2416
# 1       37852934        37863685        UID_59  INV     10751
# 1       40302119        40310285        UID_62  DEL     8166
# 1       52748187        52752297        UID_73  DEL     4110
sub get_SVtype_and_SVlen{
    my ($chr1, $pos1, $ref, $alt)=@_;


    # insertion:  1       66572   gridss0fb_85o   G       GTACTATATATTA[1:66573[  1121.28 PASS

    my ($left, $bracket, $chr2, $pos2,$right) = ($alt =~ /^(.*)([\[\]])(.+):(.+)\2(.*)$/);

    my $orient1="+";
    my $insLen=0;
    my $alt_clean=undef;
    my $svLen= undef;
    my $svType = "ITX";

    if ($left ne ""){
        $orient1="+";
        $alt_clean=$left;
    }else{
        $orient1="-";
        $alt_clean=$right
    }
    
    if(length($alt_clean)> length($ref)){
        $insLen = length($alt_clean)- length($ref);
    }

    my $orient2=($bracket eq '[')?"+":"-";



    if($chr1 eq $chr2){
        $svLen=abs($pos1-$pos2); # TODO, sv lenglth needs to be refined by SV types
        
        if($orient1 ne $orient2){
            $svType="INV";
        }else{
            
            if($insLen > $svLen*0.7){
                $svType = "INS";
                $svLen=$insLen; # use insertion size for the sv length
            }else{
                if( ($pos1 < $pos2) xor ($orient1 eq '-')){
                    $svType="DEL";
                }else{
                    $svType="DUP";
                }

            }

        }
    }

    # print "\n".join("\t", $svType, $alt, $left, $right, $orient1, $bracket, $orient2, $pos1, $pos2)."\n";

    return($svType, $svLen, $chr2, $pos2, $orient1);
}