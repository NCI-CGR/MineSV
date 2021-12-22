#!/usr/bin/env perl -w


@ARGV == 3 or die "$0 <INPUT VCF file> <Sample ID> <output Dir> \n";  

my $input=shift;
my $sample = shift;
my $output_dir = shift;
my $MIN_SV_LEN = 0;

##################################################
# To convert Meerkat to bed format 
# Wei Zhu
# 2021-09-23
##################################################

# Three files to be outputed here:
# ${outDir}${sample}_intra.bed: chr, start, end and line ID; start should be less than end. if start==end, end=start+1.
# ${outDir}${sample}_end1.bed
# ${outDir}${sample}_end2.bed 

### Sample input
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR
# chr6    27892614        bnd_1   .       .[chr6:27908860[        100     PASS    SVTYPE=BND;MATEID=bnd_2;EVENT=del_13546;GENE=;  GT:AD:DP:SS:SSC:BQ      .:.:.:.:.:.     .:5:.:2:100:.
# chr6    27908860        bnd_2   .       ]chr6:27892614].        100     PASS    SVTYPE=BND;MATEID=bnd_1;EVENT=del_13546;GENE=;  GT:AD:DP:SS:SSC:BQ      .:.:.:.:.:.     .:5:.:2:100:.
# chr8    143804684       bnd_3   .       ]chr1:78115979].        100     PASS    SVTYPE=BND;MATEID=bnd_4;EVENT=del_inss_14259_0/12324_0;GENE=;   GT:AD:DP:SS:SSC:BQ      .:.:.:.:.:.     .:5:.:2:100:.
# chr1    78115979        bnd_4   .       .[chr8:143804684[       100     PASS    SVTYPE=BND;MATEID=bnd_3;EVENT=del_inss_14259_0/12324_0;GENE=;   GT:AD:DP:SS:SSC:BQ      .:.:.:.:.:.     .:5:.:2:100:.
# chr8    143804583       bnd_5   .       .[chr1:67020324[        100     PASS    SVTYPE=BND;MATEID=bnd_6;EVENT=del_inss_14259_0/12324_0;GENE=;   GT:AD:DP:SS:SSC:BQ      .:.:.:.:.:.     .:6:.:2:100:.
# chr1    67020324        bnd_6   .       ]chr8:143804583].       100     PASS    SVTYPE=BND;MATEID=bnd_5;EVENT=del_inss_14259_0/12324_0;GENE=;   GT:AD:DP:SS:SSC:BQ      .:.:.:.:.:.     .:6:.:2:100:.
 

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
    
    my $NR = $.; 
    
    # my %info_hash = map { /([^=]+)=?(\S*)/ } split m{;}, $info; 

    
    my ($svType, $svLen, $chr2, $pos2, $orient1) = &get_SVtype_and_SVlen($chr1, $pos1,  $ref, $alt);

    #next if $svType eq "INV" & $orient1 eq "-"; # keep one pair is sufficient

    if($chr1 eq $chr2){
        # intrachromosomal
        next if $svLen < $MIN_SV_LEN; # ignore short SVs

        if($pos1 <= $pos2){
            print INTRA join("\t", $chr1, $pos1, $pos2, $NR)."\n";
        }else{
            # just hide the identical pare
            # print INTRA join("\t", $chr1, $pos2, $pos1, $NR)."\n";
        }
    }else{
        # interchromosomal
        # if( $chr1 lt $chr2){ # just show one pair
            $pos1 = ($pos1==0)?0:$pos1-1; 
            $pos2 = ($pos2==0)?0:$pos2-1; 
            print INTER1 join("\t", $chr1, $pos1, $pos1+1, $NR)."\n";
            print INTER2 join("\t", $chr2, $pos2, $pos2+1, $NR)."\n";
        # }
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
