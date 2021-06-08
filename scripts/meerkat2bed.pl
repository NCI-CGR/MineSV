#!/usr/bin/perl
# meerkat2vcf.pl
# convert meerkat output to vcf format
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
# Email: lixing_yang@hms.harvard.edu

###############################
# Modified by Wei Zhu below: 
# I have no need to use ref and vcf header 
# Input: @ARGV == 3 or die "$0 <INPUT VCF file> <Sample ID> <output Dir> \n";
# The output:
# # Three files to be outputed here:
# ${outDir}${sample}_intra.bed: chr, start, end and line ID; start should be less than end. if start==end, end=start+1.
# ${outDir}${sample}_end1.bed
# ${outDir}${sample}_end2.bed 
################################
use strict;


@ARGV == 3 or die "$0 <INPUT Meerkat file> <Sample ID> <output Dir> \n";
my $input=shift;
my $sample = shift;
my $output_dir = shift;

# open file handles
my $basefn = $output_dir.'/'.$sample;

$input = "zcat $input |" if ($input =~ 'gz$'); 

open(FIN, $input) or die $!; 
open(INTRA, '>', $basefn."_intra.bed") or die $!; 
open(INTER1, '>', $basefn."_end1.bed") or die $!;
open(INTER2, '>', $basefn."_end2.bed") or die $!;


my $newline;
my $mateid = 1;
my $NR;

while ($newline = <FIN>)
{
	$NR = $.;
	chomp $newline;
	my (undef, @data) = split (/\t/, $newline);
	my $event = $data[0].'_'.$data[2];
	my @ids;
	if ($data[3] =~ /\//)
	{
		@ids = split (/\//, $data[3]);
	}
	if ($data[0] eq 'del')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 0;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7], $data[11], $data[3], $event);
	}
	elsif ($data[0] eq 'del_ins')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 0;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7], $data[14], $data[3], $event);
	}
	elsif ($data[0] eq 'invers_f')
	{
		$event = 'transl_intra_'.$data[2];
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 1;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7], $data[11], $data[3], $event);
	}
	elsif ($data[0] eq 'invers_r')
	{
		$event = 'transl_intra_'.$data[2];
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 2;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7], $data[11], $data[3], $event);
	}
	elsif ($data[0] eq 'tandem_dup')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 3;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7], $data[11], $data[3], $event);
	}
	elsif ($data[0] eq 'invers')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 1;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6]-1, $data[5], $data[7], $data[11], $ids[0], $event);
		
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 2;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7]+1, $data[11], $ids[1], $event);
	}
	elsif ($data[0] =~ /inssd/)
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 3;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[11], $data[16], $ids[0], $event);
		
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 0;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[10], $data[17], $ids[1], $event);
	}
	elsif ($data[0] =~ /inssu/)
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 3;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[10], $data[16], $ids[0], $event);
		
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 0;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[11], $data[17], $ids[1], $event);
	}
	elsif ($data[0] =~ /inso/)
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 1;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[11], $data[16], $ids[0], $event);
		
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 2;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[10], $data[17], $ids[1], $event);
	}
	elsif ($data[0] eq 'del_invers')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 1;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[11], $data[17], $ids[0], $event);
		
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 2;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[10], $data[18], $ids[1], $event);
	}
	elsif ($data[0] =~ /inss/)
	{
		if ($data[5] lt $data[9])
		{
			my $mateida = 'bnd_'.$mateid;
			$mateid++;
			my $mateidb = 'bnd_'.$mateid;
			$mateid++;
			my $mate_type = 3;
			&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[10], $data[16], $ids[0], $event);
			
			my $mateida = 'bnd_'.$mateid;
			$mateid++;
			my $mateidb = 'bnd_'.$mateid;
			$mateid++;
			my $mate_type = 0;
			&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[11], $data[17], $ids[1], $event);
		}
		else
		{
			my $mateida = 'bnd_'.$mateid;
			$mateid++;
			my $mateidb = 'bnd_'.$mateid;
			$mateid++;
			my $mate_type = 3;
			&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[11], $data[17], $ids[0], $event);
			
			my $mateida = 'bnd_'.$mateid;
			$mateid++;
			my $mateidb = 'bnd_'.$mateid;
			$mateid++;
			my $mate_type = 0;
			&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[10], $data[16], $ids[1], $event);
		}
	}
	elsif ($data[0] eq 'transl_inter')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 0;
		$mate_type = 1 if ($data[7] == 1 and $data[10] == 1);
		$mate_type = 2 if ($data[7] == -1 and $data[10] == -1);
		$mate_type = 3 if ($data[7] == -1 and $data[10] == 1);
		# print join("\t", $data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[8], $data[9], $data[13], $data[3], $event);
		# die;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[8], $data[9], $data[13], $data[3], $event);
	}
}
close(FIN);
close(INTRA);
close(INTER1);
close(INTER2);

sub bnd{
	my $event_type = shift;
	my $mate_type = shift;
	my $mateida = shift;
	my $mateidb = shift;
	my $chra = shift;
	my $posa = shift;
	my $chrb = shift;
	my $posb = shift;
	my $gene_string = shift;
	my $ad = shift;
	my $event = shift;
	
	my ($genea, $geneb);
	if ($gene_string =~ /GN:(\S{0,20}),.*;(\S{0,20}),/)
	{
		($genea, $geneb) = ($1, $2);
		$genea = substr ($genea, 0, -1);
		$geneb = substr ($geneb, 0, -1);
	}
	#print "$gene_string\t$genea\t$geneb\n";
	
	my ($alta, $altb);
	my $nta = '.';
	my $ntb = '.';
	# $nta = $reference_db->seq($chra, $posa => $posa) if (defined($reference_path));
	# $ntb = $reference_db->seq($chrb, $posb => $posb) if (defined($reference_path));
	# if ($mate_type == 0)
	# {
	# 	$altb = $nta.'['.$chrb.':'.$posb.'[';
	# 	$alta = ']'.$chra.':'.$posa.']'.$ntb;
	# }
	# if ($mate_type == 1)
	# {
	# 	$altb = $nta.']'.$chrb.':'.$posb.']';
	# 	$alta = $ntb.']'.$chra.':'.$posa.']';
	# }
	# if ($mate_type == 2)
	# {
	# 	$altb = '['.$chrb.':'.$posb.'['.$nta;
	# 	$alta = '['.$chra.':'.$posa.'['.$ntb;
	# }
	# if ($mate_type == 3)
	# {
	# 	$altb = ']'.$chrb.':'.$posb.']'.$nta;
	# 	$alta = $ntb.'['.$chra.':'.$posa.'[';
	# }
	
	if($chra eq $chrb){
        # intrachromosomal
        if($posa <= $posb){
            print INTRA join("\t", $chra, $posa, $posb, $NR)."\n";
        }else{
            print INTRA join("\t", $chra, $posb, $posa, $NR)."\n";
        }
    }else{
        # interchromosomal
        $posa = ($posa==0)?0:$posa-1; 
        $posb = ($posb==0)?0:$posb-1; 
        print INTER1 join("\t", $chra, $posa, $posa+1, $NR)."\n";
        print INTER2 join("\t", $chrb, $posb, $posb+1, $NR)."\n";
    }
	return;
}

