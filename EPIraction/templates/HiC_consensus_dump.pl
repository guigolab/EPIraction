#!/usr/bin/perl

use strict;
use warnings;

my $data_folder  = "!{data_folder}";
my $chrom        = "!{chrom}";
my $regions      = "!{regions}";
my $genome_bed   = "!{genome_file}";
my $distance     = 2000000;
my $sample       = "Consensus.tissues.Encode.intact";

die "wrong folder $data_folder" unless -d $data_folder;
print "Working on $sample $chrom\n";
my $HiC_file = "$data_folder/HiC.files/$sample.hic";
my $prefix   = "$sample.$chrom";

if(-f "$data_folder/temp/promoter/$prefix.promoters.done")
{
    print "Already $data_folder/temp/promoter/$prefix.promoters.done\n";
    exit;
}

######################################################
####  Here I fit promoters into 2000 nt intervals ####
######################################################
my (%Promoters,%Good_starts) = ();
open PROMOTERS,">Promoters.bed";
open FILE_IN,"$data_folder/files/Human.promoters.bed";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next unless $Data[0] eq $chrom;
    my $middle = int(($Data[1] + $Data[2])/2);
    my ($start,$end) = get_positions($middle);
    print PROMOTERS "$chrom\t$start\t$end\t$Data[3]\t$middle\n";
    $Promoters{"$Data[3]"} = ["$chrom","$start","$end"];
    $start--;
    $Good_starts{$start}++;
}
close FILE_IN;
close PROMOTERS;
print "Done Promoters.bed\n";

#################################
####  Here I dump raw counts ####
#################################
my %Read_counts = ();
open FILE_IN,"java -Xmx5G -jar $data_folder/files/juicer_tools.jar dump observed NONE $HiC_file $chrom $chrom BP 2000 |";
my ($total,$good) = (0,0);
while(my $string = <FILE_IN>)
{
    unless($string=~/^\d+\t\d+/)
    {
#	print "Non informative string = $string";
	next;
    }
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[1] - $Data[0] > $distance;
    $total++;
    if(defined $Good_starts{$Data[0]} or defined $Good_starts{$Data[1]})
    {
	$Data[2] = sprintf("%d",$Data[2]);
	$Read_counts{"$Data[0]<>$Data[1]"} = $Data[2];
	$good++;
    }
}
close FILE_IN;
print "Done 2000 BP read counts\nGood starts = $good ",sprintf("%.6f\n",$good/$total);

################################################################
####  Here I dump HiC file 2000 bp SCALE, make bed and sort ####
################################################################
open FILE_IN,"java -Xmx5G -jar $data_folder/files/juicer_tools.jar dump observed SCALE $HiC_file $chrom $chrom BP 2000 | ";
open FILE_OUT,"| pigz -p2 > $prefix.not_sorted.bed.gz";
while(my $string = <FILE_IN>)
{
    unless($string=~/^\d+\t\d+/)
    {
#	print "Non informative string = $string";
	next;
    }
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[1] - $Data[0] > $distance;

    next if not defined $Good_starts{$Data[0]} and not defined $Good_starts{$Data[1]};

    $Data[2] = sprintf("%.6f",$Data[2]);
    my $real_reads = 0;
    $real_reads = $Read_counts{"$Data[0]<>$Data[1]"} if defined $Read_counts{"$Data[0]<>$Data[1]"};

    print FILE_OUT "$chrom\t",($Data[0]+1),"\t",($Data[0]+1999),"\t$Data[1]\t",($Data[1]-$Data[0]),"\t$Data[2]\t$real_reads\n";
    if($Data[0] != $Data[1])
    {
	print FILE_OUT "$chrom\t",($Data[1]+1),"\t",($Data[1]+1999),"\t$Data[0]\t",($Data[0]-$Data[1]),"\t$Data[2]\t$real_reads\n";
    }
}
close FILE_OUT;
close FILE_IN;
system "unpigz -p3 -c $prefix.not_sorted.bed.gz | LC_ALL=C sort --parallel=2 --buffer-size=5G -k 2,2n -k 4,4n - --temporary-directory='./' | pigz -p2 > $prefix.sorted.bed.gz";
print "Done 2000 BP sorted bed\n";

######################################################
####  Here I overlap this BED file with promoters ####
######################################################
open FILE_IN, "unpigz -p3 -c $prefix.sorted.bed.gz | bedtools intersect -a Promoters.bed -b stdin -wao -sorted -g $genome_bed |";
open FILE_OUT,"| LC_ALL=C sort --parallel=3 --buffer-size=5G -k 7,7n -k 4,4 -k 2,2n - --temporary-directory='./' | pigz -p3 > $prefix.promoters.gz";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[5] eq ".";
    print FILE_OUT "$Data[0]\t",($Data[8]+1),"\t",($Data[8]+1999),"\t$Data[3]\t",($Data[8]-$Data[6]+1),"\t$Data[10]\t$Data[4]\t$Data[11]\n";
}
close FILE_IN;
close FILE_OUT;
print "Done overlap promoters\n";
system "mv -f $prefix.promoters.gz $data_folder/temp/promoter/$prefix.promoters.gz";
system "touch $data_folder/temp/promoter/$prefix.promoters.done";
system "rm -f *.gz Promoters.bed";
print "Done\n";

sub get_positions
{
    my $middle = $_[0];
    my $start  = (int($middle/2000))*2000;
    my $end    = $start+2000;
    return ($start+1,$end-1) if $middle >= $start and $middle <= $end;
}
