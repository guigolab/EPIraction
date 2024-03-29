#!/usr/bin/perl

use strict;
use warnings;

my $data_folder  = "!{data_folder}";
my $sample       = "!{sample}";
my $chrom        = "!{chrom}";
my $regions      = "!{regions}";
my $genome_bed   = "!{genome_file}";
my $blanks_pl    = "!{blanks_pl}";
my $distance     = 1000000;

die "not folder $data_folder" unless -d $data_folder;
print "Working on $sample $chrom\n";

my $HiC_file = "$data_folder/HiC.files/$sample.hic";
my $prefix   = "$sample.$chrom";

if(-f "$data_folder/temp/HiC/$sample.$chrom.overlap.gz")
{
    print "Already $data_folder/temp/HiC/$sample.$chrom.overlap.gz\n\n";
    exit;
}

####################################################################
####  Here I fit promoters and enhancers into 2000 nt intervals ####
####################################################################
open FILE_IN,"$data_folder/files/$regions";
open PROMOTERS,">Promoters.bed";
open ENHANCERS,">Enhancers.bed";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next unless $Data[0] eq $chrom;
    my $middle = int(($Data[1] + $Data[2])/2);
    my ($start,$end) = get_positions($middle);
    if($Data[3] =~ /^ENSG/)
    {
	print PROMOTERS "$chrom\t$start\t$end\t$Data[3]\t$middle\n";
    }
    else
    {
	print ENHANCERS "$chrom\t$start\t$end\t$Data[3]\t$middle\n";
    }
}
close FILE_IN;
close ENHANCERS;
close PROMOTERS;
print "Done promoter and enhancer bed files\n";

###############################################################################
####  Here I dump HiC file BP 1000 and decrease it to make 2000 resolution ####
####  I ignore the digonal and the bins, located 2000 and 3000 nt away.    ####
###############################################################################
unless(-f "$prefix.2000.unsorted.done")
{
    my %Pairs = ();
    open FILE_IN,"java -Xmx15G -jar $data_folder/files/juicer_tools.jar dump observed SCALE $HiC_file $chrom $chrom BP 1000 | ";
    while(my $string = <FILE_IN>)
    {
	unless($string=~/^\d+\t\d+/)
	{
	    print "Non informative string = $string";
	    next;
	}
	chomp $string;
	my @Data = split "\t", $string;
	$Data[0] -= 1000 if $Data[0]%2000 == 1000;
	$Data[1] -= 1000 if $Data[1]%2000 == 1000;
	next if $Data[1] - $Data[0] > $distance;
	next if $Data[1] - $Data[0] <= 2000;
	$Pairs{"$Data[0]\t$Data[1]"} += $Data[2];
    }
    close FILE_IN;

    open FILE_OUT, "| pigz -p3 > $prefix.2000.unsorted.gz";
    map { print FILE_OUT $_,"\t",sprintf("%.6f",$Pairs{$_}),"\n"; } keys %Pairs;
    close FILE_OUT;
    undef(%Pairs);
    print "Done 2000 BP unsorted dump\n";
    system "touch $prefix.2000.unsorted.done";
}
###############################################################
####  Here I wrote sorted BED file from this unsorted file ####
###############################################################
unless(-f "$prefix.sorted.bed.done")
{
    open FILE_OUT,"| LC_ALL=C sort --parallel=3 --buffer-size=16G -k 2,2n -k 4,4n - --temporary-directory='./' | pigz -p3 > $prefix.sorted.bed.gz";
    open FILE_IN,"unpigz -p3 -c $prefix.2000.unsorted.gz |";
    my %Diagonal = ();
    while(my $string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$Data[2] = sprintf("%.6f",$Data[2]);
	print FILE_OUT "$chrom\t",($Data[0]+1),"\t",($Data[0]+1999),"\t$Data[1]\t",($Data[1]-$Data[0]),"\t$Data[2]\n";
	print FILE_OUT "$chrom\t",($Data[1]+1),"\t",($Data[1]+1999),"\t$Data[0]\t",($Data[0]-$Data[1]),"\t$Data[2]\n";
	$Diagonal{$Data[0]}++;
	$Diagonal{$Data[1]}++;
    }
    close FILE_IN;
    foreach my $position (keys %Diagonal)
    {
	print FILE_OUT "$chrom\t",($position+1),"\t",($position+1999),"\t",($position),"\t0\t0\n";
    }
    close FILE_OUT;
    undef(%Diagonal);
    print "Done 2000 BP sorted bed\n";
    system "touch $prefix.sorted.bed.done";
}
#########################################################
####  Here I overlap this BED file it with promoters ####
#########################################################
unless(-f "$prefix.promoters.done")
{
    open FILE_IN, "unpigz -p3 -c $prefix.sorted.bed.gz | bedtools intersect -a Promoters.bed -b stdin -wao -sorted -g $genome_bed |";
    open FILE_OUT,"| LC_ALL=C sort --parallel=3 --buffer-size=16G -k 4,4 -k 2,2n - --temporary-directory='./' | pigz -p3 > $prefix.promoters.gz";
    while(my $string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	next if $Data[5] eq ".";
	print FILE_OUT "$Data[0]\t",($Data[8]+1),"\t",($Data[8]+1999),"\t$Data[3]\t",($Data[8]-$Data[6]+1),"\t$Data[10]\n";
    }
    close FILE_IN;
    close FILE_OUT;
    print "Done Promoters\n";
    system "touch $prefix.promoters.done";
}

######################################################################################################
####  Here I calculate contact probabilities. I load Consensus.tissues.Encode.intact HiC data and ####
####  form the zero and +/- 2000 diagonal values promotionally to Consensus ones for every gene   ####
######################################################################################################
open FILE_IN,"unpigz -p3 -c $data_folder/temp/promoter/Consensus.tissues.Encode.intact.$chrom.probabilities.gz |";
my %Consensus = ();
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    push @{$Consensus{$Data[3]}},\@Data;
}
close FILE_IN;
print "Load consensus\n";

open FILE_OUT,"| pigz -p3 > $prefix.probabilities.gz";
open FILE_IN,"unpigz -p3 -c $prefix.promoters.gz |";
my @Lines = ();
my $string = <FILE_IN>;
chomp $string;
my @Data = split "\t", $string;
push @Lines,\@Data;

while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    if($Data[3] eq $Lines[-1]->[3])
    {
	push @Lines,\@Data;
    }
    else
    {
	&report_probability(\@Lines);
	@Lines = ();
	push @Lines,\@Data;
    }
}
&report_probability(\@Lines);
@Lines = ();
close FILE_IN;
close FILE_OUT;
system "touch $prefix.probabilities.done";
print "Done probabilities\n";

#########################################
####  Here I overlap with enhancers  ####
#########################################
open FILE_IN,"unpigz -p3 -c $prefix.probabilities.gz | bedSort stdin stdout | bedtools intersect -a Enhancers.bed -b stdin -wao -sorted -g $genome_bed |";
open FILE_OUT,"| LC_ALL=C sort --parallel=3 --buffer-size=16G -k 1,1 -k 3,3n - --temporary-directory='./' | $blanks_pl | pigz -p3 > $prefix.overlap.gz";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[5] eq ".";
    print FILE_OUT "$Data[8]\t$Data[3]\t$Data[9]\t$Data[10]\t$Data[11]\n";
}
close FILE_IN;
close FILE_OUT;
system "touch $prefix.overlap.done";
print "Done Enhancers\n";

system "mv -f $prefix.overlap.gz $data_folder/temp/HiC/$sample.$chrom.overlap.gz";
system "mv -f $prefix.probabilities.gz $data_folder/temp/promoter/$sample.$chrom.probabilities.gz";
system "rm -f *.gz Enhancers.bed Promoters.bed *.done";
print "Done\n";

sub get_positions
{
    my $middle = $_[0];
    my $start  = (int($middle/2000))*2000;
    my $end    = $start+2000;
    return ($start+1,$end-1) if $middle >= $start and $middle <= $end;
}

sub report_probability
{
    my $list = $_[0];
    my $gene_id = $list->[0]->[3];
    return unless defined $Consensus{$gene_id};

    my @Sorted = sort { $b <=> $a } map { "$_->[5]" } grep { abs($_->[4]) <= 25000 } @{$list};
    return if $Sorted[0] < 5;
    return if $#Sorted < 5;
    my $dist_2000 = 2*$Sorted[0]-$Sorted[1];
    my $dist_0000 = 3*$Sorted[0]-$Sorted[1]-$Sorted[2];

    my (%Current,%Complete,%Altogether)  = ();
    map { $Current{$_->[1]}  = ["$_->[0]","$_->[1]","$_->[2]","$_->[3]","$_->[4]","$_->[5]"] } @{$list};
    map { $Complete{$_->[1]} = ["$_->[0]","$_->[1]","$_->[2]","$_->[3]","$_->[4]","$_->[5]","$_->[6]"] } @{$Consensus{$gene_id}};
    map { $Altogether{$_}++ } keys %Current;
    map { $Altogether{$_}++ } keys %Complete;

    my @Combined = ();
    foreach my $position (sort { $a <=> $b } keys %Altogether)
    {
	if(defined $Current{$position} and defined $Complete{$position})
	{
	    my @Report = (@{$Complete{$position}},$Current{$position}->[5]);
	    push @Combined,\@Report;
	}
	elsif(not defined $Current{$position} and defined $Complete{$position})
	{
	    my @Report = (@{$Complete{$position}},0);
	    push @Combined,\@Report;
	}
	elsif(defined $Current{$position} and not defined $Complete{$position})
	{
	    my @Report = ($Current{$position}->[0],$Current{$position}->[1],$Current{$position}->[2],$Current{$position}->[3],$Current{$position}->[4],0,0,$Current{$position}->[5]);
	    push @Combined,\@Report;
	}
    }

    my ($non_diagonal_consensus,$non_diagonal_current) = (0,0);
    foreach my $line (@Combined)
    {
	next if abs($line->[4]) <= 2000;
	$non_diagonal_consensus += $line->[5];
	$non_diagonal_current   += $line->[7];
    }
    return if $non_diagonal_current < 10;
    my $rescale = $non_diagonal_current/$non_diagonal_consensus;
    print "$gene_id:\n";
    print "non_diagonal_consensus\t$non_diagonal_consensus\n";
    print "non_diagonal_current\t$non_diagonal_current\n";
    print "rescale = ",sprintf("%.6f",$rescale),"\n";

    foreach my $index (0..$#Combined)
    {
	next if abs($Combined[$index]->[4]) > 2000;
	if(abs($Combined[$index]->[4]) == 2000)
	{
	    $dist_2000 = $Combined[$index]->[5]*$rescale if $dist_2000 < $Combined[$index]->[5]*$rescale;
	    $Combined[$index]->[7] = sprintf("%.6f",$dist_2000);
	}
	if($Combined[$index]->[4]==0)
	{
	    $dist_0000 = $Combined[$index]->[5]*$rescale if $dist_0000 < $Combined[$index]->[5]*$rescale;
	    $Combined[$index]->[7] = sprintf("%.6f",$dist_0000);
	}
    }

    foreach my $line (@Combined)
    {
	next if $line->[7] == 0;
	my $contact = $line->[7]/$dist_0000;
	$contact = 1 if $contact > 1;
	print FILE_OUT join "\t",@{$line}[0..4];
	print FILE_OUT sprintf("\t%.6f\t%.8f\n",$line->[7],$contact);
    }
}
