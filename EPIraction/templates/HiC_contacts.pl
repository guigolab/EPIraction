#!/usr/bin/perl

use strict;
use warnings;

my $data_folder  = "!{data_folder}";
my $chromosomes  = "!{chromosomes}";
my $regions      = "!{regions}";
my $genome_bed   = "!{genome_file}";
my $baseline_hic = "!{HiC_baseline}";
my $distance     = 2000000;
my $sample       = "!{sample}";
my $baseline_min = 7;

die "not folder $data_folder" unless -d $data_folder;

my $HiC_baseline = load_contacts($baseline_hic);
$chromosomes =~ s/\[|\]//g;

foreach my $chrom (sort split ", ", $chromosomes)
{
    my $prefix = "$sample.$chrom";
    print "Do $sample $chrom\n";
    open LOG,">Work.log";

########################################################
####  Here I fit all regions into 2000 nt intervals ####
########################################################
    open ENHANCERS,">Enhancers.bed";
    open FILE_IN,"$data_folder/files/EPIraction.regions.bed";
    while(my $string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	next unless $Data[0] eq $chrom;
	my $middle = int(($Data[1] + $Data[2])/2);
	my ($start,$end) = get_positions($middle);
	print ENHANCERS "$chrom\t$start\t$end\t$Data[3]\t$middle\n";
    }
    close FILE_IN;
    close ENHANCERS;
    print "Done enhancer bed file\n";

##################################################
####  Here I calculate contact probabilities. ####
##################################################
    open FILE_OUT,"| pigz -p3 > $prefix.probabilities.gz";
    open FILE_IN,"unpigz -p3 -c $data_folder/temp/promoter/$prefix.promoters.gz |";
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
    print "Done probabilities\n";

#########################################
####  Here I overlap with enhancers  ####
#########################################
    open FILE_IN,"unpigz -p3 -c $prefix.probabilities.gz | bedSort stdin stdout | bedtools intersect -a Enhancers.bed -b stdin -wao -sorted -g $genome_bed |";
    open FILE_OUT,"| LC_ALL=C sort --parallel=6 --buffer-size=24G -k 6,6n -k 3,3n - --temporary-directory='./' | pigz -p3 > $prefix.overlap.gz";
    while(my $string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	next if $Data[5] eq ".";
	print FILE_OUT "$Data[8]\t$Data[3]\t$Data[9]\t$Data[10]\t$Data[11]\t$Data[12]\t$Data[13]\n";
    }
    close FILE_IN;
    close FILE_OUT;
    system "touch $prefix.overlap.done";
    system "mv -f $prefix.overlap.gz $data_folder/temp/contacts/$sample.$chrom.overlap.gz";
    system "mv -f $prefix.probabilities.gz $data_folder/temp/promoter/$sample.$chrom.probabilities.gz";
    system "rm -f *.gz Enhancers.bed *.done";
    system "mv Work.log $data_folder/temp/contacts/$sample.$chrom.overlap.log";
    print "Done Enhancers\n";
    print "Done $chrom\n";
}

sub report_probability
{
    my $list = $_[0];
    my $gene_id = $list->[0]->[3];

    my %Distances = ();
    foreach my $line (@{$list})
    {
	my $dist       = $line->[4];
	my $contacts   = $line->[5];
	my $promoter   = $line->[6];
	my $real_reads = $line->[7];
	my $prefix     = "$line->[0]\t$line->[1]\t$line->[2]\t$line->[3]";
	$Distances{$dist} = [$prefix,$contacts,0,$promoter,$real_reads];
    }
    return if not defined $Distances{"-2000"};
    return if not defined $Distances{"2000"};
    return if not defined $Distances{"-4000"};
    return if not defined $Distances{"4000"};

    my $baseline = ($Distances{"-4000"}->[1] + $Distances{"4000"}->[1])/2;
    if($baseline < $baseline_min)
    {
	$baseline = (sort { $b <=> $a } ($Distances{"-4000"}->[1],$Distances{"4000"}->[1]))[0];
    }
    return if $baseline < $baseline_min;

    my ($corrected,$total)= (0,0);
    foreach my $distance (keys %Distances)
    {
	my $probability = $Distances{$distance}->[1]/$baseline;
	$probability    = 1    if $distance == 0;
	$probability    = 0.99 if abs($distance) == 2000;
	$probability    = 0.95 if abs($distance) == 4000;
	$probability    = 1    if $probability > 1;
	$probability    = sprintf("%.8f",$probability);
	my $real_reads  = $Distances{$distance}->[4];

	my $baseline_probability = $HiC_baseline->{abs($distance)};
	$Distances{$distance}->[4] = "$real_reads<>passed";
	$Distances{$distance}->[2] = $probability;
	if($real_reads == 1)
	{
#########------ I cannot trust if we have one reads supporting contact, thus I limit the probability by 1.0 baseline ------######
	    if($probability > $baseline_probability)
	    {
		$Distances{$distance}->[4] = "$real_reads<>corrected";
		$Distances{$distance}->[2] = sprintf("%.8f",$baseline_probability);
		$corrected++;
	    }
	}
	elsif($real_reads == 2)
	{
#########------ I still cannot trust if we have two reads, thus I limit the probability by 1.5 baseline ------######
	    if($probability > $baseline_probability*1.5)
	    {
		$Distances{$distance}->[4] = "$real_reads<>corrected";
		$Distances{$distance}->[2] = sprintf("%.8f",$baseline_probability*1.5);
		$corrected++;
	    }
	}
	else
	{
	    $Distances{$distance}->[4] = "$real_reads<>substantial";
	}
	$total++;
    }
    print LOG sprintf("%s --  baseline = %10.4f -- corrected = %3d   %.4f\n",$gene_id,$baseline,$corrected,$corrected/$total);
    map { print FILE_OUT $Distances{$_}->[0],"\t",$_,"\t",$Distances{$_}->[1],"\t",$Distances{$_}->[2],"\t",$Distances{$_}->[3],"\t",$Distances{$_}->[4],"\n" } sort { $a <=> $b } keys %Distances;
    print FILE_OUT "###\n";
}

sub get_positions
{
    my $middle = $_[0];
    my $start  = (int($middle/2000))*2000;
    my $end    = $start+2000;
    return ($start+1,$end-1) if $middle >= $start and $middle <= $end;
}

sub load_contacts
{
    my $file  = $_[0];
    my %ret_h = ();
    open FILE_IN,$file or die "Absent $file";
    my $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$ret_h{$Data[0]} = $Data[1];
    }
    close FILE_IN;
    return \%ret_h;
}
