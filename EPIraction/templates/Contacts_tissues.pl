#!/usr/bin/perl

use strict;
use warnings;

my $data_folder  = "!{data_folder}";
my $tissue       = "!{tissue}";
my $regions      = "!{regions}";
my $HiC_samples  = "!{HiC_samples}";
my $chromosomes  = "!{chromosomes}";
my $baseline_hic = "!{HiC_baseline}";
my $upper_hic    = "!{HiC_upper}";     ### not in use eventually

die "wrong folder $data_folder\n" unless -d $data_folder;

my $HiC_baseline = load_contacts($baseline_hic);
my $HiC_upper    = load_contacts($upper_hic);
$chromosomes     =~ s/\[|\]//g;

my %Regions = ();
open FILE_IN,"unpigz -p3 -c $data_folder/data/$tissue.regions.data.gz |";
my $string = <FILE_IN>;
chomp $string;
my @Header = split "\t", $string;

while($string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my ($start,$end,$activity) = @Data[1,2,11];
    $Regions{"$Data[3]"} = ["$start","$end","$activity"];
}
close FILE_IN;

open FILE_OUT,"| pigz -p3 > Tissue.contacts.data.gz";
print FILE_OUT "chrom\tpromoter\tenhancer\tdistance\tactivity\tcontact\trescaled\tHiC-fold\n";
open FILE_LOG,">Log.txt";
print FILE_LOG "ensembl_gene_id\tfrom_baseline\tfrom_consensus\tfrom_tissue\tnonbaseline_fraction\n";
foreach my $chrom (sort split ", ", $chromosomes)
{
    my %Contacts = ();
    foreach my $HiC_sample (split ";", $HiC_samples)
    {
	open FILE_IN,"zcat $data_folder/temp/contacts/$HiC_sample.$chrom.overlap.gz |";
	while($string = <FILE_IN>)
	{
	    chomp $string;
	    next if $string eq "";
	    my ($promoter,$enhancer,$distance,$reads,$contact,undef,undef) = split "\t", $string;
	    $enhancer =~ s/^\s+//;
	    unless(defined $Contacts{"$promoter<>$enhancer"})
	    {
		$Contacts{"$promoter<>$enhancer"} = $contact;
	    }
	    else
	    {
		$Contacts{"$promoter<>$enhancer"} = $contact if $Contacts{"$promoter<>$enhancer"} < $contact;
	    }
	}
	close FILE_IN;
    }
    print "Load tissues HiC for $chrom\n";

    my %Consensus = ();
    open FILE_IN,"unpigz -p3 -c $data_folder/temp/contacts/Consensus.tissues.Encode.intact.$chrom.overlap.gz |";
    while($string = <FILE_IN>)
    {
	chomp $string;
	next if $string eq "";
	my ($promoter,$enhancer,$distance,$reads,$contact,undef,undef) = split "\t", $string;
	$enhancer =~ s/^\s+//;
	$Consensus{"$promoter<>$enhancer"} = $contact;
    }
    close FILE_IN;
    print "Load consensus HiC for $chrom\n";

    open FILE_IN,"unpigz -p3 -c $data_folder/data/zzz.Pairs.$chrom.gz |";
    $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my ($promoter,undef,$list) = split "\t", $string;
	next unless defined $Regions{$promoter};
	my ($p_start,$p_end,$p_activity) = @{$Regions{$promoter}};
	my ($better_tissue,$better_consensus,$better_baseline) = (0,0,0);
	foreach my $enhancer (split ";", $list)
	{
	    next unless defined $Regions{$enhancer};
	    my ($e_start,$e_end,$e_activity) = @{$Regions{$enhancer}};
	    next if $e_activity == 0 and $enhancer !~/ENSG/;

	    my @Sorted       = sort { $a <=> $b } ($p_start,$p_end,$e_start,$e_end);
	    my $distance     = $Sorted[2] - $Sorted[1];
	    $distance        = 0 if $enhancer eq $promoter;
	    my $distance_200 = get_200_distance($distance);
	    my $baseline     = $HiC_baseline->{$distance_200};
############--- set deafault contacts to baseline ---############
	    my $contact      = $baseline;

	    $Consensus{"$promoter<>$enhancer"} = $baseline unless defined $Consensus{"$promoter<>$enhancer"};
	    $Contacts{"$promoter<>$enhancer"}  = $baseline unless defined $Contacts{"$promoter<>$enhancer"};
	    my $consensus_contacts = $Consensus{"$promoter<>$enhancer"};
	    my $tissue_contacts    = $Contacts{"$promoter<>$enhancer"};

	    if($baseline >= $consensus_contacts and $baseline >= $tissue_contacts)
	    {
		$better_baseline++;
	    }
	    elsif($consensus_contacts > $baseline and $consensus_contacts >= $tissue_contacts)
	    {
############--- consensus is higher than baseline and tissue ---############
		$contact = $consensus_contacts;
		$better_consensus++;
	    }
	    else
	    {
############--- tissue is higher than baseline and consensus ---############
		$contact = $tissue_contacts;
		$better_tissue++;
	    }
	    $contact = 1 if $contact > 1;

	    my $hic_fold = sprintf("%.3f",$contact/$baseline);
	    $contact     = sprintf("%.8f",$contact);
	    my $rescaled = sprintf("%.9f",$e_activity*$contact);
	    $e_activity  = sprintf("%.8f",$e_activity);
	    $enhancer    = sprintf("%27s",$enhancer);
	    print FILE_OUT "$chrom\t$promoter\t$enhancer\t$distance\t$e_activity\t$contact\t$rescaled\t$hic_fold\n";
	}
	my $fraction = sprintf("%.6f",($better_consensus+$better_tissue)/($better_consensus + $better_tissue + $better_baseline));
	print FILE_LOG "$promoter\t$better_baseline\t$better_consensus\t$better_tissue\t$fraction\n";
	print FILE_OUT "###\n";
    }
    close FILE_IN;
    print "Done $tissue.$chrom\n";
}
close FILE_OUT;
close FILE_LOG;
system "mv -f Tissue.contacts.data.gz $data_folder/data/$tissue.contacts.data.gz";
system "mv -f Log.txt $data_folder/data/$tissue.contacts.log";
print "Done all\n";

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

sub get_200_distance
{
    my $value = $_[0];
    return 0 if $value == 0;
    return 2000000 if $value >= 1999900;

    my $times = int(($value+100)/200);
    return $times*200;
}
