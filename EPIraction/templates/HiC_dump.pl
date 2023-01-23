#!/usr/bin/perl

use strict;

my $data_folder = "!{data_folder}";
my $sample      = "!{sample}";
my $chrom       = "!{chrom}";
my $regions     = "!{regions}";
my $distance    = 1000000;

die "not folder $data_folder" unless -d $data_folder;
print "Working on $sample $chrom\n";

my $HiC_file = "$data_folder/HiC.files/$sample.hic";
my $prefix   = "Current";

if(-f "$data_folder/temp/HiC/$sample.$chrom.overlap.gz")
{
    print "$data_folder/temp/HiC/$sample.$chrom.overlap.gz";
    exit;
}

open FILE_IN,"java -jar $data_folder/files/juicer_tools.jar dump observed SCALE $HiC_file $chrom $chrom BP 5000 | ";

open FILE_OUT,"| LC_ALL=C sort --parallel=2 --buffer-size=12G -k 1,1n -k 2,2n - --temporary-directory='./' | gzip > $prefix.dump.gz";

while(my $string = <FILE_IN>)
{
    chomp $string;
    next unless $string=~/^\d+\t\d+/;
    my @Data = split "\t", $string;
    print FILE_OUT "$Data[0]\t$Data[1]\t$Data[2]\n";
    last;
}

while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[1] - $Data[0] > $distance;
    print FILE_OUT "$Data[0]\t$Data[1]\t$Data[2]\n";
}
close FILE_IN;
close FILE_OUT;

my %Diagonal = ();
open FILE_IN,"zcat $prefix.dump.gz |";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[1] == $Data[0];
    if(not defined $Diagonal{$Data[0]})
    {
	$Diagonal{$Data[0]} = $Data[2];
    }
    else
    {
	$Diagonal{$Data[0]} = $Data[2] if $Diagonal{$Data[0]} < $Data[2];
    }

    if(not defined $Diagonal{$Data[1]})
    {
	$Diagonal{$Data[1]} = $Data[2];
    }
    else
    {
	$Diagonal{$Data[1]} = $Data[2] if $Diagonal{$Data[1]} < $Data[2];
    }
}
close FILE_IN;

open FILE_IN,"zcat $prefix.dump.gz |";
open FILE_OUT,"| gzip > $prefix.diagonal.gz";

while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    if($Data[0] != $Data[1])
    {
	print FILE_OUT $string,"\n";
    }
    else
    {
	if(not defined $Diagonal{$Data[0]})
	{
	    my $scale = $Data[2]/4;
	    print FILE_OUT "$Data[0]\t$Data[0]\t$scale\t$Data[2]\n";
	}
	else
	{
	    my $scale = $Diagonal{$Data[0]};
	    print FILE_OUT "$Data[0]\t$Data[0]\t$scale\t$Data[2]\n";
	}
    }
}
close FILE_IN;
close FILE_OUT;

system "zcat $prefix.diagonal.gz | annotate_pairs.awk -v chrom=$chrom | LC_ALL=C sort --parallel=2 --buffer-size=12G -k 2,2n -k 4,4n - --temporary-directory='./' | gzip > $prefix.intersect.gz";

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
	print PROMOTERS "$chrom\t$start\t$end\t$Data[3]\t$Data[1]\n";
    }
    else
    {
	print ENHANCERS "$chrom\t$start\t$end\t$Data[3]\t$Data[1]\n";
    }
}
close FILE_IN;
close ENHANCERS;
close PROMOTERS;

open FILE_IN, "zcat $prefix.intersect.gz | bedtools intersect -a Promoters.bed -b stdin -wao -sorted -g $data_folder/files/Human.genome |";
my %Contact_data = ();
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[5] eq ".";
    push @{$Contact_data{$Data[3]}}, ["$Data[8]",($Data[8]-$Data[6]),"$Data[10]","$Data[11]"];
}
close FILE_IN;

open FILE_OUT,"| gzip > $prefix.promoters.gz";
open FILE_IN,"Promoters.bed";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my $promoter = $Data[3];
    my $p_start  = $Data[4];
    next unless $Contact_data{$promoter};
    my @Report_data = @{$Contact_data{$promoter}};

    my @Center = grep { $Report_data[$_]->[1] == 0 } (0..$#Report_data);
    next if $#Center == -1;
    my $center_index = $Center[0];
    my $norm = $Report_data[$center_index]->[2];

    foreach my $line (@Report_data)
    {
	my $observed = $line->[2];
	$observed    = $norm if $observed > $norm;
	my $contact  = sprintf("%.7f",$observed/$norm);
	print FILE_OUT "$chrom\t$line->[0]\t",($line->[0]+4999),"\t$promoter\t$line->[0]\t$line->[1]\t$contact\t$observed\t$p_start\n";
    }
}
close FILE_IN;
close FILE_OUT;

open FILE_IN,"zcat $prefix.promoters.gz | bedSort stdin stdout | bedtools intersect -a Enhancers.bed -b stdin -wao -sorted -g $data_folder/files/Human.genome |";
open FILE_OUT,"| LC_ALL=C sort --parallel=2 --buffer-size=15G -k 5,5n -k 6,6n - --temporary-directory='./' | blanks.pl | gzip > $prefix.overlap.gz";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[5] eq ".";
#    map { print $_,"\t",$Data[$_],"\n"; } (0..$#Data);exit;
    print FILE_OUT "$Data[8]\t$Data[3]\t$Data[10]\t$Data[11]\t$Data[13]\t$Data[4]\n";
}
close FILE_IN;
close FILE_OUT;
system "mv $prefix.overlap.gz $data_folder/temp/HiC/$sample.$chrom.overlap.gz";
system "rm -f Current.diagonal.gz Current.dump.gz Current.intersect.gz Current.promoters.gz Enhancers.bed Promoters.bed";

sub get_positions
{
    my $middle = $_[0];
    my $start  = int($middle/5000);
    my $before = $start*5000+1;
    my $next   = $before+5000-2;
    return ($before,$next);
}
