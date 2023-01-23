#!/usr/bin/env perl
use strict;

##################################################
#### If EPIraction.config variables are valid ####
##################################################
open FILE_IN,"EPIraction.config" or die "no EPIraction.config file present in current folder\nExiting\n";
my @Lines = <FILE_IN>;
close FILE_IN;
chomp @Lines;

my @Good = grep { /data_folder\s*=\s*\'(.+)\'/;$_ = $1 } map { "$_" } @Lines;
die "No data folder provided in EPIraction.config\nExpecting \"data_folder = '/users/project/encode_005982/EPIraction.Nextflow'\" in param block\nExiting\n" if $#Good == -1;
die "Found many instances of data_folder:\n".(join "\n", @Good)."\nExiting\n" if $#Good > 0;

my $data_folder = $Good[0];
die "$data_folder is not a folder\nExiting\n" unless -d $data_folder;

#### Good ####
print "1.1 $data_folder                                    Good so far\n";
##############

@Good = grep { /tissues_index\s*=\s*\'(.+)\'/;$_ = $1 } map { "$_" } @Lines;
die "\nNo tissues_index file provided in EPIraction.config\nExpecting \"tissues_index = 'EPIraction.tissues.data'\" in param block\nExiting\n" if $#Good == -1;
die "Found many instances of tissues_index:\n".(join "\n", @Good)."\nExiting\n" if $#Good > 0;

my $tissues_index = $Good[0];
die "$data_folder/files/$tissues_index is not a file\nExiting\n" unless -f "$data_folder/files/$tissues_index";

#### Good ####
print "1.2 $data_folder/files/$tissues_index      Good so far\n";
##############

@Good = grep { /samples_index\s*=\s*\'(.+)\'/;$_ = $1 } map { "$_" } @Lines;
die "\nNo samples_index file provided in EPIraction.config\nExpecting \"samples_index = 'EPIraction.H3K27ac.data'\" in param block\nExiting\n" if $#Good == -1;
die "Found many instances of samples_index:\n".(join "\n", @Good)."\nExiting\n" if $#Good > 0;

my $samples_index = $Good[0];
die "$data_folder/files/$samples_index is not a file\nExiting\n" unless -f "$data_folder/files/$samples_index";

#### Good ####
print "1.3 $data_folder/files/$samples_index      Good so far\n";
##############

@Good = grep { /regions\s*=\s*\'(.+)\'/;$_ = $1 } map { "$_" } @Lines;
die "\nNo regions file provided in EPIraction.config\nExpecting \"regions       = 'EPIraction.regions.bed'\" in param block\nExiting\n" if $#Good == -1;
die "Found many instances of regions:\n".(join "\n", @Good)."\nExiting\n" if $#Good > 0;

my $regions = $Good[0];
die "$data_folder/files/$regions is not a file\nExiting\n" unless -f "$data_folder/files/$regions";

#### Good ####
print "1.4 $data_folder/files/$regions       Good so far\n\n";
##############

##############################################
#### If all nessesary files are available ####
##############################################
die "\nThere is no Human.genome file in $data_folder/files\nExiting\n" unless -f "$data_folder/files/Human.genome";
print "2.1 $data_folder/files/Human.genome                 Good so far\n";

die "\nThere is no Gencode.v40.promoters.data file in $data_folder/files\nExiting\n" unless -f "$data_folder/files/Gencode.v40.promoters.data";
print "2.2 $data_folder/files/Gencode.v40.promoters.data   Good so far\n";

die "\nThere is no HiC.baseline file in $data_folder/files\nExiting\n" unless -f "$data_folder/files/HiC.baseline";
print "2.3 $data_folder/files/HiC.baseline                 Good so far\n";

die "\nThere is no juicer_tools.jar file in $data_folder/files\nExiting\n" unless -f "$data_folder/files/juicer_tools.jar";
print "2.4 $data_folder/files/juicer_tools.jar             Good so far\n";

die "\nThere is no Downstream.penalty file in $data_folder/files\nExiting\n" unless -f "$data_folder/files/Downstream.penalty";
print "2.5 $data_folder/files/Downstream.penalty           Good so far\n\n";

die "\nThere is no Upstream.penalty file in $data_folder/files\nExiting\n" unless -f "$data_folder/files/Upstream.penalty";
print "2.6 $data_folder/files/Upstream.penalty             Good so far\n\n";


####################################################
#### If all BigWigs and Hic files are available ####
####################################################


my ($tissues,undef) = load_data("$data_folder/files/$tissues_index");
my ($samples,undef) = load_data("$data_folder/files/$samples_index");
foreach my $tissue_id (sort { $a <=> $b } keys %{$tissues})
{
    my $data   = $tissues->{$tissue_id};
    my $tissue = $data->{"tissue"};

    my $sample = $data->{"Open.Chromatin"};
    die "\n$data_folder/BigWig/$sample.bigWig (Open.Chromatin) absent\nExiting\n" unless -f "$data_folder/BigWig/$sample.bigWig";

    $sample = $data->{"Cofactor"};
    die "\n$data_folder/BigWig/$sample.bigWig (Cofactor) absent\nExiting\n" unless -f "$data_folder/BigWig/$sample.bigWig";

    die "\n$data_folder/BigWig/$tissue.H3K27ac.consensus.bigWig absent\nExiting\n" unless -f "$data_folder/BigWig/$tissue.H3K27ac.consensus.bigWig";

    foreach my $HiC (split ";", $data->{"HiC"})
    {
	die "\n$data_folder/HiC.files/$HiC.hic absent\nExiting\n" unless -f "$data_folder/HiC.files/$HiC.hic";
    }

    foreach my $sample_id (grep { $samples->{$_}->{"tissue"} eq $tissue } sort { $a <=> $b } keys %{$samples})
    {
	my $sample_data = $samples->{$sample_id};
	my $sample = $sample_data->{"file"};
	die "\n$data_folder/BigWig/$sample.H3K27ac.bigWig absent\nExiting\n" unless -f "$data_folder/BigWig/$sample.H3K27ac.bigWig";
    }
    my $report = "3.$tissue_id";
    $report.=" " while length $report < 6;
    $report.="tissue=\"$tissue\", data files";
    $report.=" " while length $report < 87;
    print  "$report Good so far\n";
}
print "\nAll files are there\n";

sub load_data
{
    open FILE_L,$_[0] or die "$! $_[0]";
    my $string = <FILE_L>;chomp $string;$string=~s/\r//;
    my @Head = split "\t",$string;
    my $num_rows=$#Head;
    my %Hash_data=();
    while($string = <FILE_L>)
    {
        chomp $string;$string=~s/\r//;
        $string.="NA" if $string=~/\t$/;
        my @Content = split "\t",$string;
        %{$Hash_data{$Content[0]}}=();
        my $hash_ref = $Hash_data{$Content[0]};
        map { $hash_ref->{$Head[$_]}= $Content[$_]} (1..$num_rows);
    }
    close FILE_L;
    shift @Head;
    return (\%Hash_data,\@Head);
}
