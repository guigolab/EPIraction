#!/usr/bin/env perl
use strict;

my $work_folder = `pwd`;
chomp $work_folder;

##########################################
#### If nessesary folders are present ####
##########################################
print "If nessesary folders are present:\n";
my $index = 1;
foreach my $folder ("BigWig","data","files","Glmnet","HiC.files","report","RSEM","temp","temp/Glmnet","temp/HiC","temp/promoter","temp/signal")
{
    if(-d "$work_folder/$folder")
    {
	print "1.$index Present $work_folder/$folder\n";
    }
    else
    {
	print "1.$index Absent $work_folder/$folder\n";
	exit;
    }
    $index++;
}
print "All folders are here\n\n";

##########################################################
#### If nessesary index and support files are present ####
##########################################################
print "If nessesary index and support files are present:\n";
$index = 1;
foreach my $file ("EPIraction.H3K27ac.data","EPIraction.regions.bed","EPIraction.tissues.data","Gencode.v40.promoters.data","Human.promoters.bed","juicer_tools.jar")
{
    if(-f "$work_folder/files/$file")
    {
	print "2.$index Present $work_folder/files/$file\n";
    }
    else
    {
	print "2.$index Absent $work_folder/files/$file\n";
	exit;
    }
    $index++;
}
print "All such files are here\n\n";

#####################################################
#### If bigWigs and expression files are present ####
#####################################################
my ($tissues,undef) = load_data("$work_folder/files/EPIraction.tissues.data");
my ($samples,undef) = load_data("$work_folder/files/EPIraction.H3K27ac.data");
my @Samples = sort { $a <=> $b } keys %{$samples};

print "If all sample-specific H3K27ac files are present:\n";
foreach my $tissue_index (sort { $a <=> $b } keys %{$tissues})
{
    my $tissue = $tissues->{$tissue_index}->{"tissue"};
    my $total = 0;
    foreach my $sample_id ( grep { $samples->{$_}->{"tissue"} eq $tissue } @Samples)
    {
	my $file = $samples->{$sample_id}->{"file"};
	unless(-f "$work_folder/BigWig/$file.H3K27ac.bigWig")
	{
	    print "\n$work_folder/BigWig/$file.H3K27ac.bigWig absent\n";
	    exit;
	}
	$total++;
    }
    print "3.$tissue_index Found $total *.H3K27ac.bigWig files for $tissue\n";
}
print "All sample-specific H3K27ac files are here\n\n";

print "If all tissue-consensus files are present:\n";
foreach my $tissue_index (sort { $a <=> $b } keys %{$tissues})
{
    my $tissue = $tissues->{$tissue_index}->{"tissue"};

#####  RNAseq   #####
    my $label  = $tissues->{$tissue_index}->{"RNA-seq"};
    unless(-f "$work_folder/RSEM/$label.genes.results")
    {
	print "\n$work_folder/RSEM/$label.genes.results absent\n";
	exit;
    }

#####  Open   #####
    $label  = $tissues->{$tissue_index}->{"Open.Chromatin"};
    unless(-f "$work_folder/BigWig/$label.bigWig")
    {
	print "\n$work_folder/BigWig/$label.bigWig absent\n";
	exit;
    }

#####  Cofactor   #####
    $label  = $tissues->{$tissue_index}->{"Cofactor"};
    unless(-f "$work_folder/BigWig/$label.bigWig")
    {
	print "\n$work_folder/BigWig/$label.bigWig absent\n";
	exit;
    }

#####  Hi-C   #####
    foreach my $label (split ";", $tissues->{$tissue_index}->{"HiC"})
    {
	unless(-f "$work_folder/HiC.files/$label.hic")
	{
	    print "\n$work_folder/HiC.files/$label.hic absent\n";
	    exit;
	}
    }
    print "4.$tissue_index Found all tissue-consensus files for $tissue\n";
}
unless(-f "$work_folder/HiC.files/Consensus.tissues.Encode.intact.hic")
{
    print "\n$work_folder/HiC.files/Consensus.tissues.Encode.intact.hic absent\n";
    exit;
}
print "4.78 Consensus.tissues.Encode.intact.hic present\n";
print "All tissue-consensus files are here\n\n";

print "It looks like all files are present.\nDo not forget to put the current folder '$work_folder'\ninto 'EPIraction.config' file within params section\nDone QC!!!!\n";

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
