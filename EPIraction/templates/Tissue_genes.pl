#!/usr/bin/perl -w
use strict;

my $data_folder = "!{data_folder}";
my $tissues     = "!{tissues}";
my $blocks      = "!{blocks}";
my $chromosomes = "!{chromosomes}";

die "wrong folder $data_folder\n" unless -d $data_folder;
$tissues =~ s/\[|\]//g;
$chromosomes =~ s/\[|\]//g;

my $valid_genes = load_genes();

my @Tissues = sort split ", ", $tissues;
my $index = 1;

open FILE_OUT,">Tissue.genes.tsv";
print FILE_OUT "index\ttissue\tchromosome\tlist\n";
foreach my $tissue (@Tissues)
{
#    next unless $tissue eq "Blood.Myeloid.K562";
    open FILE_IN,"$data_folder/data/$tissue.regions.data";
    my $string = <FILE_IN>;
    my %Genes = ();

    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	if($Data[3]=~/^ENSG/)
	{
	    next unless defined $valid_genes->{$Data[3]};
	    push @{$Genes{$Data[0]}},$Data[3];
	}
    }
    close FILE_IN;
    foreach my $chrom (sort keys %Genes)
    {
	my @Report = ();
	foreach my $gene_id (@{$Genes{$chrom}})
	{
	    push @Report, $gene_id;
	    if($#Report >= $blocks-1)
	    {
		print FILE_OUT "$index\t$tissue\t$chrom\t",(join ";", @Report),"\n";
		$index++;
		@Report = ();
	    }
	}
	if($#Report >= 0)
	{
	    print FILE_OUT "$index\t$tissue\t$chrom\t",(join ";", @Report),"\n";
	    $index++;
	    @Report = ();
	}
    }
}

sub load_genes
{
    my %ret_h = ();
    foreach my $chrom (split ", ", $chromosomes)
    {
	opendir DIR,"$data_folder/Glmnet/$chrom";
	foreach my $file (grep { /rotations.matrix.gz/ } readdir DIR)
	{
	    $ret_h{(split /\./, $file)[0]}++;
	}
	closedir DIR;
    }
    return \%ret_h;
}
