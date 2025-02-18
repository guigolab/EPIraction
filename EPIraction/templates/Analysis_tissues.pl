#!/usr/bin/perl -w
use strict;

my $data_folder = "!{data_folder}";
my $tissues     = "!{tissues}";
my $blocks      = "!{blocks}";
my $chromosomes = "!{chromosomes}";
my $minimum_exp = "!{minimum_exp}";

die "wrong folder $data_folder\n" unless -d $data_folder;

$tissues =~ s/\[|\]//g;
$chromosomes =~ s/\[|\]//g;

my $valid_genes = load_genes();

my @Tissues = sort split ", ", $tissues;

open FILE_OUT,">Tissue.genes.tsv";
print FILE_OUT "index\ttissue\tchromosome\tlist\n";
my ($jobs,$total)=(0,0);
foreach my $tissue (sort @Tissues)
{
    my %Genes = ();
    open FILE_IN,"$data_folder/data/$tissue.expression.data";
    my $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my ($gene_id,$chrom,$type,$expression,$symbol) = split "\t", $string;
	next if $expression < $minimum_exp;
	next unless defined $valid_genes->{$gene_id};
	push @{$Genes{$chrom}},$gene_id;
    }
    close FILE_IN;
    my $index = 1;
    foreach my $chrom (sort keys %Genes)
    {
	my @Report = ();
	foreach my $gene_id (@{$Genes{$chrom}})
	{
	    push @Report, $gene_id;
	    if($#Report >= $blocks-1)
	    {
		unless(-f "$data_folder/temp/batches/$tissue.$index.pairs.data.gz")
		{
		    print FILE_OUT "$index\t$tissue\t$chrom\t",(join ";", @Report),"\n";
		    $jobs++;
		}
		$total++;
		$index++;
		@Report = ();
	    }
	}
	if($#Report >= 0)
	{
	    unless(-f "$data_folder/temp/batches/$tissue.$index.pairs.data.gz")
	    {
		print FILE_OUT "$index\t$tissue\t$chrom\t",(join ";", @Report),"\n";
		$jobs++;
	    }
	    $index++;
	    $total++;
	    @Report = ();
	}
    }
}
close FILE_OUT;
sleep(5);
print "Done ",($total-1)," batches\n$jobs jobs to run\n";

sub load_genes
{
    my %ret_h = ();
    foreach my $chrom (split ", ", $chromosomes)
    {
	opendir DIR,"$data_folder/Genes/$chrom";
	foreach my $file (grep { /matrix.gz/ } readdir DIR)
	{
	    $ret_h{(split /\./, $file)[0]}++;
	}
	closedir DIR;
    }
    return \%ret_h;
}
