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
#@Tissues = ("Blood.Leukemia.Jurkat","Blood.Lymphoid.B_cell","Blood.Lymphoid.CD4_conventional","Blood.Lymphoid.CD4_memory","Blood.Lymphoid.CD4_regulatory","Blood.Lymphoid.CD8_conventional","Blood.Lymphoid.GM12878","Blood.Lymphoid.NK_cell","Blood.Myeloid.K562","Blood.Myeloid.Monocyte","Brain.Cortex.Primary","Brain.Neurons","Breast.Adenocarcinoma.MCF7","Breast.MCF10A","Digestive.Colon.Primary.Mix","Digestive.Stomach.Primary","Kidney.Primary","Liver.Carcinoma.HepG2","Liver.Primary","Lung.Fibroblast","Lung.Primary","Mesenchyme.Adipose","Mesenchyme.Endothelial.HAEC","Muscle.Heart","Muscle.Skeletal","Muscle.Smooth.Artery","Pancreas.Islets","Pancreas.Primary.Body","Pluripotent.hESCs.H9","Pluripotent.hESCs.Rest","Pluripotent.iPSCs","Prostate.Primary","Skin.Fibroblast","Skin.Keratinocyte","Spleen.Primary");

open FILE_OUT,">Tissue.genes.tsv";
print FILE_OUT "index\ttissue\tchromosome\tlist\n";
my ($jobs,$total)=(0,0);
foreach my $tissue (sort @Tissues)
{
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
	    next if $Data[10] < 0.5;
	    push @{$Genes{$Data[0]}},$Data[3];
	}
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
		unless(-f "$data_folder/temp/Glmnet/$tissue.$index.glmnet.log")
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
	    unless(-f "$data_folder/temp/Glmnet/$tissue.$index.glmnet.log")
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
	opendir DIR,"$data_folder/Glmnet/$chrom";
	foreach my $file (grep { /matrix.gz/ } readdir DIR)
	{
	    $ret_h{(split /\./, $file)[0]}++;
	}
	closedir DIR;
    }
    return \%ret_h;
}
