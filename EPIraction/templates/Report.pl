#!/usr/bin/perl -w
use strict;

my $data_folder = "!{data_folder}";
my $tissue      = "!{tissue}";

die "wrong folder $data_folder\n" unless -d $data_folder;

my $HiC_baseline  = load_baseline();
my ($genes,undef) = load_data("$data_folder/files/Gencode.v40.promoters.data");

open PROMOTERS,">Genes.bed";
print PROMOTERS "#chrom\tstart\tend\tgene_id\tscore\tstrand\tSymbol\tTSS\tTPMs\tH3K27ac\tOpen\tSample\n";

open ENHANCERS,">Enhancers.bed";
print ENHANCERS "#chrom\tstart\tend\tenhancer\tscore\tstrand\tH3K27ac\tOpen\tCofactor\tSample\n";

my (%Promoters,%Enhancers) = ();
open FILE_IN,"$data_folder/data/$tissue.regions.data";
my $string = <FILE_IN>;
while($string = <FILE_IN>)
{
    my @Data = split "\t", $string;
    if($Data[3]=~/^ENSG/)
    {
	my $gene_data = $genes->{$Data[3]};
	print PROMOTERS join "\t", @Data[0..5];
	print PROMOTERS "\t",$gene_data->{"Symbol"},"\t",$gene_data->{"TSS"},"\t",(join "\t", @Data[10,7,8,12]),"\n";
    }
    else
    {
	my $signals  = join "\t",@Data[7..9];
	my $position = join "\t",@Data[0..3];
	print ENHANCERS join "\t", @Data[0..5];
	print ENHANCERS "\t$signals\t$tissue\n";
	$Enhancers{$Data[3]} = [$position,$signals];
    }
}
close FILE_IN;
close PROMOTERS;
close ENHANCERS;

open FILE_IN,"zcat $data_folder/data/$tissue.scalings.data.gz |";
my %HiC = ();
$string = <FILE_IN>;
while($string = <FILE_IN>)
{
    chomp $string;
    next if $string eq "";
    my @Data = split "\t", $string;
    my $HiC_distance = get_HiC_distance($Data[4]);
    $HiC{"$Data[2]<>$Data[3]"} = [$Data[4],sprintf("%.5f\t%.5f\t%.5f",$Data[6],$Data[6]/$HiC_baseline->{$HiC_distance},$Data[8])];
}
close FILE_IN;

open FILE_OUT,">Pairs.bed";
print FILE_OUT "#chr\tstart\tend\tname\tclass\tTargetGene\tTargetGeneEnsemblID\tTargetGeneTSS\tCellType\tScore\tDistanceToTSS\tH3K27ac\tOpen\tCofactor\tHiC_contacts\tHiC_foldchange\tscalings\n";

opendir DIR,"$data_folder/temp/Glmnet";
foreach my $file (sort { compare_files($a,$b) } grep { /$tissue/ }  grep{ /pairs.data.gz/} readdir DIR)
{
    open FILE_IN,"zcat $data_folder/temp/Glmnet/$file |";
    $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	if($string =~/^#/)
	{
	    print FILE_OUT "#\n";
	    next;
	}
	my ($promoter,$enhancer,$impact) = split "\t", $string;
	next if $enhancer =~ /promoter/;
	next unless defined $Enhancers{$enhancer};
	my $gene_data = $genes->{$promoter};
	print FILE_OUT $Enhancers{$enhancer}->[0],"\tenhancer\t",$gene_data->{"Symbol"},"\t$promoter\t",$gene_data->{"TSS"},"\t$tissue\t$impact\t";
	print FILE_OUT $HiC{"$promoter<>$enhancer"}->[0],"\t",$Enhancers{$enhancer}->[1],"\t",$HiC{"$promoter<>$enhancer"}->[1],"\n";
    }
    close FILE_IN;
    print "Load $file\n";
}
close FILE_OUT;
system "mv Genes.bed $data_folder/report/EPIraction.$tissue.genes.bed";
system "mv Enhancers.bed $data_folder/report/EPIraction.$tissue.enhancers.bed";
system "mv Pairs.bed $data_folder/report/EPIraction.$tissue.pairs.bed";

sub compare_files
{
    my ($AAA,$BBB) = @_;
    my @AAA = split /\./, $AAA;
    my @BBB = split /\./, $BBB;
    return $AAA[-4] <=> $BBB[-4];
}

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

sub load_baseline
{
    my %ret_h = ();
    open FILE_IN,"$data_folder/files/HiC.baseline";
    $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$ret_h{$Data[0]} = $Data[1];
    }
    close FILE_IN;
    return \%ret_h;
}
sub get_HiC_distance
{
    my $value = $_[0];
    return 1000000 if $value >= 997500;
    return 0       if $value <=   2500;
    my $times = int(($value+2500)/5000);
    return $times*5000;
}
