#!/usr/bin/perl -w

use strict;
use warnings;

my $data_folder   = "!{data_folder}";
my $regions       = "!{regions}";
my $tissue        = "!{tissue}";
my $RNA_seq       = "!{RNA_seq}";
my $Open          = "!{Open}";
my $Cofactor      = "!{Cofactor}";

die "wrong folder $data_folder\n" unless -d $data_folder;

my $expression    = get_expression($RNA_seq);
my ($genes,undef) = load_data("$data_folder/files/Gencode.v40.promoters.data");

my %Open_Chrom = ();
my $bigWig = "$data_folder/BigWig/$Open.bigWig";
open FILE_IN, "bigWigAverageOverBed $bigWig $data_folder/files/$regions -minMax stdout |";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my $max  = sprintf("%.4f",$Data[7]);
    $Open_Chrom{$Data[0]} = "$max";
}
close FILE_IN;
print "Load $Open\n";

my %Cofactor = ();
$bigWig     = "$data_folder/BigWig/$Cofactor.bigWig";
open FILE_IN, "bigWigAverageOverBed $bigWig $data_folder/files/$regions -minMax stdout |";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my $max  = sprintf("%.4f",$Data[7]);
    $Cofactor{$Data[0]} = "$max";
}
close FILE_IN;
print "Load $Cofactor\n";

my %Consensus = ();
$bigWig = "$data_folder/BigWig/$tissue.H3K27ac.consensus.bigWig";
open FILE_IN, "bigWigAverageOverBed $bigWig $data_folder/files/$regions -minMax stdout |";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my $max  = sprintf("%.4f",$Data[7]);
    my $mean = sprintf("%.4f",$Data[4]);
    $Consensus{$Data[0]} = "$max";
}
close FILE_IN;
print "Load $tissue.H3K27ac.consensus\n";

my ($promoters,$enhancers) = (0,0);
open FILE_OUT,">tissue.regions.data";
print FILE_OUT "#chrom\tstart\tend\tname\tscore\tstrand\ttype\tH3K27ac\tOpen\tCofactor\tTPMs\tSymbol\tSample\n";
open FILE_IN,"$data_folder/files/$regions";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    if($Data[3] =~ /^ENSG/)
    {
	my $gene_id = $Data[3];
	next unless defined $expression->{$gene_id};
	next if $expression->{$gene_id} < 0.5;
	next unless $genes->{$gene_id};
	
	$Open_Chrom{$gene_id} = sprintf("%.3f",$Open_Chrom{$gene_id}) if $Open_Chrom{$gene_id} > 10;
	$Open_Chrom{$gene_id} = sprintf("%.2f",$Open_Chrom{$gene_id}) if $Open_Chrom{$gene_id} > 100;
	$Consensus{$gene_id}  = sprintf("%.3f",$Consensus{$gene_id})  if $Consensus{$gene_id}  > 10;
	$Consensus{$gene_id}  = sprintf("%.2f",$Consensus{$gene_id})  if $Consensus{$gene_id}  > 100;
	$Cofactor{$gene_id}   = sprintf("%.3f",$Cofactor{$gene_id})   if $Cofactor{$gene_id}   > 10;
	$Cofactor{$gene_id}   = sprintf("%.2f",$Cofactor{$gene_id})   if $Cofactor{$gene_id}   > 100;

	print FILE_OUT "$string\t1000\t",$genes->{$gene_id}->{"strand"},"\tpromoter\t",$Consensus{$gene_id},"\t",$Open_Chrom{$gene_id},"\t",$Cofactor{$gene_id},"\t";
	print FILE_OUT $expression->{$gene_id},"\t",$genes->{$gene_id}->{"Symbol"},"\t$tissue\n";
	$promoters++;
    }
    else
    {
	my $enhancer = $Data[3];
	next if $Open_Chrom{$enhancer} < 1.5;
	next if $Consensus{$enhancer}  < 1.5;
	next if $Cofactor{$enhancer}   < 1.5;

	$Open_Chrom{$enhancer} = sprintf("%.3f",$Open_Chrom{$enhancer}) if $Open_Chrom{$enhancer} > 10;
	$Open_Chrom{$enhancer} = sprintf("%.2f",$Open_Chrom{$enhancer}) if $Open_Chrom{$enhancer} > 100;
	$Consensus{$enhancer}  = sprintf("%.3f",$Consensus{$enhancer})  if $Consensus{$enhancer}  > 10;
	$Consensus{$enhancer}  = sprintf("%.2f",$Consensus{$enhancer})  if $Consensus{$enhancer}  > 100;
	$Cofactor{$enhancer}   = sprintf("%.3f",$Cofactor{$enhancer})   if $Cofactor{$enhancer}   > 10;
	$Cofactor{$enhancer}   = sprintf("%.2f",$Cofactor{$enhancer})   if $Cofactor{$enhancer}   > 100;

	print FILE_OUT "$string\t1000\t.\tenhancer\t",$Consensus{$enhancer},"\t",$Open_Chrom{$enhancer},"\t",$Cofactor{$enhancer},"\t0.0000\tenhancer\t$tissue\n";
	$enhancers++;
    }
}
close FILE_IN;
close FILE_OUT;
print "Signal $tissue\npromoters = $promoters\nenhancers = $enhancers\n";
system "mv tissue.regions.data $data_folder/temp/tissues/$tissue.regions.data";

sub get_nearest
{
    my $value = $_[0];
    my $times = int(($value+25)/50);
    return $times*50;
}

sub get_expression
{
    my $accession = $_[0];
    my %ret_h = ();
    open FILE_IN,"$data_folder/RSEM/$accession.genes.results";
    my $string = <FILE_IN>;
    chomp $string;
    my @Header = split "\t", $string;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$Data[0] =~ s/\.\d\d?//;
	$ret_h{$Data[0]} = $Data[5];
#	map { print $_,"\t",$Header[$_],"\t",$Data[$_],"\n" } (0..$#Data);exit;
    }
    close FILE_IN;
    return \%ret_h;
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
