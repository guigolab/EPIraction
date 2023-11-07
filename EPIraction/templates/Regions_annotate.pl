#!/usr/bin/perl

use strict;
use warnings;

my $data_folder   = "!{data_folder}";
my $regions       = "!{regions}";
my $tissue        = "!{tissue}";
my $RNA_seq       = "!{RNA_seq}";
my $Open          = "!{Open}";
my $Cofactor      = "!{Cofactor}";
my $transform_aaa = "!{transform_aaa}";
my $transform_bbb = "!{transform_bbb}";
my $log_two       = 1/log(2);

die "wrong folder $data_folder\n" unless -d $data_folder;

my $expression    = get_expression($RNA_seq);
my ($genes,undef) = load_data("$data_folder/files/Gencode.v40.promoters.data");

open FILE_LOG,">Annotate_regions.log";

######################
### Open Chromatin ###
######################
my %Open_Chrom = ();
my ($mean,$count,$sd,$max) = (0,0,0,0);
my $bigWig = "$data_folder/BigWig/$Open.bigWig";
open FILE_IN, "bigWigAverageOverBed $bigWig $data_folder/files/$regions -minMax stdout |";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data   = split "\t", $string;
    my $signal = $Data[7];
    $mean +=$signal;
    $count++;
    $max = $signal if $signal > $max;

    $signal = sprintf("%.4f",$signal);
    $Open_Chrom{$Data[0]} = "$signal";
}
close FILE_IN;
$mean = sprintf("%.6f",$mean/$count);
map { $sd += ($_ - $mean)**2 } values %Open_Chrom;
$sd   = sprintf("%.6f",sqrt($sd/$count));
$max  = sprintf("%.6f",$max);

print "Load $Open\n$mean\t$sd\t$max\n";
print FILE_LOG "Open chromatin = $Open\n$mean\t$sd\t$max\n\n";

######################
###    Cofactor    ###
######################
my %Cofactor = ();
($mean,$count,$sd,$max) = (0,0,0,0);
$bigWig     = "$data_folder/BigWig/$Cofactor.bigWig";
open FILE_IN, "bigWigAverageOverBed $bigWig $data_folder/files/$regions -minMax stdout |";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my $signal = $Data[7];
    $mean +=$signal;
    $count++;
    $max = $signal if $signal > $max;

    $signal = sprintf("%.4f",$signal);
    $Cofactor{$Data[0]} = "$signal";
}
close FILE_IN;
$mean = sprintf("%.6f",$mean/$count);
map { $sd += ($_ - $mean)**2 } values %Cofactor;
$sd   = sprintf("%.6f",sqrt($sd/$count));
$max  = sprintf("%.6f",$max);

print "Load $Cofactor\n$mean\t$sd\t$max\n";
print FILE_LOG "Cofactor = $Cofactor\n$mean\t$sd\t$max\n\n";

######################
###    H3K27ac     ###
######################
my %Consensus = ();
($mean,$count,$sd,$max) = (0,0,0,0);
$bigWig = "$data_folder/BigWig/$tissue.H3K27ac.consensus.bigWig";
open FILE_IN, "bigWigAverageOverBed $bigWig $data_folder/files/$regions -minMax stdout |";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my $signal = $Data[7];
    $signal = 90 if $signal > 90;
    $mean +=$signal;
    $count++;
    $max = $signal if $signal > $max;

    $signal = sprintf("%.4f",$signal);
    $Consensus{$Data[0]} = "$signal";
}
close FILE_IN;
$mean = sprintf("%.6f",$mean/$count);
map { $sd += ($_ - $mean)**2 } values %Cofactor;
$sd   = sprintf("%.6f",sqrt($sd/$count));
$max  = sprintf("%.6f",$max);

print "Load $tissue.H3K27ac.consensus\n$mean\t$sd\t$max\n";
print FILE_LOG "Consensus = $tissue.H3K27ac.consensus\n$mean\t$sd\t$max\n\n";

############################
### Altogether. Activity ###
############################
my ($promoters,$enhancers) = (0,0);
($mean,$count,$sd,$max) = (0,0,0,0);
my @List = ();

open FILE_OUT,">tissue.regions.data";
print FILE_OUT "#chrom\tstart\tend\tname\tscore\tstrand\ttype\tH3K27ac\tOpen\tCofactor\tTPMs\tSymbol\tSample\tActivity\n";

open FILE_IN,"$data_folder/files/$regions";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    if($Data[3] =~ /^ENSG/)
    {
	my $gene_id = $Data[3];
	next unless defined $expression->{$gene_id};
	next unless $genes->{$gene_id};
	my $strand   = $genes->{$gene_id}->{"strand"};
	my $H3K27ac  = $Consensus{$gene_id};
	my $Open     = $Open_Chrom{$gene_id};
	my $Cofactor = $Cofactor{$gene_id};

	my $gene_exp = $expression->{$gene_id};
	my $response = transform_expression($gene_exp);
	my $symbol   = $genes->{$gene_id}->{"Symbol"};

	my $Activity = sprintf("%.6f",(12*$Open*$H3K27ac)**0.5);
	$Activity    = sprintf("%.6f",$response*0.25)  if $Activity < $response*0.25;
	$Activity    = sprintf("%.6f",100)  if $Activity > 100;
	print FILE_OUT "$string\t1000\t$strand\tpromoter\t$H3K27ac\t$Open\t$Cofactor\t$gene_exp\t$symbol\t$tissue\t$Activity\n";

	$mean += $Activity;
	$count++;
	$max = $Activity if $max < $Activity;
	push @List,$Activity;
	$promoters++;
    }
    else
    {
	my $enhancer = $Data[3];
	my $strand   = ".";
	my $H3K27ac  = $Consensus{$enhancer};
	my $Open     = $Open_Chrom{$enhancer};

	my $Cofactor = $Cofactor{$enhancer};
	my $gene_exp = "0";
	my $symbol   = "enhancer";

#####################################################################################################
### Minimum 2.5 value of H3K27ac and 1.5 of Open chromatin requered. Pseudocount of 4 to Cofactor ###
#####################################################################################################
	if($H3K27ac < 2.5 or $Open < 1.5)
	{
	    print FILE_OUT "$string\t1000\t$strand\tenhancer\t$H3K27ac\t$Open\t$Cofactor\t$gene_exp\t$symbol\t$tissue\t0\n";
	    $count++;
	    next;
	}
	my $Activity = sprintf("%.6f",($Open*$H3K27ac*($Cofactor+4))**0.5);
	$Activity    = sprintf("%.6f",100)  if $Activity > 100;
	print FILE_OUT "$string\t1000\t$strand\tenhancer\t$H3K27ac\t$Open\t$Cofactor\t$gene_exp\t$symbol\t$tissue\t$Activity\n";

	$mean += $Activity;
	$count++;
	$max = $Activity if $max < $Activity;
	push @List,$Activity;
	$enhancers++;
    }
}
close FILE_IN;
close FILE_OUT;
$mean = sprintf("%.6f",$mean/$count);
map { $sd += ($_ - $mean)**2 } values @List;
$sd   = sprintf("%.6f",sqrt($sd/$count));
$max  = sprintf("%.6f",$max);

print "Signal $tissue\npromoters = $promoters\nenhancers = $enhancers\nActivity  = $mean\t$sd\t$max\n";
print FILE_LOG "Activity\n$mean\t$sd\t$max\n";
close FILE_LOG;

system "pigz -p2 -f tissue.regions.data";
system "mv -f tissue.regions.data.gz $data_folder/temp/tissues/$tissue.regions.data.gz";
system "mv -f Annotate_regions.log   $data_folder/temp/tissues/$tissue.signals.log";

sub transform_expression
{
    my $TPM_value = $_[0];
    my $log_expression = $log_two*log(0.01 + $TPM_value)*$transform_aaa + $transform_bbb;
    return $log_expression;
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
