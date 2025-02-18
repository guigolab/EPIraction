#!/usr/bin/perl

use strict;
use warnings;

my $data_folder = "!{data_folder}";
my $regions     = "!{regions}";
my $tissue      = "!{tissue}";
my $RNA_seq     = "!{RNA_seq}";
my $Open_bw     = "!{Open}";
my $Cofactor_bw = "!{Cofactor}";
my $CTCF_bw     = "!{CTCF}";
my $genome_bed  = "!{genome_file}";
my $min_exp     = "!{expression}";
my $H3K27ac_min = "!{H3K27ac_min}";
my $Open_min    = "!{Open_min}";
my $CTCF_min    = "!{CTCF_min}";
my $ctcf_impact = 0.5;
srand(127);

die "wrong folder $data_folder\n" unless -d $data_folder;
my $expression  = get_expression($RNA_seq);

my ($genes,undef) = load_data("$data_folder/files/Gencode.v40.genes.data");

open FILE_LOG,">Annotate_regions.log";
print FILE_LOG "H3K27ac_min = $H3K27ac_min\nOpen_min    = $Open_min\nCTCF_min    = $CTCF_min\n\n";

######################
### Open Chromatin ###
######################
my %Open_Chrom = ();
my ($mean,$count,$sd,$max,$peaks) = (0,0,0,0,0);
my $bigWig = "$data_folder/BigWig/$Open_bw.bigWig";
open FILE_IN, "bigWigAverageOverBed $bigWig $data_folder/files/$regions -minMax stdout |";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data   = split "\t", $string;
    my $signal = $Data[7];
    $mean +=$signal;
    $count++;
    $max = $signal if $signal > $max;

    $signal = sprintf("%.6f",$signal);
    $Open_Chrom{$Data[0]} = "$signal";
    $peaks++ if $signal >= $Open_min;
}
close FILE_IN;
$mean = sprintf("%.8f",$mean/$count);
map { $sd += ($_ - $mean)**2 } values %Open_Chrom;
$sd   = sprintf("%.8f",sqrt($sd/$count));
$max  = sprintf("%.8f",$max);

print "Load $Open_bw\n$mean\t$sd\t$max\n\n";
print FILE_LOG "Open chromatin = $Open_bw\n$mean\t$sd\t$max\npassed = $peaks\n\n";

######################
###    Cofactor    ###
######################
my %Cofactor = ();
my ($cof_mean,$cof_count,$cof_sd,$cof_max) = (0,0,0,0);
$bigWig     = "$data_folder/BigWig/$Cofactor_bw.bigWig";
open FILE_IN, "bigWigAverageOverBed $bigWig $data_folder/files/$regions -minMax stdout |";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my $signal = $Data[7];
    $cof_mean += $signal;
    $cof_count++;
    $cof_max = $signal if $signal > $cof_max;

    $signal = sprintf("%.8f",$signal);
    $Cofactor{$Data[0]} = "$signal";
}
close FILE_IN;
$cof_mean = sprintf("%.8f",$cof_mean/$cof_count);
map { $cof_sd += ($_ - $cof_mean)**2 } values %Cofactor;
$cof_sd   = sprintf("%.8f",sqrt($cof_sd/$cof_count));
$cof_max  = sprintf("%.8f",$cof_max);

print "Load $Cofactor_bw\n$cof_mean\t$cof_sd\t$cof_max\n\n";
print FILE_LOG "Cofactor = $Cofactor_bw\n$cof_mean\t$cof_sd\t$cof_max\n\n";

###################
###    CTCF     ###
###################
my %CTCF_signal = ();
($mean,$count,$sd,$max) = (0,0,0,0);
$bigWig = "$data_folder/BigWig/$CTCF_bw.bigWig";
open FILE_IN, "bigWigAverageOverBed $bigWig $data_folder/files/$regions -minMax stdout |";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my $signal = $Data[7];
    $mean +=$signal;
    $count++;
    $max = $signal if $signal > $max;

    $signal = sprintf("%.8f",$signal);
    $CTCF_signal{$Data[0]} = "$signal";
}
close FILE_IN;
$mean = sprintf("%.8f",$mean/$count);
map { $sd += ($_ - $mean)**2 } values %CTCF_signal;
$sd   = sprintf("%.8f",sqrt($sd/$count));
$max  = sprintf("%.8f",$max);

print "Load $CTCF_bw\n$mean\t$sd\t$max\n\n";
print FILE_LOG "CTCF = $CTCF_bw\n$mean\t$sd\t$max\n\n";

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
    $mean +=$signal;
    $count++;
    $max = $signal if $signal > $max;

    $signal = sprintf("%.8f",$signal);
    $Consensus{$Data[0]} = "$signal";
}
close FILE_IN;
$mean = sprintf("%.8f",$mean/$count);
map { $sd += ($_ - $mean)**2 } values %Cofactor;
$sd   = sprintf("%.8f",sqrt($sd/$count));
$max  = sprintf("%.8f",$max);

print "Load $tissue.H3K27ac.consensus\n$mean\t$sd\t$max\n\n";
print FILE_LOG "Consensus = $tissue.H3K27ac.consensus\n$mean\t$sd\t$max\n\n";

############################
### Altogether. Activity ###
############################
($mean,$count,$sd,$max) = (0,0,0,0);
my $non_conventional = 0;
my @List = ();

open FILE_EXPR,">Expression.data";
print FILE_EXPR "ensembl_gene_id\tchrom\ttype\texpression\tsymbol\n";

open FILE_OUT,"| pigz -p2 > tissue.regions.data.gz";
print FILE_OUT "#chrom\tstart\tend\tname\tscore\tstrand\ttype\tH3K27ac\tOpen\tCofactor\tCTCF\tActivity\tStatus\n";
my %Status = ();

open FILE_IN,"$data_folder/files/$regions";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;

    my $H3K27ac  = $Consensus{$Data[3]};
    my $Open     = $Open_Chrom{$Data[3]};
    my $CTCF     = $CTCF_signal{$Data[3]};
    my $Cofactor = $Cofactor{$Data[3]};

###################################
####  Protein-coding promoter. ####
###################################
    if($Data[3] =~ /^ENSG/ and $genes->{$Data[3]}->{"type"} eq "protein")
    {
	my $gene_id  = $Data[3];
	next unless defined $expression->{$gene_id};
	next unless $genes->{$gene_id};
	my $strand   = $genes->{$gene_id}->{"strand"};
	my $gene_exp = $expression->{$gene_id};
	my $symbol   = $genes->{$gene_id}->{"Symbol"};
	print FILE_EXPR "$gene_id\t$Data[0]\tprotein\t$gene_exp\t$symbol\n";

	if($gene_exp >= $min_exp)
	{
	    my $pseudo_H3K27ac = $H3K27ac;
	    my $pseudo_Open    = $Open;
	    $pseudo_H3K27ac    = 5 if $pseudo_H3K27ac < 5;
	    $pseudo_Open       = 5 if $pseudo_Open < 5;
	    my $Combined       = $pseudo_H3K27ac + $Cofactor;
	    my $Activity       = sprintf("%.8f",($pseudo_Open*$Combined)**0.5);
	    print FILE_OUT join "\t", @Data[0..3];
	    print FILE_OUT "\t1000\t$strand\tprotein\t$H3K27ac\t$Open\t$Cofactor\t$CTCF\t$Activity\texpressed\n";
	    $mean += $Activity;
	    push @List,$Activity;
	    $max = $Activity if $max < $Activity;
	    $count++;
	    $Status{"expressed"}++;
	}
	else
	{
	    my $Combined = $H3K27ac + $Cofactor;
	    my $Activity = sprintf("%.8f",($Open*$Combined)**0.5);
	    print FILE_OUT join "\t", @Data[0..3];
	    print FILE_OUT "\t1000\t$strand\tprotein\t$H3K27ac\t$Open\t$Cofactor\t$CTCF\t$Activity\tnot_expressed\n";
	    $mean += $Activity;
	    push @List,$Activity;
	    $max = $Activity if $max < $Activity;
	    $count++;
	    $Status{"not_expressed"}++;
	}
	next;
    }

##################
####  lncRNA. ####
##################
    if($Data[3] =~ /^ENSG/ and $genes->{$Data[3]}->{"type"} eq "lncRNA")
    {
	my $gene_id  = $Data[3];
	next unless defined $expression->{$gene_id};
	next unless $genes->{$gene_id};
	my $strand   = $genes->{$gene_id}->{"strand"};
	my $gene_exp = $expression->{$gene_id};
	my $symbol   = $genes->{$gene_id}->{"Symbol"};
	print FILE_EXPR "$gene_id\t$Data[0]\tlncRNA\t$gene_exp\t$symbol\n";

	if($gene_exp >= $min_exp)
	{
	    my $pseudo_H3K27ac = $H3K27ac;
	    my $pseudo_Open    = $Open;
	    $pseudo_H3K27ac    = 5 if $pseudo_H3K27ac < 5;
	    $pseudo_Open       = 5 if $pseudo_Open < 5;
	    my $Combined       = $pseudo_H3K27ac + $Cofactor;
	    my $Activity       = sprintf("%.8f",($pseudo_Open*$Combined)**0.5);
	    print FILE_OUT join "\t", @Data[0..3];
	    print FILE_OUT "\t1000\t$strand\tlncRNA\t$H3K27ac\t$Open\t$Cofactor\t$CTCF\t$Activity\texpressed\n";
	    $mean += $Activity;
	    push @List,$Activity;
	    $max = $Activity if $max < $Activity;
	    $count++;
	    my $status = get_enhancer_status($H3K27ac,$CTCF);
	    $Status{"expressed"}++;
	    $Status{"enhancer"}++;
	    $Status{$status}++;
	}
	else
	{
	    my $status   = get_status($H3K27ac,$CTCF,$Open);
	    if($status eq "active")
	    {
#######---------Promoter of lncRNA is active enhancer------
		my $Combined = $H3K27ac + $Cofactor + $ctcf_impact*$CTCF;
		my $Run_Open = $Open;
		$Run_Open    = 1 if $Run_Open < 1; #### CTCF peaks may have low Open chromatin since they are difficult to displace. Normal enhancers have high Open by default.
		my $Activity = sprintf("%.8f",($Run_Open*$Combined)**0.5);
		my $status   = get_enhancer_status($H3K27ac,$CTCF);
		print FILE_OUT join "\t", @Data[0..3];
		print FILE_OUT "\t1000\t$strand\tlncRNA\t$H3K27ac\t$Open\t$Cofactor\t$CTCF\t$Activity\t$status\n";
		$mean += $Activity;
		push @List,$Activity;
		$max = $Activity if $max < $Activity;
		$Status{$status}++;
		$count++;
		$Status{"enhancer"}++;
		$non_conventional++ if $H3K27ac < $H3K27ac_min or $Open < $Open_min;
	    }
	    else
	    {
#######---------Promoter of lncRNA has not enough Open and/or (H3K27ac or CTCF)------
		my $Activity    = "0.0";
		print FILE_OUT join "\t", @Data[0..3];
		print FILE_OUT "\t1000\t$strand\tlncRNA\t$H3K27ac\t$Open\t$Cofactor\t$CTCF\t$Activity\tinactive\n";
		$count++;
		$Status{"inactive"}++;
	    }
	}
	next;
    }

#####################
####  Enhancers  ####
#####################
    my $status   = get_status($H3K27ac,$CTCF,$Open);
    if($status eq "active")
    {
#######---------enhancer is active------
	my $Combined = $H3K27ac + $Cofactor + $ctcf_impact*$CTCF;
	my $Run_Open = $Open;
	$Run_Open    = 1 if $Run_Open < 1;  #### CTCF peaks may have low Open chromatin since they are difficult to displace. Normal enhancers have high Open by default.
	my $Activity = sprintf("%.8f",($Run_Open*$Combined)**0.5);
	my $status   = get_enhancer_status($H3K27ac,$CTCF);
	print FILE_OUT join "\t", @Data[0..3];
	print FILE_OUT "\t1000\t.\tenhancer\t$H3K27ac\t$Open\t$Cofactor\t$CTCF\t$Activity\t$status\n";
	$mean += $Activity;
	push @List,$Activity;
	$max = $Activity if $max < $Activity;
	$count++;
	$Status{$status}++;
	$Status{"enhancer"}++;
	$non_conventional++ if $H3K27ac < $H3K27ac_min or $Open < $Open_min;
    }
    else
    {
#######-----Enhancer has not enough Open and/or (H3K27ac or CTCF)------
	my $Activity = "0.0";
	print FILE_OUT join "\t", @Data[0..3];
	print FILE_OUT "\t1000\t.\tenhancer\t$H3K27ac\t$Open\t$Cofactor\t$CTCF\t$Activity\tinactive\n";
	$count++;
	$Status{"inactive"}++;
    }
}
close FILE_IN;
close FILE_OUT;
close FILE_EXPR;
$mean = sprintf("%.6f",$mean/$count);
map { $sd += ($_ - $mean)**2 } values @List;
$sd   = sprintf("%.6f",sqrt($sd/$count));
$max  = sprintf("%.6f",$max);

print FILE_LOG "Activity = $tissue\n$mean\t$sd\t$max\n\n";
print FILE_LOG "expressed = ",$Status{"expressed"},"\n";
print FILE_LOG "enhancers = ",$Status{"enhancer"},"\n\n";

print FILE_LOG "H3K27ac_enhancers = ",$Status{"H3K27ac_enhancer"},"\n";
print FILE_LOG "CTCF_enhancers    = ",$Status{"CTCF_enhancer"},"\n";
print FILE_LOG "Interim_enhancers = ",$Status{"Interim_enhancer"},"\n\n";
print FILE_LOG "Non_conventional  = $non_conventional\n";

close FILE_LOG;

system "mv -f tissue.regions.data.gz $data_folder/data/$tissue.regions.data.gz";
system "mv -f Annotate_regions.log   $data_folder/data/$tissue.signals.log";
system "mv -f Expression.data        $data_folder/data/$tissue.expression.data";

sub get_enhancer_status
{
    my ($H3K27ac,$CTCF) = @_;
    return "H3K27ac_enhancer" if ($H3K27ac+1) / ($CTCF+1)    >= 2.5;
    return "CTCF_enhancer"    if ($CTCF+1)    / ($H3K27ac+1) >= 2.5;
    return "Interim_enhancer";
}

sub get_status
{
    my ($H3K27ac,$CTCF,$Open) = @_;
    return "active"  if $CTCF >= $CTCF_min;

    return "nothing" if $Open < $Open_min;
    return "nothing" if $H3K27ac < $H3K27ac_min;

    return "active";
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
	$ret_h{$Data[0]} = $Data[6];
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
