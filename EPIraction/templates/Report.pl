#!/usr/bin/perl -w
use strict;

my $data_folder  = "!{data_folder}";
my $tissue       = "!{tissue}";
my $genome_file  = "!{genome_file}";
my $interact     = "!{interact}";
my $baseline_hic = "!{HiC_baseline}";

die "wrong folder $data_folder\n" unless -d $data_folder;

my $HiC_baseline  = load_baseline();
my ($genes,undef) = load_data("$data_folder/files/Gencode.v40.promoters.data");

open PROMOTERS,"| gzip > Genes.bed.gz";
print PROMOTERS "#chrom\tstart\tend\tgene_id\tscore\tstrand\tSymbol\tTSS\tTPMs\tH3K27ac\tOpen\tActivity\tsample\n";

open ENHANCERS,"| gzip > Enhancers.bed.gz";
print ENHANCERS "#chrom\tstart\tend\tenhancer\tscore\tstrand\tH3K27ac\tOpen\tCofactor\tActivity\tsample\n";

my (%Promoters,%Enhancers) = ();
open FILE_IN,"$data_folder/data/$tissue.regions.data";
my $string = <FILE_IN>;
while($string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    if($Data[3]=~/^ENSG/)
    {
	my $gene_data = $genes->{$Data[3]};
	print PROMOTERS join "\t", @Data[0..5];
	print PROMOTERS "\t",$gene_data->{"Symbol"},"\t",$gene_data->{"TSS"},"\t",(join "\t", @Data[10,7,8,13,12]),"\n";
    }
    else
    {
	next if $Data[13] == 0;
	my $signals  = join "\t",@Data[7,8,9,13];
	my $position = join "\t",@Data[0..3];
	print ENHANCERS join "\t", @Data[0..5];
	print ENHANCERS "\t$signals\t$tissue\n";
	$Enhancers{$Data[3]} = [$position,$signals];
    }
}
close FILE_IN;
close PROMOTERS;
close ENHANCERS;

open FILE_IN,"zcat $data_folder/data/$tissue.contacts.data.gz |";
my %HiC = ();
$string = <FILE_IN>;
while($string = <FILE_IN>)
{
    chomp $string;
    next if $string eq "";
    my @Data = split "\t", $string;
    $HiC{"$Data[2]<>$Data[3]"} = [$Data[4],sprintf("%.5f\t%.5f",$Data[6],$Data[8])];
}
close FILE_IN;

open FILE_OUT,"| gzip > Pairs.bed.gz";
print FILE_OUT "#chr\tstart\tend\tname\tclass\tTargetGene\tTargetGeneEnsemblID\tTargetGeneTSS\tCellType\tScore\tDistanceToTSS\tH3K27ac\tOpen\tCofactor\tActivity\tHiC_Contacts\tHiC_FoldChange\n";

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

my $promoters = load_promoters();

open FILE_IN,"zcat Pairs.bed.gz |";
$string = <FILE_IN>;
chomp $string;
my @Header = split "\t", $string;
open FILE_OUT,"| gzip > Threshold.bed.gz";
print FILE_OUT $string,"\n";
open FILE_INTERACT,"| LC_ALL=C sort --parallel=1 --buffer-size=1G -k 1,1d -k 14,14d -k 2,2n -k 15,15n - --temporary-directory='./' > Pairs.interact";
my @Pairs = ();
while($string = <FILE_IN>)
{
    chomp $string;
    if($string =~/^\#/)
    {
	&report_pair();
	@Pairs = ();
    }
    else
    {
	my @Data = split "\t", $string;
	push @Pairs,\@Data;
    }
}
close FILE_IN;
&report_pair();
close FILE_OUT;
close FILE_INTERACT;

system "bedToBigBed -as=$interact -type=bed5+13 Pairs.interact $genome_file Threshold.bb";

open FILE_IN,"Pairs.interact";
open FILE_OUT,"| gzip > Threshold.bedpe.gz";
print FILE_OUT "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tEPIraction_score\tcolor\n";
while($string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    print FILE_OUT join "\t", @Data[8,9,10,13,14,15];
    my $score = sprintf("%d",$Data[5]*5000);
    $score = 1000 if $score > 1000;
    if($Data[16] =~/^ENSG/)
    {
	print FILE_OUT "\t$Data[11]:vs:$Data[16]\t$score\t.\t.\t$Data[5]\t$Data[7]\n";
    }
    else
    {
	print FILE_OUT "\t$Data[16]:vs:$Data[11]\t$score\t.\t.\t$Data[5]\t$Data[7]\n";
    }
}
close FILE_IN;
close FILE_OUT;

system "mv Genes.bed.gz      $data_folder/report/EPIraction.$tissue.genes.bed.gz";
system "mv Enhancers.bed.gz  $data_folder/report/EPIraction.$tissue.enhancers.bed.gz";
system "mv Pairs.bed.gz      $data_folder/report/EPIraction.$tissue.complete.pairs.bed.gz";
system "mv Threshold.bed.gz  $data_folder/report/EPIraction.$tissue.threshold.pairs.bed.gz";

system "mv Threshold.bedpe.gz $data_folder/report/EPIraction.$tissue.threshold.pairs.bedpe.gz";
system "mv Threshold.bb       $data_folder/report/EPIraction.$tissue.threshold.pairs.bb";

system "rm -f Pairs.interact";

sub report_pair
{
    return if $#Pairs == -1;
    my @Subsantial = grep { $_->[9] >= 0.01 } @Pairs;
    return if $#Pairs == -1;
    foreach my $pair_data (@Subsantial)
    {
	print FILE_OUT join "\t", @{$pair_data};
	print FILE_OUT "\n";

	my $color = "#cacfd3";
	my $score = sprintf("%d",$pair_data->[9]*5000);
	$score = 1000 if $score > 1000;

	if($pair_data->[9] > 0.15)
	{
	    $color = "#000000";
	}
	elsif($pair_data->[9] > 0.05)
	{
	    $color = "#1e8450";
	}
	
	if($pair_data->[1] < $promoters->{$pair_data->[6]}->[1])
	{
	    print FILE_INTERACT "$pair_data->[0]\t$pair_data->[1]\t$pair_data->[2]\t$pair_data->[3]\t$score\t$pair_data->[9]\t$tissue\t$color\t";
	    print FILE_INTERACT "$pair_data->[0]\t$pair_data->[1]\t$pair_data->[2]\t$pair_data->[3]\t.\t",$promoters->{$pair_data->[6]}->[0],"\t$pair_data->[6]\t.\n";
	}
	else
	{
	    print FILE_INTERACT $promoters->{$pair_data->[6]}->[0],"\t$pair_data->[6]\t$score\t$pair_data->[9]\t$tissue\t$color\t";
	    print FILE_INTERACT $promoters->{$pair_data->[6]}->[0],"\t$pair_data->[6]\t.\t$pair_data->[0]\t$pair_data->[1]\t$pair_data->[2]\t$pair_data->[3]\t.\n";
	}
    }
    print FILE_OUT "#\n";
}

sub load_promoters
{
    my %ret_h = ();
    open FILE_IN,"$data_folder/files/Gencode.v40.promoters.data";
    my $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$ret_h{"$Data[0]"} = ["$Data[1]\t".($Data[3]-200)."\t".($Data[3]+200),$Data[3]];
    }
    close FILE_IN;
    return \%ret_h;
}

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
    open FILE_IN,$baseline_hic or die "Absent $baseline_hic";
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

sub get_200_distance
{
    my $value = $_[0];
    return 1000000 if $value >= 999900;
    my $times = int(($value+100)/200);
    return $times*200;
}
