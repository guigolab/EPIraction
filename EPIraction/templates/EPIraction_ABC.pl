#!/usr/bin/perl -w
use strict;

my $data_folder  = "!{data_folder}";
my $tissue       = "!{tissue}";
my $minimum_exp  = "!{minimum_exp}";
my $promoter_minimum = 0.05;

die "wrong folder $data_folder\n" unless -d $data_folder;

my (%Regions,%Contacts,@Genes) = ();
my ($genes,undef) = load_data("$data_folder/files/Gencode.v40.genes.data");

open FILE_IN,"$data_folder/data/$tissue.expression.data";
my $string = <FILE_IN>;
chomp $string;
my @Header = split "\t", $string;
while($string = <FILE_IN>)
{
    chomp $string;
    my ($gene_id,$chrom,$type,$expression,$symbol) = split "\t", $string;
    next if $expression < $minimum_exp;
    push @Genes,$gene_id;
}
close FILE_IN;
print "Load expressed\n";

open FILE_IN,"unpigz -p3 -c $data_folder/data/$tissue.regions.data.gz |";
$string = <FILE_IN>;
chomp $string;
@Header = split "\t", $string;

while($string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[11] == 0;
    map { $Data[$_] = sprintf("%.4f",$Data[$_]); } (7..10);
    $Regions{$Data[3]} = \@Data;
}
close FILE_IN;
print "Load regions\n";

open FILE_IN,"unpigz -p3 -c $data_folder/data/$tissue.contacts.data.gz |";
$string = <FILE_IN>;
chomp $string;
@Header = split "\t", $string;

while($string = <FILE_IN>)
{
    next if $string =~/^##/;
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[6] == 0;
    $Data[2] =~s/^\s+//;
    next if $Data[2] ne $Data[1] and $Data[2] =~ /^ENSG/ and $genes->{$Data[2]}->{"type"} eq "protein";

    push @{$Contacts{"$Data[1]"}}, ["$Data[2]","$Data[3]","$Data[6]","$Data[5]","$Data[7]"];
}
close FILE_IN;
print "Load contacts\n";

open FILE_OUT,"| gzip > Pairs.bed.gz";
print FILE_OUT "#chr\tstart\tend\tname\tclass\tTargetGene\tTargetGeneEnsemblID\tTargetGeneTSS\tCellType\tScore\tDistanceToTSS\tH3K27ac\tOpen\tCofactor\tCTCF\tActivity\tHiC_contacts\tHiC_foldchange\n";

foreach my $gene_id (@Genes)
{
    next unless $Contacts{"$gene_id"};
    my $gene_data = $Regions{$gene_id};
    my $TSS       = $genes->{$gene_id}->{"TSS"};
    my $symbol    = $genes->{$gene_id}->{"Symbol"};
    my @Report    = ();
    my $ABC_sum   = 0;
    my ($promoter_index,$promoter_abc) = ();

    foreach my $enhancer_data (@{$Contacts{"$gene_id"}})
    {
	my ($enhancer_id,$distance,$rescaled,$contact,$fold_change) = @{$enhancer_data};
	$ABC_sum += $rescaled;

	my $region = $Regions{$enhancer_id};
	my $type   = $region->[6];
	if($enhancer_id eq $gene_id)
	{
	    my $report_before = "$region->[0]\t$region->[1]\t$region->[2]\t$gene_id\tpromoter\t$symbol\t$gene_id\t$TSS\t$tissue";
	    my $report_after  = "$distance\t$region->[7]\t$region->[8]\t$region->[9]\t$region->[10]\t$region->[11]\t$contact\t$fold_change\n";
	    push @Report,[$report_before,$rescaled,$report_after];
	    $promoter_index = $#Report;
	    $promoter_abc = $rescaled;
	}
	else
	{
	    my $report_before = "$region->[0]\t$region->[1]\t$region->[2]\t$enhancer_id\t$type\t$symbol\t$gene_id\t$TSS\t$tissue";
	    my $report_after  = "$distance\t$region->[7]\t$region->[8]\t$region->[9]\t$region->[10]\t$region->[11]\t$contact\t$fold_change\n";
	    push @Report,[$report_before,$rescaled,$report_after];
	}
    }
    my $promoter_fraction = $promoter_abc/$ABC_sum;
    my $shift = ($promoter_minimum*$ABC_sum - $promoter_abc)/(1-$promoter_minimum);
    $shift = 0 if $shift < 0;
    $Report[$promoter_index]->[1] += $shift;
    $ABC_sum += $shift;
#    print "$gene_id promoter fraction: ",sprintf("%.6f",$promoter_fraction)," --- shift = $shift\n";

    foreach my $line (sort { $b->[1] <=> $a->[1] } @Report)
    {
	next if $line->[1] == 0;
	print FILE_OUT $line->[0],"\t",sprintf("%.12f",$line->[1]/$ABC_sum),"\t",$line->[2];
    }
    print FILE_OUT "#\n";
}
close FILE_OUT;
print "Done EPIraction.ABC\n";
system "mv Pairs.bed.gz $data_folder/report/EPIraction_ABC.$tissue.complete.pairs.bed.gz";

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
