#!/usr/bin/perl -w
use strict;

my $data_folder  = "!{data_folder}";
my $tissue       = "!{tissue}";
my $genome_file  = "!{genome_file}";
my $interact     = "!{interact}";
my $version      = "!{version}";
my $minimum_exp  = "!{minimum_exp}";
my $promoter_min = 0.05;
my $substantial  = 0.01;
srand(127);

die "wrong folder $data_folder\n" unless -d $data_folder;

my (%Regions,%Contacts,@Genes,%Promoters) = ();
my ($genes,undef) = load_data("$data_folder/files/Gencode.v40.genes.data");

open FILE_IN,"$data_folder/data/$tissue.expression.data";
my $string = <FILE_IN>;
chomp $string;
my @Header = split "\t", $string;
my %Expression = ();
while($string = <FILE_IN>)
{
    chomp $string;
    my ($gene_id,$chrom,$type,$expression,$symbol) = split "\t", $string;
    $Expression{$gene_id} = $expression;
    next if $expression < $minimum_exp;
    push @Genes,$gene_id;
}
close FILE_IN;
print "Load expressed\n";

open PROMOTERS,"> Genes.bed";
print PROMOTERS "#chrom\tstart\tend\tgene_id\tscore\tstrand\ttype\tSymbol\tTSS\tTPMs\tH3K27ac\tOpen\tCofactor\tCTCF\tActivity\tsample\tversion\n";

open ENHANCERS,">Enhancers.bed";
print ENHANCERS "#chrom\tstart\tend\tenhancer\tscore\tstrand\ttype\tH3K27ac\tOpen\tCofactor\tCTCF\tActivity\tsample\tversion\n";

open FILE_IN,"unpigz -p3 -c $data_folder/data/$tissue.regions.data.gz |";
$string = <FILE_IN>;
chomp $string;
@Header = split "\t", $string;

while($string = <FILE_IN>)
{
    chomp $string;
    my @Data      = split "\t", $string;
    my $signals   = sprintf("%.4f\t%.4f\t%.4f\t%.4f\t%.8f",$Data[7],$Data[8],$Data[9],$Data[10],$Data[11]);
    my $position  = join "\t", @Data[0,1,2,3,6];
    my $region_id = $Data[3];
    $Regions{$region_id} = [$position,$signals];

    if($region_id =~ /^ENSG/)
    {
	my $gene_data  = $genes->{$region_id};
	my $expression = $Expression{$region_id};
	if($gene_data->{"type"} eq "protein" and $expression >= 1)
	{
	    print PROMOTERS join "\t", @Data[0..5];
	    print PROMOTERS "\t",$gene_data->{"type"},"\t",$gene_data->{"Symbol"},"\t",$gene_data->{"TSS"},"\t$expression\t$signals\t$tissue\t$version\n";
	    $Promoters{$Data[3]} = ["$Data[0]","$Data[1]","$Data[2]"];
	}
	elsif($gene_data->{"type"} eq "lncRNA" and $expression >= 1)
	{
	    print PROMOTERS join "\t", @Data[0..5];
	    print PROMOTERS "\t",$gene_data->{"type"},"\t",$gene_data->{"Symbol"},"\t",$gene_data->{"TSS"},"\t$expression\t$signals\t$tissue\t$version\n";
	    $Promoters{$Data[3]} = ["$Data[0]","$Data[1]","$Data[2]"];

	    print ENHANCERS join "\t", @Data[0..6];
	    print ENHANCERS "\t$signals\t$tissue\t$version\n";
	}
	elsif($gene_data->{"type"} eq "lncRNA" and $Data[11] > 0)
	{
	    print ENHANCERS join "\t", @Data[0..6];
	    print ENHANCERS "\t$signals\t$tissue\t$version\n";
	}
    }
    else
    {
	next if $Data[11] == 0;
	print ENHANCERS join "\t", @Data[0..6];
	print ENHANCERS "\t$signals\t$tissue\t$version\n";
	$Regions{$Data[3]} = [$position,$signals];
    }
}
close FILE_IN;
close PROMOTERS;
close ENHANCERS;
print "Load regions\n";

open FILE_IN,"unpigz -p3 -c $data_folder/data/$tissue.contacts.data.gz |";
$string = <FILE_IN>;
chomp $string;
@Header = split "\t", $string;

while($string = <FILE_IN>)
{
    next if $string =~/^#/;
    chomp $string;
    my @Data = split "\t", $string;
    next if $Expression{$Data[1]} < $minimum_exp;
    next if $Data[4] == 0;
    $Data[2] =~s/^\s+//;
    $Contacts{"$Data[1]<>$Data[2]"} = ["$Data[3]","$Data[5]\t$Data[7]"];
}
close FILE_IN;
print "Load contacts\n";

open FILE_OUT,">Pairs.bed";
print FILE_OUT "#chr\tstart\tend\tname\tclass\tTargetGene\tTargetGeneEnsemblID\tTargetGeneTSS\tCellType\tScore\tDistanceToTSS\tH3K27ac\tOpen\tCofactor\tCTCF\tActivity\treward\tabc_tissue\tabc_closest\tabc_complete\tHiC_contacts\tHiC_foldchange\tVersion\n";

opendir DIR,"$data_folder/temp/batches";
foreach my $file (sort { compare_files($a,$b) } grep { /$tissue/ }  grep{ /pairs.data.gz/} readdir DIR)
{
    open FILE_IN,"zcat $data_folder/temp/batches/$file |";
    $string = <FILE_IN>;
    my @List = ();
    while($string = <FILE_IN>)
    {
	chomp $string;
	if($string =~/^#/)
	{
	    report_gene(\@List);
	    @List = ();
	    next;
	}
	else
	{
	    my @Data = split "\t", $string;
	    $Data[1] =~ s/^\s+//;
	    push @List,\@Data;
	}
    }
    close FILE_IN;
    print "Done $file\n";
}
close FILE_OUT;
print "Done EPIraction complete\n";

open FILE_IN,"Pairs.bed";
$string = <FILE_IN>;
chomp $string;
@Header = split "\t", $string;
open FILE_OUT,">Threshold.bed";
print FILE_OUT "#chr\tstart\tend\tname\tclass\tTargetGene\tTargetGeneEnsemblID\tTargetGeneTSS\tCellType\tScore\tDistanceToTSS\tH3K27ac\tOpen\tCofactor\tCTCF\tActivity\tHiC_contacts\tVersion\n";

open FILE_INTERACT,"| LC_ALL=C sort --parallel=1 --buffer-size=1G -k 1,1d -k 14,14d -k 2,2n -k 15,15n - --temporary-directory='./' > Pairs.interact";
my @Pairs = ();
while($string = <FILE_IN>)
{
    chomp $string;
    if($string =~/^\#/)
    {
	&report_threshold(\@Pairs);
	@Pairs = ();
    }
    else
    {
	my @Data = split "\t", $string;
	push @Pairs,\@Data;
    }
}
close FILE_IN;
close FILE_OUT;
close FILE_INTERACT;
print "Done EPIraction threshold\n";

system "bedToBigBed -as=$interact -type=bed5+13 Pairs.interact $genome_file Threshold.bb";

open FILE_IN,"Pairs.interact";
open FILE_OUT,">Threshold.bedpe";
open FILE_STRONG,">Pairs.strong";
print FILE_OUT "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tEPIraction_score\tcolor\n";
while($string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    print FILE_OUT join "\t", @Data[8,9,10,13,14,15,3,4];
    print FILE_OUT "\t.\t.\t$Data[5]\t$Data[7]\n";
    print FILE_STRONG $string,"\n" if $Data[5] >= 0.03;
}
close FILE_IN;
close FILE_OUT;
close FILE_STRONG;

system "bedToBigBed -as=$interact -type=bed5+13 Pairs.strong $genome_file Strong.bb";

system "pigz -p2 Pairs.bed";
system "pigz -p2 Genes.bed";
system "pigz -p2 Enhancers.bed";
system "pigz -p2 Threshold.bed";
system "pigz -p2 Threshold.bedpe";

system "mv Pairs.bed.gz       $data_folder/report/EPIraction.$tissue.complete.pairs.bed.gz";
system "mv Genes.bed.gz       $data_folder/report/EPIraction.$tissue.genes.bed.gz";
system "mv Enhancers.bed.gz   $data_folder/report/EPIraction.$tissue.enhancers.bed.gz";
system "mv Threshold.bed.gz   $data_folder/report/EPIraction.$tissue.threshold.pairs.bed.gz";
system "mv Threshold.bedpe.gz $data_folder/report/EPIraction.$tissue.threshold.pairs.bedpe.gz";
system "mv Threshold.bb       $data_folder/report/EPIraction.$tissue.threshold.pairs.bb";
system "mv Strong.bb          $data_folder/report/EPIraction.$tissue.strong.pairs.bb";
system "rm -f Pairs.interact";

sub report_threshold
{
    my $array_ref = $_[0];
    return if $#$array_ref == -1;
    my @Substantial = grep { $_->[9] >= $substantial } @{$array_ref};
    return if $#Substantial == -1;
    my $gene_id = $array_ref->[0]->[6];

    foreach my $pair_data (@Substantial)
    {
	print FILE_OUT join "\t", @{$pair_data}[(0..15),21];print FILE_OUT "\t$version\n";
	next if $pair_data->[4] eq "promoter";

	my $score = sprintf("%d",$pair_data->[9]*10000);
	$score = 1000 if $score > 1000;
	my $color = "";
	if($pair_data->[9] > 0.10)
	{
	    $color = "#00b400";
	}
	elsif($pair_data->[9] > 0.05)
	{
	    $color = "#1e8450";
	}
	else
	{
	    $color = "#8e7cc3";
	}
	my $promoter_bed4 = $Promoters{$gene_id}->[0]."\t".$Promoters{$gene_id}->[1]."\t".$Promoters{$gene_id}->[2]."\t$gene_id";
	my $enhancer_bed4 = "$pair_data->[0]\t$pair_data->[1]\t$pair_data->[2]\t$pair_data->[3]";
	
	if($pair_data->[1] < $Promoters{$gene_id}->[1])
	{
	    my $complete_bed = "$pair_data->[0]\t$pair_data->[1]\t".$Promoters{$gene_id}->[2]."\t$pair_data->[3]:$gene_id";
	    print FILE_INTERACT "$complete_bed\t$score\t$pair_data->[9]\t$tissue\t$color\t$enhancer_bed4\t.\t$promoter_bed4\t.\n";
	}
	else
	{
	    my $complete_bed = "$pair_data->[0]\t".$Promoters{$gene_id}->[1]."\t$pair_data->[2]\t$pair_data->[3]:$gene_id";
	    print FILE_INTERACT "$complete_bed\t$score\t$pair_data->[9]\t$tissue\t$color\t$promoter_bed4\t.\t$enhancer_bed4\t.\n";
	}
    }
    print FILE_OUT "###\n";
}

sub report_gene
{
    my $array_ref = $_[0];
    return if $#$array_ref == -1;
    my $gene_id   = $array_ref->[0]->[0];
    my $TSS       = $genes->{$gene_id}->{"TSS"};
    my $symbol    = $genes->{$gene_id}->{"Symbol"};
    my $chrom     = $genes->{$gene_id}->{"chrom"};

    my ($promoter_index,$promoter_abc,$sum_abc) = ();
    my $shift_promoter = "not";

    foreach my $index (0..$#$array_ref)
    {
	my $feature_id   = $array_ref->[$index]->[1];
	my $abc_tissue   = $array_ref->[$index]->[2];
	my $abc_closest  = $array_ref->[$index]->[3];
	my $abc_complete = $array_ref->[$index]->[4];
	$sum_abc        += $abc_tissue;

	if($feature_id eq $gene_id)
	{
	    $promoter_index = $index;
	    $promoter_abc   = $abc_tissue;
	    push @{$array_ref->[$index]},($abc_tissue,0);
	}
	else
	{
	    my $reward     = $abc_closest - $abc_complete;
	    $reward        = 0 if $reward < 0;
	    $reward        = $abc_tissue if $reward > $abc_tissue;
	    push @{$array_ref->[$index]},($abc_tissue,$reward);
	}
    }

    my $promoter_fraction = $promoter_abc/$sum_abc;
    if($promoter_fraction < $promoter_min)
    {
	my $shift = ($promoter_min*$sum_abc - $promoter_abc)/(1-$promoter_min);
	$array_ref->[$promoter_index]->[5] += $shift;
	$sum_abc += $shift;
    }

    foreach my $line (sort { $b->[5] <=> $a->[5] } @{$array_ref})
    {
	my $enhancer   = $line->[1];
	my $pos_report = $Regions{$enhancer}->[0];
	if($enhancer eq $gene_id)
	{
	    my @Position_bed = split "\t",$Regions{$gene_id}->[0];
	    $Position_bed[-1] = "promoter";
	    $pos_report = join "\t", @Position_bed;
	}
	my $abc_tissue   = $line->[2];
	my $abc_closest  = $line->[3];
	my $abc_complete = $line->[4];

	my $base_score   = $line->[5]/$sum_abc;
	my $reward       = $line->[6]/$sum_abc;
	$reward          = 0.1 if $reward > 0.1;
	$reward          = sprintf("%.10f",$reward);
	my $score        = sprintf("%.10f",$base_score+$reward);

	next if $score == 0;
	my $contact_data = [0,"1.000000\t1.00"];
	$contact_data = $Contacts{"$gene_id<>$enhancer"} if defined $Contacts{"$gene_id<>$enhancer"};
	print FILE_OUT "$pos_report\t$symbol\t$gene_id\t$TSS\t$tissue\t$score\t$contact_data->[0]\t",$Regions{$enhancer}->[1],"\t$reward\t$abc_tissue\t$abc_closest\t$abc_complete\t$contact_data->[1]\t$version\n";
    }
    print FILE_OUT "###\n";
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

sub get_200_distance
{
    my $value = $_[0];
    return 1000000 if $value >= 999900;
    my $times = int(($value+100)/200);
    return $times*200;
}

sub compare_files
{
    my ($AAA,$BBB) = @_;
    my @AAA = split /\./, $AAA;
    my @BBB = split /\./, $BBB;
    return $AAA[-4] <=> $BBB[-4];
}
