#!/usr/bin/perl
use strict;

my $data_folder = "!{data_folder}";
my $tissues     = "!{tissues}";
my $blocks      = "!{blocks}";
my $minimum_exp = "!{minimum_exp}";

die "wrong folder $data_folder\n" unless -d $data_folder;

$tissues =~ s/\[|\]//g;
my @Tissues = sort split ", ", $tissues;
my (%Expressed,%Total,%Summary) = ();
foreach my $tissue (@Tissues)
{
    open FILE_IN,"$data_folder/data/$tissue.expression.data";
    my $string = <FILE_IN>;
    chomp $string;
    my @Header = split "\t", $string;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my ($gene_id,$chrom,$type,$expression,$symbol) = split "\t", $string;
	next if $expression < $minimum_exp;
	$Expressed{"$tissue<>$gene_id"} = $expression;
	$Total{$gene_id}++;
	$Summary{$gene_id} += $expression;
    }
    close FILE_IN;
    print "Load $tissue\n";
}

my ($genes,undef) = load_data("$data_folder/files/Gencode.v40.genes.data");
my %Chromosomes = ();
map { push @{$Chromosomes{$genes->{$_}->{"chrom"}}},$_ } sort keys %{$genes};

open FILE_MATRIX,">Expressed.genes.matrix";
print FILE_MATRIX "ensembl_gene\tsymbol\texpressed\t",(join "\t", @Tissues),"\n";
my @List = ();
foreach my $chrom (sort keys %Chromosomes)
{
    my @Report = ();
    my $sum = 0;
    foreach my $gene_id (@{$Chromosomes{$chrom}})
    {
	next unless $Total{$gene_id};
	print FILE_MATRIX $gene_id,"\t",$genes->{$gene_id}->{"Symbol"},"\t",$Total{$gene_id},"\t",(join "\t", map { defined $Expressed{"$_<>$gene_id"} ? $Expressed{"$_<>$gene_id"} : 0 } @Tissues),"\n";

	push @Report,$gene_id;
	$sum+=$Summary{$gene_id};
	if($#Report >= $blocks-1)
	{
	    my $report = "$chrom\t".(join ";", @Report);
	    push @List,[$report,$sum];
	    @Report = ();
	    $sum=0;
	}
    }
    if($#Report >= 0)
    {
	my $report = "$chrom\t".(join ";", @Report);
	push @List,[$report,$sum];
    }
}
close FILE_MATRIX;
system "mv Expressed.genes.matrix $data_folder/data/www.Expressed.genes.matrix";

open FILE_LIST,">Expressed.genes.tsv";
print FILE_LIST "group\tlist\tscore\n";
map { print FILE_LIST $_->[0],"\t",$_->[1],"\n" } @List;
close FILE_LIST;

print "Done\n";

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
