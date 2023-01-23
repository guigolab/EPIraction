#!/usr/bin/perl -w
use strict;

my $data_folder = "!{data_folder}";
my $tissues     = "!{tissues}";
my $blocks      = "!{blocks}";

die "wrong folder $data_folder\n" unless -d $data_folder;
$tissues =~ s/\[|\]//g;
my @Tissues = sort split ", ", $tissues;
my (%Expressed,%Total) = ();
foreach my $tissue (@Tissues)
{
    open FILE_IN,"$data_folder/data/$tissue.regions.data";
    my $string = <FILE_IN>;
    chomp $string;
    my @Header = split "\t", $string;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	if($Data[3]=~/^ENSG/)
	{
	    $Expressed{"$tissue<>$Data[3]"} = $Data[10];
	    $Total{"$Data[3]"}++;
#	    map { print $_,"\t",$Header[$_],"\t",$Data[$_],"\n" } (0..$#Data);exit;
	}
    }
    close FILE_IN;
    print "Load $tissue\n";
}

my ($genes,undef) = load_data("$data_folder/files/Gencode.v40.promoters.data");
my %Chromosomes = ();
map { push @{$Chromosomes{$genes->{$_}->{"chrom"}}},$_ } sort keys %{$genes};

open FILE_MATRIX,">Expressed.genes.matrix";
print FILE_MATRIX "ensembl_gene\tsymbol\texpressed\t",(join "\t", @Tissues),"\n";

open FILE_LIST,">Expressed.genes.tsv";
print FILE_LIST "group\tlist\n";

foreach my $chrom (sort keys %Chromosomes)
{
    my @Report = ();
    foreach my $gene_id (sort @{$Chromosomes{$chrom}})
    {
	next unless $Total{$gene_id};
	next if     $Total{$gene_id} < 2;
	print FILE_MATRIX $gene_id,"\t",$genes->{$gene_id}->{"Symbol"},"\t",$Total{$gene_id},"\t",(join "\t", map { defined $Expressed{"$_<>$gene_id"} ? $Expressed{"$_<>$gene_id"} : 0 } @Tissues),"\n";

	push @Report,$gene_id;
	if($#Report >= $blocks-1)
	{
	    print FILE_LIST $chrom,"\t",(join ";", @Report),"\n";
	    @Report = ();
	}
    }
    if($#Report >= 0)
    {
	print FILE_LIST $chrom,"\t",(join ";", @Report),"\n";
    }
}
close FILE_LIST;
close FILE_MATRIX;
system "mv Expressed.genes.matrix $data_folder/data/yyy.Expressed.genes.matrix";
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
