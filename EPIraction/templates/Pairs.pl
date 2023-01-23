#!/usr/bin/perl -w
use strict;

my $data_folder   = "!{data_folder}";
my $regions       = "!{regions}";
my $chrom         = "!{chrom}";
die "wrong folder $data_folder\n" unless -d $data_folder;

if(-f "$data_folder/data/zzz.Pairs.$chrom.gz")
{
    print "already $data_folder/data/zzz.Pairs.$chrom.gz\n";
    exit;
}

open FILE_IN,"$data_folder/files/$regions";
open PROMOTERS,">Promoters.bed";
open ENHANCERS,">Enhancers.bed";
while(my $string = <FILE_IN>)
{
    my @Data = split "\t", $string;
    next unless $Data[0] eq $chrom;
    if($Data[3] =~ /^ENSG/)
    {
	print PROMOTERS $string;
    }
    else
    {
	print ENHANCERS $string;
    }
}
close FILE_IN;
close ENHANCERS;
close PROMOTERS;
system "bedtools window -a Promoters.bed -b Enhancers.bed -w 1000000 | gzip > Overlap.gz";

open FILE_OUT,"| gzip > Pairs.gz";
print FILE_OUT "promoter\tcount\tlist\n";

open FILE_IN,"zcat Overlap.gz |";
my $string = <FILE_IN>;chomp $string;
my @Data = split "\t", $string;
my ($genes,$pairs,$all_genes) = (0,0,0);
my %All_genes = ();

my @Pairs = ();push @Pairs,\@Data;

while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Current = split "\t", $string;
    $All_genes{$Current[3]}++;
    if($Current[3] eq $Pairs[0]->[3])
    {
        push @Pairs,\@Current;
    }
    else
    {
	$all_genes++;
	&report_pairs();
	@Pairs = ();
	push @Pairs,\@Current;
    }
}
close FILE_IN;
&report_pairs();
close FILE_OUT;
$all_genes++;
my @AAA = keys %All_genes;
print "Done $chrom\nGenes=$genes ($all_genes) (",(1+$#AAA),")\n";
print "Pairs=$pairs\n";
print sprintf("Average=%.1f\n",$pairs/$genes);
system "mv Pairs.gz $data_folder/data/zzz.Pairs.$chrom.gz";
system "rm -f Promoters.bed Enhancers.bed Overlap.gz";

sub report_pairs
{
    @Pairs = sort { $a->[6] <=> $b->[6] } @Pairs;
    my $promoter = $Pairs[0]->[3];
    print FILE_OUT "$promoter\t",(1+$#Pairs),"\t",(join ";", map { "$_->[7]" } @Pairs),"\n";
    $genes++;
    $pairs+=1+$#Pairs;
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
