#!/usr/bin/perl -w
use strict;

my $data_folder = "!{data_folder}";
my $chromosome  = "!{chrom}";
my $tissues     = "!{tissues}";
my $regions     = "!{regions}";
die "wrong folder $data_folder\n" unless -d $data_folder;

my ($genes,undef) = load_data("$data_folder/files/Gencode.v40.genes.data");

if(-f "$data_folder/data/xxx.Contacts.$chromosome.data.gz")
{
    print "Already $data_folder/data/xxx.Contacts.$chromosome.data.gz\n";
}

$tissues =~ s/\[|\]//g;
my @Tissues = sort split ", ", $tissues;
my $positions   = load_regions("$data_folder/files/$regions");

my %Contacts = ();
foreach my $tissue (@Tissues)
{
    open FILE_IN,"unpigz -p3 -c $data_folder/data/$tissue.contacts.data.gz |";
    my $string = <FILE_IN>;
    while(my $string = <FILE_IN>)
    {
	chomp $string;
	next if $string eq "";
	next if $string =~ /^##/;
	my ($chrom,$promoter,$enhancer,$distance,$activity,$contact,$rescaled,$foldchange) = split "\t", $string;
	next unless $chrom eq $chromosome;
	$enhancer =~s/\s+//;
	$Contacts{"$tissue<>$promoter<>$enhancer"} = sprintf("%.8f",$contact);
    }
    close FILE_IN;
}

open FILE_OUT,"| pigz -p3 > Contacts.data.gz";
print FILE_OUT "promoter\tenhancer\tdistance\t",(join "\t", @Tissues),"\n";

open FILE_IN,"unpigz -p3 -c $data_folder/data/zzz.Pairs.$chromosome.gz |";
my $string = <FILE_IN>;
while($string = <FILE_IN>)
{
    chomp $string;
    my ($promoter,undef,$list) = split "\t", $string;
    my ($p_start,$p_end) = @{$positions->{$promoter}};
    foreach my $enhancer (split ";", $list)
    {
	my ($e_start,$e_end) = @{$positions->{$enhancer}};
	my @Sorted       = sort { $a <=> $b } ($e_start,$e_end,$p_start,$p_end);
	my $distance     = $Sorted[2] - $Sorted[1];
	$distance        = 0 if $enhancer eq $promoter;

#########################################################################################
### Contrary to conventional ABC I ignore the promoters of other protein-coding genes ###
#########################################################################################
	next if $enhancer =~/^ENSG/ and $enhancer ne $promoter and $genes->{$enhancer}->{"type"} eq "protein";

	my @Report = ();
	foreach my $tissue (@Tissues)
	{
	    if(defined $Contacts{"$tissue<>$promoter<>$enhancer"})
	    {
		push @Report,$Contacts{"$tissue<>$promoter<>$enhancer"};
	    }
	    else
	    {
		push @Report,"0.000000";
	    }
	}
	print FILE_OUT "$promoter\t$enhancer\t$distance\t",(join "\t", @Report),"\n";
    }
    print FILE_OUT "\n";
}
close FILE_IN;
close FILE_OUT;
system "mv Contacts.data.gz $data_folder/data/xxx.Contacts.$chromosome.data.gz";
print "Done $chromosome\n";

sub load_regions
{
    my %ret_h = ();
    open FILE_IN,"$_[0]";
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$ret_h{$Data[3]} = [$Data[1],$Data[2]];
    }
    close FILE_IN;
    return\%ret_h;
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
