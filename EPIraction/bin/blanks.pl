#!/usr/bin/perl -w
use strict;

my $string = <STDIN>;
my @Data = split "\t", $string;
print $string;
my $promoter = $Data[0];

while($string = <STDIN>)
{
    my @Data = split "\t", $string;
    if($Data[0] ne $promoter)
    {
	print "\n";
	$promoter = $Data[0];
    }
    print $string;
}
