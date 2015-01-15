#!/usr/bin/perl

use lib "/home/msamanta/Z-Codes";
use Useful;

open(IN,"kmers-21");
while(<IN>)
{
	@s=split(/\s+/,$_);
	print "$s[0]\n";
	$x=Useful::reverse($s[0]);
	print "$x\n";
}

