#!/usr/bin/perl

use lib "/home/msamanta/Z-Codes";
use Useful;

open(IN,$ARGV[0]);
while(<IN>)
{
	@s=split(/\s+/,$_);
	$x=Useful::reverse($s[0]);
	$rs0=$s[0]; $rx=$x;
	$s[0]=~s/A/0/g; $s[0]=~s/C/1/g; $s[0]=~s/T/2/g; $s[0]=~s/G/3/g;
	$x=~s/A/0/g; $x=~s/C/1/g; $x=~s/T/2/g; $x=~s/G/3/g;
	$small=$rs0;
	$small = $rx if($x lt $s[0]);
	print "$small\n";
}

