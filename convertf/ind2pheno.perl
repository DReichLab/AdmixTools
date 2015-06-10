#!/usr/bin/perl

$in = $ARGV[0]; # .ind file
$out = $ARGV[1]; # .pheno file

open(IN,$in) || die("COF");
open(OUT,">$out") || die("COF");

while($line = <IN>)
{
  if($line =~ /Case/) { print OUT ("1"); $case=1; }
  elsif($line =~ /Control/) { print OUT ("0"); $control=1; }
  else { print OUT ("9"); }
}
print OUT ("\n");
unless($case) { print("WARNING: no cases\n"); }
unless($control) { print("WARNING: no controls\n"); }

