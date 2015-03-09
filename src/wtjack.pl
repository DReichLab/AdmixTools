#!/usr/bin/perl
use Getopt::Std ;

getopts('i:m:b:',\%opts) ;

$BIN ="/groups/reich/pm82/dev/ADMIXTOOLS/bin";
 
# filename of rolloff output file, in parfile output:
if (defined $opts{"i"}) {
        $pop = $opts{"i"} ;
}
else {
        die "input parameter compulsory\n"; 
}

# SNP file
if (defined $opts{"m"}) {
   	$map = $opts{"m"} ;
}
else {
    	die "m parameter compulsory- Enter map file name\n" ;
}

if (defined $opts{"b"}) {
        $bin = $opts{"b"} ;
}
else {
        die "path to bin directory required\n";
}


open (IN, "$map") or die "COF mapfile\n";
@count = ();
$count_total =0;

#count the number of SNPs on each chromosome. This is used as weights for the weighted jackknife calculation
while (my $line = <IN>)
{
	chomp $line; 
	@data =split(" ",$line);
	for (my $i = 1; $i <= 22; $i++)
	{
		if ($data[1] == $i)
		{
			$count[$i]++;
			$count_total++;
		}
	}
}

# open expfit output
$log = "expfit_${pop}.log";
open (LOG, $log) or die "cannot open $log\n";
@log = <LOG>;

$flog = "expfit_${pop}.flog";
open (FLOG, $flog) or die "cannot open $flog\n";

# open output file
open (JIN, ">${pop}.jin") or die "COF outfile\n";

$chr  = 1;
for (my $i = 0; $i < scalar(@log); $i++)
{
        $line = $log[$i];
        my $string = "${pop}:${chr}";
        unless ($line =~ /$string/) {next;}
        $start[$chr] = $i;
        $chr++;
}
push(@start, scalar(@log));

printf JIN "%3s  %12s  %12s\n",   #CHR,    SNPs,   Estimated_date;
for my $chr(1..22)
{
        for (my $i = $start[$chr]; $i <  $start[$chr+1]; $i++)
        {
                my $line = $log[$i];
                if ($line =~ /mean/) {
                        chomp $line;
                        @array = split(/[\t ]+/,$line);
                        $this = $array[2];
                        printf JIN "%3d  %12.6f  %12.6f\n",   $chr,    $count[$chr],   $this;

                }
        }
}

# Find mean for the full run
while (my $line = <FLOG>)
{
	if ($line =~ /mean/) { 
		chomp $line;
        	@array = split(/[\t ]+/,$line);
        	$mean = $array[2];
	}
}
#print "mean: $mean\n";
# Run jackknife
$jin = $pop.".jin";
$jout = $pop.".jout";
print "\nJackknife summary: ${pop}.jin\n";
system "$bin/dowtjack -i $jin -m $mean -o $jout";
$ans = `tail $jout` ;  
($jmean, $serr) = split " ", $ans ;

printf "Jackknife mean:      %9.3f\n", $jmean ;
printf "Jackknife std. err:  %9.3f", $serr ;
print "\n" ;
