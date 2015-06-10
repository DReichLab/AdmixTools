#!/usr/bin/perl

# $parfile = "par.ANCESTRYMAP.EIGENSTRAT";
  $parfile = "par.EIGENSTRAT.PED"; 
# $parfile = "par.PED.EIGENSTRAT";                    
# $parfile = "par.PED.PACKEDPED";                     
# $parfile = "par.PACKEDPED.PACKEDANCESTRYMAP";       
# $parfile = "par.PACKEDANCESTRYMAP.ANCESTRYMAP";     

system("../bin/convertf -p $parfile");
