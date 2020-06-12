#!/usr/bin/perl
use warnings;


$file = $ARGV[0]; # Input CMtblout format hits from cmsearch
$out = $ARGV[1]; # Output filename
$e = $ARGV[2] || 0.0001; # E-value cutoff - default of 0.0001

if (not defined $file){die "usage: <bed_from_cmtblout.pl>\n bed_from_cmtblout.pl input cmsearch-tblout output.bed\n"};

open(OUT,">$out");
open(CMs, $file);
while($line = <CMs>){
        chomp $line;
    if($line =~ m/^N/){
        my @chunks;
        @chunks = split /\s+/, $line;
        if($chunks[0]=~ /N/){
		$len=abs($chunks[8]-$chunks[7]);
		$chunks[0] =~ s/\.\d//g;
	   if($chunks[15] < $e){	
                if($chunks[9]=~/-/){
        print OUT "$chunks[0]\t$chunks[8]\t$chunks[7]\t$len\t$chunks[15]\t$chunks[9]\t$chunks[2]\n";
        } else {
                 print OUT "$chunks[0]\t$chunks[7]\t$chunks[8]\t$len\t$chunks[15]\t$chunks[9]\t$chunks[2]\n";
                }} else { next; }
	}else {next;}
	
    }
else {next;}
}
close CMs;

system("sort -k1,1 -k2,2n $out -o $out");
exit;


