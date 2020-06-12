#!/usr/bin/perl -w
use warnings;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Basename;
use diagnostics;

# VERSION v1.0

$file='';
$score=0.4;

my $basic_usage =  "usage:\n<filter-sto.pl -in .sto -score (dissimilarity value; optional) -order (reorder sequences by closest neighbor joining; optional) -minimal> requires a .sto file to filter out sequences that have only fraction of the predicted base pairing (Default: 0.4)\n";
if(!@ARGV){print "$basic_usage\n";
           exit;}

GetOptions(
        'in=s' => \$file,
	'score:f' => \$score,
	'order' => \$order_sto,
	'minimal' => \$minimal,
	'open_at_end' => \$run,
	'bed' => \$bed,
        'help|h' => \$help) or die ("$basic_usage\n");

if (defined $help) {
print "$basic_usage\n";
exit;
}


@dotbracket = ("(", ")", "[", "]", "<", ">", "{", "}");
$position="";
$str="";
$cons="";


open(STO, $file);
while($line = <STO>){ 
	if($line=~/# STOCKHOLM 1.0/){next;}
	elsif($line=~/#=GF/){next;}
	elsif($line=~/^\w/){
	($species, $seq) = $line =~ m/(^.+?\s+)(.+)/g;
#	$seq_str{ $species } = $seq;
	$species =~s/\(.\)/   /;
	$seq_str{ $species } .= "" if exists $seq_str{$species};
	$seq_str{$species} .= $seq;
	
}
    elsif($line=~ /#=GC SS_cons/){
	$line =~ s/~/\./g;
	($gcstr, $str) = $line =~ m/(#=GC SS_cons\s+)(.+)/g;
	$gcss{$gcstr} .= "" if exists $gcss{$gcstr};
	$gcss{$gcstr} .= $str;
	$struc = $gcss{$gcstr};
	} elsif($line=~/#=GC RF/) {
	$line =~ s/~/\./g;
	($gccons, $cons) = $line =~ m/(#=GC RF\s+)(.+)/g;
	$gcrf{$gccons} .= "" if exists $gcrf{$gccons};
	$gcrf{$gccons} .= $cons;
	} elsif($line=~/\/\//){next;}
	else { next; }
}
print "Number of sequences in file: ", scalar keys %seq_str, "\n";

foreach $b(@dotbracket){
$offset=0;
$position = index($struc, $b, $offset);
push @positions, $position;
 while ($position != -1) {

    #print "Found $b at $position\n";

    $offset = $position + 1;
    $position = index($str, $b, $offset);
	push @positions, $position;

  }
pop @positions;
}

$nopos = scalar @positions;
$out = basename($file, ".new.sto");
$inbed = "$out\.out.bed";
$outfilter = "$out\.new.filtered_" . "$score\.sto";
$tmp = "tmp";
#$outdeleted = "$out\.deleted_" . "$score\.sto";
open (TMP, '>', $tmp);
open ($fh, '>', $outfilter);
print $fh "# STOCKHOLM 1.0\n#=GF AU Infernal 1.1.1\n";
#open ($dh, '>', $outdeleted);
#print $dh "# STOCKHOLM 1.0\n#=GF AU Infernal 1.1.1\n";

while (($species, $seq) = each(%seq_str)) {
	$count=0;
	foreach $a(@positions){
	    $base="";
		$base = substr $seq, $a, 1; 
	        if ($base =~ m/-/){
			$count++;} 
		else {next;}
}
	print "Found $count -s out of $nopos structured positions in $species\n";
	
	if($count/$nopos <=$score){
	print $fh "$species$seq\n";
	($tag,$from) = $species =~ /(.+):(.+)-.+/;
	print TMP "$tag\t$from\n";
	} else {
#	print $dh "$species$seq\n";
	delete $seq_str{$species};
}
}

foreach $gcstr(keys %gcss){
print $fh "$gcstr$gcss{$gcstr}\n";
#print $dh "$gcstr$gcss{$gcstr}\n";
}
foreach $gccons(keys %gcrf){
print $fh "$gccons$gcrf{$gccons}\n";
#print $dh "$gccons$gcrf{$gccons}\n";
}
print $fh "//\n";
print "Sequences left after deleting ones with only $score of predicted base pairing: ", scalar keys %seq_str, "\n";


close STO;
close $fh;
#close $dh;
#system("sed -i 's/(.)/   /' $outfilter");
close $tmp;

if (defined $order_sto) {
order();
}

if (defined $minimal){
	unlink $outfilter;
	unlink $file;

}

if (defined $run){
	$stofile = $outStoFileName;
	system("emacs $stofile");
}

if (defined $tmp){
	$outbed = $outStoFileName . ".bed";
	system("grep -Ff tmp $inbed >$outbed");
	unlink $tmp;
}


sub order {
$inStoFileName=$outfilter;
$outStoFileName="$outfilter\.ordered.sto";

if (!open(IN,$inStoFileName)) {
    die "cannot open $inStoFileName";
}
if(!open(OUT,">$outStoFileName")){
    die "cannot open $outStoFileName";
}
print OUT "# STOCKHOLM 1.0\n";
#%idToSeq;
$idNum=0;
while (<IN>) {
    s/[\r\n]//g;

    if (/^\#/) {
        if (/^\#=GC ([-_A-Za-z0-9]+)([ \t]+)(.*)$/) {
            $field=$1;
            $spaces=$2;
            $seq=$3;
            if (length($gcFieldToSeq{$field})==0) {
                push @fieldList,$field;
            }
            $gcFieldToSpaces{$field}=$spaces;
            $gcFieldToSeq{$field}=$gcFieldToSeq{$field}.$seq;
        }
        elsif (/^\#=GR/){next;}
        else {
            if (/^\#=G.*$/) {
                print OUT "$_\n";
            }
        }
    }
    else {
        if (/^([^ \t]+)([ \t]+)(.*)$/) {
            # original: /^([-._A-Za-z0-9]+\/[0-9]+[-][0-9]+)([ \t]+)(.*)$/) {
            $id=$1;
            $spaces=$2;
            $seq=$3;
            $idToSpaces{$id}=$spaces;
            $isNew=(length($idToSeq{$id})==0);
            $idToSeq{$id}=$idToSeq{$id}.$seq;
            if ($isNew) {
                $idNumToId{$idNum}=$id;
                $idNum++;
            }
        }
    }
}
$numSeqs=$idNum;

$SS_cons=$gcFieldToSeq{"SS_cons"};
if (!$SS_cons) {
    die "no SS_cons field -- this is not wrong, but I'm not going to deal with it.";
}
$length=length($SS_cons);
$dashes=$SS_cons;
$dashes =~ s/./-/g;

for ($idNum=0; $idNum<$numSeqs; $idNum++) {
    $id=$idNumToId{$idNum};
    $seq=$idToSeq{$id};
    $seq =~ s/-/./g;
    #print $seq;
    @{$idNumToSeqArray{$idNum}}=split //,$seq;
}

$numClusters=$numSeqs;
for ($c=0; $c<$numClusters; $c++) {
    push @{$cluster[$c]},$c;
}

while ($numClusters>1) {

    DumpClusters();

    $bestC1=-1;
    $bestC2=-1;
    $bestD=$length; # worst possible
    for ($c1=0; $c1<$numClusters-1; $c1++) {
        for ($c2=$c1+1; $c2<$numClusters; $c2++) {
            for $s1 (@{$cluster[$c1]}) {
                for $s2 (@{$cluster[$c2]}) {

                    if ($didDist{$s1}{$s2}) {
                        $d=$dist{$s1}{$s2};
                    }
                    else {
                        $d=0;
                        for ($i=0; $i<$length; $i++) {
                            $ch1=uc($idNumToSeqArray{$s1}[$i]);
                            $ch2=uc($idNumToSeqArray{$s2}[$i]);
                            if (!($ch1 eq $ch2)) {
                                $d++;
                            }

                            #print "$ch1,$ch2 ";
                        }
                        #print "$c1,$c2,$d ::: ";
                        $didDist{$s1}{$s2}=1;
                        $dist{$s1}{$s2}=$d;
                    }
                    if ($d<$bestD) {
                        $bestD=$d;
                        $bestC1=$c1;
                        $bestC2=$c2;
                    }
                }
            }
        }
    }
    push @{$cluster[$bestC1]},@{$cluster[$bestC2]};
    @{$cluster[$bestC2]}=@{$cluster[$numClusters-1]};
    $numClusters--;
}
for $idNum (@{$cluster[0]}) {
    $id=$idNumToId{$idNum};
    $seq=$idToSeq{$id};
    $spaces=$idToSpaces{$id};
    print OUT "$id$spaces$seq\n";
    #print "$id$spaces$seq\n";
}
for $field (@fieldList) {
    $spaces=$gcFieldToSpaces{$field};
    $seq=$gcFieldToSeq{$field};
    print OUT "#=GC $field$spaces$seq\n";
}
print OUT "//\n";
close(OUT);
close(IN);

print "Done.\n";


sub DumpClusters {
    print "clusters: ";
    for ($c=0; $c<$numClusters; $c++) {
        print "$c:{";
        for $idNum (@{$cluster[$c]}) {
            print "$idNum,";
        }
        print "}  ";
    }
    print "\n";

}
}

unlink $tmp;
exit;

