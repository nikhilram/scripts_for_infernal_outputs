#!/usr/bin/perl -w
use warnings;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Basename;
use diagnostics;
use Term::ANSIColor 4.00 qw(RESET :constants);
{$Term::ANSIColor::AUTORESET=1;}

# Version v1.0

my $start_time = time();

#########################################################################################################################################################################################$
#
#
#									Absolute path to required programs
#
#
$samtools = "/usr/local/bin/samtools";
$cmsearch = "/usr/local/infernal-1.1.1/bin/cmsearch";
$cmconvert = "/usr/local/infernal-1.1.1/bin/cmconvert";
$cmcalibrate = "/usr/local/infernal-1.1.1/bin/cmcalibrate";
$cmbuild = "/usr/local/infernal-1.1.1/bin/cmbuild";
$cmalign = "/usr/local/infernal-1.1.1/bin/cmalign";
$bedtools = "/usr/bin/bedtools";
#
$db_default ="/xxx/xxx/refseq77/refseq77_bacteria_complete.fna"; ### Path to default database 
#
#
#
#								  Substitute with their absolute paths in your system
#
#########################################################################################################################################################################################$





#########################################################################################################################################################################################$
#
#										 Code below here
#
#########################################################################################################################################################################################$

my $file = '';
my $d = '';
my $cpu = 12;
my $out='';
my $e=1;
my $align;
my $minimal;
my $compile;
my $help;

my $basic_usage = qq^

################################################################################################################################################################################################
#
#
#	  									 ncRNA Pipeline
#
#
################################################################################################################################################################################################
#
# Required:
#
# -in <string>			:Either .cm or .sto file
#
# -prefix <int/string>          :Prefix for output files
#
#
#
##############################################################################
#
# Optional:
#
# -db <string>			:Fasta sequence database to run searches against. (.fna) 
#					Default - /xxx/xxx/refseq77/refseq77_bacteria_complete.fna
#
# -cpu <int>			:Number of cores to run the pipeline on. Default - 12
#
# -e                            :e value for cmsearch. Default is 1. 
#
# -realign			:Realigns sequences resulting from the cmsearch back to the original covariance model to generate a new .sto file
#
# -compile			:Creates and compiles the coordinate info for each hit to screen for overlaps later
#
# -minimal			:Removes intermediary files from the disk
#
# -help				:Lists the basic usage of the pipeline
#
#
# Typical usage:
#
#		pipeline.pl -in <.sto or .cm> -db <fasta file> -prefix <parent directory name> -realign -compile
###############################################################################################################################################################################################
^; 

if(!@ARGV){print BOLD "$basic_usage\n";
	   exit;}

GetOptions(
        'in=s' => \$file,
        'db:s' => \$db,
        'cpu:i' => \$cpu,
        'prefix:s' => \$out,
        'realign' => \$align,
	'e:i' => \$e,
        'compile' => \$compile,
        'minimal' => \$minimal,
        'help|h' => \$help) or die("$basic_usage\n");


if (defined $help) {
print("$basic_usage\n");
exit;
}


$db=$db_default unless defined $db;
$dbfai = $db . '.fai';
open(my $fasta, $db);
if (<$fasta> !~ /^>/){
	print BOLD ("\n\n********************************************************\n\nSubmitted database isn't a fasta file! Please check!\n\n********************************************************\n\n");
	exit;
}
close $fasta;
print BOLD ("\n\n********************************************************\n\nScanning database directory for index..\n\n********************************************************\n\n");
if (-e $dbfai){
	print BOLD ("\n\n********************************************************\n\nIndex found, database is all set, lets proceed!\n\n********************************************************\n\n");
} else {
	print BOLD ("\n\n********************************************************\n\nDatabase file is not indexed, will be indexed first using Samtools\n\n********************************************************\n\n");
	system("$samtools faidx $db");
}



$dbname = basename($db, ".fna");
$filename_input = basename($file, ".cm.gz");
#@suffixes=(".sto", ".sth", ".cm");
#$filename = basename($file, @suffixes);


#$outcm = $filename . '.' . $dbname . '.new.cm';
#$outref = $filename . '.' . $dbname . '.out';
#$outrefcsv = $filename . '.' . $dbname . '.out.csv';
#$outrefbed = $filename . '.' . $dbname . '.out.bed';
#$outreffasta = $filename . '.' . $dbname . '.new.fasta';
#$outrefsto = $filename . '.' . $dbname . '.new.sto';
$model = $file;


if (defined $out) {
	$outcm = $out . '.' . $filename_input . '.' .  $dbname . '.new.cm';
	$outref = $out . '.' . $filename_input . '.' .  $dbname . '.out';
	$outrefcsv = $out . '.' . $filename_input . '.' .  $dbname . '.out.csv';
	$outrefbed = $out . '.' . $filename_input . '.' .  $dbname . '.out.bed';
	$outreffasta = $out . '.' . $filename_input . '.' .  $dbname . '.new.fasta';
	$outrefsto = $out . '.' . $filename_input . '.' .  $dbname . '.new.sto';
	$model = $out;
}

$newcm = "./$outcm";
$refout = "./$outref";
if ($file !~ /sto|cm/){
print BOLD ("\n\n********************************************************\n\nInput file is neither CM nor STO, please check\n\n********************************************************\n\n");
exit;
}

elsif ($file =~ /sto/){
                if (-e $refout){
                mkdir(next_round);
                system("$cmbuild next_round/$outcm $file");
                print BOLD ("\n\n********************************************************\n\nCM file built from sto!\n\n********************************************************\n\n");
                $new_dir = "./next_round/";
                chdir $new_dir or die "Cannot switch into new directory\n";
                system("$cmcalibrate --cpu $cpu $outcm");
                print BOLD ("\n\n********************************************************\n\nCM file calibrated\n\n********************************************************\n\n");
                system("$cmsearch --tblout $outref --cpu $cpu -E $e $outcm $db");
                gen_csv();
                } else {
                system("$cmbuild $outcm $file");
                print BOLD ("\n\n********************************************************\n\nCM file built from sto!\n\n********************************************************\n\n");
                system("$cmcalibrate --cpu $cpu $outcm");
                print BOLD ("\n\n********************************************************\n\nCM file calibrated, will search against: \n$db\n\n********************************************************\n\n");
                system("$cmsearch --tblout $outref --cpu $cpu -E $e $outcm $db");
                gen_csv();
                        }
        }

elsif ($file =~ /cm/){
open(my $fh,"<", $file);
        while($row = <$fh>){
                if($row =~ /INFERNAL/){ 
                $firstline = $row;
                } elsif ($row =~ /cmcalibrate/) {
                $COM = $row;
                } else {next;} 
                }
                close $fh;
if (-e $refout){
gen_csv();
}
elsif (-e $newcm){
	print BOLD ("\n\n********************************************************\n\nNo refseq out file found but calibrated model exists, performing search\n\n********************************************************\n\n");
	system("$cmsearch --tblout $outref --cpu $cpu -E $e $outcm $db");
	gen_csv();
} else {
	if ($firstline =~ /1.1/){
		if ($COM =~ /cmcalibrate/){
			print BOLD ("\n\n********************************************************\n\nCalibrated Infernal 1.1.1 CM file given, will search against: \n$db\n\n********************************************************\n\n");
        		system("$cmsearch --tblout $outref --cpu $cpu -E $e $file $db");
        		gen_csv();
			} else {
				print BOLD ("\n\n********************************************************\n\nCM file will be calibrated!\n\n********************************************************\n\n");
			        system("$cmcalibrate --cpu $cpu $outcm");
			        print BOLD ("\n\n********************************************************\n\nCM file calibrated, will search against: \n$db\n\n********************************************************\n\n");
			        system("$cmsearch --tblout $outref --cpu $cpu -E $e $outcm $db");
        			gen_csv();
				}
			} elsif ($firstline !~ /1.1/) {	
					print BOLD ("\n\n********************************************************\n\nCM file submitted not of Infernal 1.1.1 , will be converted to 1.1.1!\n\n********************************************************\n\n");
					system("$cmconvert -o $outcm $file"); 
					print BOLD ("\n\n********************************************************\n\nCM file converted, will be calibrated!\n\n********************************************************\n\n");
					system("$cmcalibrate --cpu $cpu $outcm");
					print BOLD ("\n\n********************************************************\n\nCM file calibrated, will search against: \n$db\n\n********************************************************\n\n");
					system("$cmsearch --tblout $outref --cpu $cpu -E $e $outcm $db");	
					gen_csv();
		} else {print BOLD "\n\n********************************************************\n\nUnknown CM file origin. Requires CM file from either Infernal 1.0.2 or 1.1.1\n\n********************************************************\n\n";}

	}
	} 
#elsif ($file =~ /sto/){
#		if (-e $refout){
#		mkdir(next_round);
#		system("$cmbuild next_round/$outcm $file");
 #       	print BOLD ("\n\n********************************************************\n\nCM file built from sto!\n\n********************************************************\n\n");
  #  	        $new_dir = "./next_round/";
#		chdir $new_dir or die "Cannot switch into new directory\n";
#		system("$cmcalibrate --cpu $cpu $outcm");
 #       	print BOLD ("\n\n********************************************************\n\nCM file calibrated\n\n********************************************************\n\n");
  #      	system("$cmsearch --tblout $outref --cpu $cpu -E $e $outcm $db");
   #     	gen_csv();
#		} else {
#		system("$cmbuild $outcm $file");
 #       	print BOLD ("\n\n********************************************************\n\nCM file built from sto!\n\n********************************************************\n\n");
  #      	system("$cmcalibrate --cpu $cpu $outcm");
   #     	print BOLD ("\n\n********************************************************\n\nCM file calibrated, will search against: \n$db\n\n********************************************************\n\n");
    #    	system("$cmsearch --tblout $outref --cpu $cpu -E $e $outcm $db");
     #   	gen_csv();
#			}
#	}

if (defined $align) {
	system("$bedtools getfasta -fi $db -bed $outrefbed -fo $outreffasta");
	print BOLD ("\n\n*****************************************\n\nFasta file generated\n\n*****************************************\n\n");
	if (-e $outcm){
	system("$cmalign -o $outrefsto --cpu $cpu $outcm $outreffasta");
	} else { 
	        system("$cmalign -o $outrefsto --cpu $cpu $file $outreffasta");
		}
	}


if (defined $compile) {
	if (-e $outrefbed){
	system("cut -f1,2,3,4,5 $outrefbed >>../compiled.bed");
	} else {
	print BOLD ("\n\n********************************************************\n\nBed file from run not found in the directory, run pipeline first.\n\n********************************************************\n\n");
		}
	}


if (defined $minimal) { 
	unlink $outreffasta;
	unlink $outrefbed;
	}

sub gen_csv {
$cmsearch_out = $outref;
open(CMs, $cmsearch_out);
open($out, '>', $outrefcsv);
open($bed, '>', $outrefbed);
print $out "gaisr name,orig-start,orig-end,hit-start,hit-end,hit-strand,hit-rel-start,hit-rel-end,score,evalue,pvalue,gc,structure,tgt-seq,/cluster/home/nmr07003/metatranscriptomes\n";
while($line = <CMs>){
	chomp $line;
    if($line !~ m/^#/){
	my @chunks;
	my $org="";
	@chunks = split /\s+/, $line;
	$chunks[0] =~ s/\.\d//;
        if($chunks[0]=~ /N/){
                if($chunks[9]=~/-/){
	$org = $chunks[17] . " " .  $chunks[18];
        print $out "$chunks[0],$chunks[5],$chunks[6],$chunks[7],$chunks[8],-1,$chunks[7],$chunks[8],$chunks[14],$chunks[15],0,$model,$org\n";
        print $bed "$chunks[0]\t$chunks[8]\t$chunks[7]\t$chunks[2]\t$chunks[15]\t$chunks[9]\t$model\n" if ($chunks[15] <= 0.0001);
	} else {
                 print $out "$chunks[0],$chunks[5],$chunks[6],$chunks[7],$chunks[8],1,$chunks[7],$chunks[8],$chunks[14],$chunks[15],0,$model,$org\n";
		 print $bed "$chunks[0]\t$chunks[7]\t$chunks[8]\t$chunks[2]\t$chunks[15]\t$chunks[9]\t$model\n" if ($chunks[15] <= 0.0001);
		}} else { next; }

    }
else {next;}
}
close $out;
close $bed;
close CMs;

}
system("sort -k1,1 -k2,2n $outrefbed -o $outrefbed");
my $end_time = time();
my $run_time = $end_time - $start_time;
my $run_time_min = $run_time/60;
print BOLD ("\n\n********************************************************\n\nPipeline took $run_time_min minutes to complete running on $cpu cores\n\n********************************************************\n\n");
exit;

