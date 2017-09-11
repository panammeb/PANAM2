#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Cwd 'abs_path';
use FindBin;

use lib "$FindBin::RealBin/modules";

use clustering_64;
use parse_ini;

##################################################
my $command; # Commande linux
my $abundance;
my $pourcentage; my $flag_pourcentage=0;
my $nb_total_sequences;
my $fasta_clean; # Fichier à traiter
my $ct; 
my $log_file;
my $script="";

##################################################


my ($USAGE) =  "\n\n\t****** USAGE of $0 PROGRAM ******\n\n\n\n\tUSAGE: perl $0 <panam.ini file> \n\n\n\n";
die "$USAGE" if ((scalar(@ARGV))< 1);
my $option_file = $ARGV[0];
chomp($option_file);

die "\n\n\t=> Cannot find configuration file: $option_file.\n\n" unless (-e $option_file);
die "\n\n\t=> Configuration file, $option_file, appears to be empty!\n\n" if (-z $option_file);


my @parse = &parse_ini($option_file);

my  $path_results = $parse[0]; 
die "\n\n\tOutput directory for NGS analyses is not defined. Check $option_file.\n\n" unless ( defined  $path_results);
if (! -e $path_results){
`mkdir $path_results`;
}

my $ident = $parse[20];
die "\n\n\t Clustering cutoff missed. Check panam.ini.\n\n" unless ($ident <=1);
die "\n\n\t Invalid input! clustering cutoff must be specified as a fractional identity in the range 0.0 to 1.0\n\n" unless (($ident ne "") and ($ident >=0));

############# Infile ############################################################################

if (! -e $path_results."/quality_output/seqAll.fasta") # généré lorsque un dossier est spécifié
{

    $fasta_clean=$parse[33];
    if ( $fasta_clean eq ""){

    $command="cat ".$path_results."/quality_output/seqAll_* > ".$path_results."/quality_output/seqAll.fasta";
    qx($command);
    $fasta_clean=$path_results."/quality_output/seqAll.fasta";
    }
    else{
    die "\n\n\t=> Cannot find file for clustering process: $fasta_clean\n\n" unless (-e $fasta_clean);

    }

}
else
{
 $fasta_clean=$path_results."/quality_output/seqAll.fasta";
}
###################################################################################################

my $filter_sequences = $parse[29];# Abondance ou pourcentage

if ($filter_sequences=~/(\d*)%/){
$pourcentage=$1/100;
$flag_pourcentage=1;
$nb_total_sequences=qx(grep -c ">" $fasta_clean);
}
else{
$abundance=$filter_sequences;
}

#Sauvegarde script 13/04/2016
# $script=$0;
# $command="cp ".$script." ".$path_results."/.".$script;
# qx($command);


# Fichier LOG

$log_file=">>".$path_results."/panam.log";

open (LOG, $log_file);
open (LOG, $log_file);
print "\n\n* STEP 2 - Clustering by PANAM2 \n\n";
print LOG "\n\n* STEP 2 - Clustering by PANAM2\n\n";
print LOG "Clustering with vsearch\nIdentity = $ident\n Results in: $path_results\n\n";
print LOG "File: $fasta_clean\n";
print LOG "Abundance filter: remove clusters lower than $filter_sequences\n";
print "Clustering with vsearch\nIdentity = $ident\n Results in: $path_results\n";
print "File: $fasta_clean\n";
print "Abundance filter: remove clusters lower than $filter_sequences\n";





#############" SCRIPT

qx(mkdir $path_results/preprocess_output);
qx(mkdir $path_results/preprocess_output/pooled_sample);



&vsearch_64($fasta_clean, $path_results, $ident);




######### Génération pooled_sample_OTU pour lancer panam.pl ################

open (F1, $path_results."/preprocess_output/mapping.uc") || die "Don't found ".$path_results."/preprocess_output/mapping.uc";


my $seed;my %hseed;my $seqseed; my %hseqseed;
while (my $l =<F1>){
        chomp ($l);
	my @tab = split ("\t", $l);
        if ($l =~ m/^H/){
               
                $seed = $tab[9];
                $seqseed = $tab[8];
                $hseed{$seed} ++;
                $hseqseed{$seed} .= $seqseed.", ";
	}


}



# Génération du pool_OTU
open (F1, ">$path_results/preprocess_output/pooled_sample/pooled_sample_OTU");
open (BAD, ">$path_results/preprocess_output/pooled_sample/BAD_OTU");
print F1 "OTU\tSeed sequence\t#Sequences\tSequences\n";
print BAD "OTU\tSeed sequence\t#Sequences\tSequences\n";

$ct=1;
my $bad_sequences=0;

foreach my $k (keys(%hseqseed)) { 
# $k=~/(.*)_.*/;# extraction de l'environnement 11/04/2016

################## 4/05/2016

  if($flag_pourcentage==0){

    if($hseed{$k} >= $abundance)
    {
    print F1 "OTU_$ct\t$k\t$hseed{$k}\t$hseqseed{$k}\n"; # 11/04/2016
	$ct++;
    }
    else
    {
    print BAD "$k\t$hseed{$k}\t$hseqseed{$k}\n";
    delete($hseed{$k}); 
    delete($hseqseed{$k});
    $bad_sequences++;
    }
  }

  else{

   if(($hseed{$k}/100) >= $pourcentage)
    {
    print F1 "OTU_$ct\t$k\t$hseed{$k}\t$hseqseed{$k}\n"; # 11/04/2016
	$ct++;
    }
    else
    {
    print BAD "$k\t$hseed{$k}\t$hseqseed{$k}\n";
    delete($hseed{$k}); 
    delete($hseqseed{$k});
    $bad_sequences++;
    }

    }
  

################## 4/05/2016

}

print "Abundance filter:  $bad_sequences removed\n";


open (F1, $path_results."/preprocess_output/Seqs_seed".$ident.".fasta");
my $nomSeq;my %hseq_fw;
while (my $l = <F1>){
        chomp ($l);
        if ($l =~ m/^>(.*)/){
                $nomSeq = $1;
        }
        else {
                $hseq_fw{$nomSeq} .= $l;
        }
}
close F1;

open (OTU, "$path_results/preprocess_output/pooled_sample/pooled_sample_OTU");
open (F1, ">$path_results/preprocess_output/pooled_sample/pooled_sample_SEQ_OTU");

while (my $l =<OTU>){
        chomp ($l);
	  if ($l !~ m/.*Seed sequence.*/)
	  {
	  my @tab = split ("\t", $l);
                print F1 ">$tab[0]\n$hseq_fw{$tab[1]}\n";
	  }

        }

close F1;

#13/04/2016
# $command="rm $fasta_clean";
# qx($command);
$command="rm ".$path_results."/preprocess_output/*.fasta";
qx($command);

