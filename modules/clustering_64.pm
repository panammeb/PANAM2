#!/usr/bin/perl -w
use strict;
use warnings;

#############################
# Le 13/04/2016 DD
#
#############################

use Cwd 'abs_path'; # permet d'extraire le path de l'exécutable

sub vsearch_64 
{
my $fasta_clean;
my $path_results;
my $ident;
my $command;
my $uclust;
my $path;

($fasta_clean, $path_results, $ident)=@_;



$path=abs_path($0);
$path =~ /(.*)\/preprocess.*/; # déterminer le path du bin
$uclust="$1/bin/vsearch";



# print "Path $uclust\n";

my $log_file=">>".$path_results."/panam.log";
open (LOG, $log_file);



# Clustérisation SMALLMEM
$command=$uclust." -derep_fulllength ".$fasta_clean." -output ".$path_results."/preprocess_output/Seqs_derep.fasta";
print LOG $command."\n";
qx($command);

$command=$uclust." -sortbylength ".$path_results."/preprocess_output/Seqs_derep.fasta -output ".$path_results."/preprocess_output/Seqs_sorted.fasta";
print LOG $command."\n";
qx($command);

$command=$uclust."  -cluster_smallmem ".$path_results."/preprocess_output/Seqs_sorted.fasta -id ".$ident." -centroids ".$path_results."/preprocess_output/Seqs_seed".$ident.".fasta -uc ".$path_results."/preprocess_output/mapident.uc ";
print LOG $command."\n";
qx($command);

$command=$uclust." -usearch_global ".$fasta_clean." -db ".$path_results."/preprocess_output/Seqs_seed".$ident.".fasta -id ".$ident." -strand plus -uc ".$path_results."/preprocess_output/mapping.uc";
print LOG $command."\n";
qx($command);


close LOG;

}
1;