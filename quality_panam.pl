#!/usr/bin/perl
use warnings;
use strict;
use Cwd 'abs_path';
use FindBin;

use lib "$FindBin::RealBin/modules"; 
use quality;
use parse_ini;

###########################
my $command1="" ; # commande linux
my $command2="" ; # commande linux
my $nb_sequences; # Nb de séquances après nettoyage
my $sample; # Echantillons à demultiplexer
my $log_file; #fichier panam.log en écriture
my @barIN; # Tags
my @ligne; #tableau générique
my $demul; #booleen demultiplexage ou non
my $demul_folder;
my $first_option_file;
##########################


my ($USAGE) =  "\n\n\t****** USAGE of $0 PROGRAM ******\n\n\n\n\tUSAGE: perl $0 <panam.ini file> \n\n\n\n";
die "$USAGE" if ((scalar(@ARGV))< 1);
my $option_file = $ARGV[0];
chomp($option_file);

die "\n\n\t=> Cannot find configuration file: $option_file.\n\n" unless (-e $option_file);
die "\n\n\t=> Configuration file, $option_file, appears to be empty!\n\n" if (-z $option_file);

my $path=abs_path($0);
$path =~ /(.*)\/quality.*/; # déterminer le path du bin
my $path_panam="$FindBin::RealBin/"; # 


my @parse = &parse_ini($option_file);

my $NGS_id_Results = $parse[0];
die "\n\n\tOutput directory for NGS analyses is not defined. Check $option_file.\n\n" unless ( defined $NGS_id_Results);
unless (-d $NGS_id_Results) { mkdir $NGS_id_Results || die "Could not create Directory $NGS_id_Results !\n"; }




$demul_folder=$parse[34];
$demul=$parse[35];

if ($demul_folder eq "")
{
    if ($demul eq "Y")
    {
    print "\t\t- Processing script quality for demultiplexing files:\n";
    }
    else
    {
    print "\n\n*  STEP 1 - Processing script quality PANAM2 for Illumina sequencing:\n";
    }
quality_total();
}
else
{
print "\n\n* STEP 1 - Processing script quality without demultiplexing PANAM2 for Illumina sequencing:\n\n";

# Stockage des valeurs initiales du INI
$first_option_file="cp ".$option_file." ".$NGS_id_Results."/.first.ini";
system($first_option_file);
quality_withoutdemultiplexing();
}







sub quality_total
{

if (-e "$NGS_id_Results/quality_output/") {qx(rm -R "$NGS_id_Results/quality_output/");}# Sinon ajout de séquences à l'existant
system("mkdir $NGS_id_Results/quality_output/");
system("mkdir $NGS_id_Results/quality_output/BAD");
system("mkdir $NGS_id_Results/quality_output/tmp");


$log_file=">".$NGS_id_Results."/panam.log";# Initialisation du fichier/début d'analyse

open (LOG, $log_file);
print LOG "Script Quality PANAM2\n";

print LOG "PANAM2 version: ";
system("cat $path_panam/.version >> $NGS_id_Results/panam.log");
print LOG "\n";

close LOG;

# my $dataType = $parse[1];
# die "\n\n\t Data Type value missed. Check panam.ini.\n\n" unless ($dataType ne "");
# die "\n\n\t Data Type value is not correct. Check $option_file.\n\n" unless (($dataType eq "454Roche") or ($dataType eq "IlluminaMiSeq"));
my $dataType="IlluminaMiSeq";

my $MinSeqLength = $parse[4];
die "\n\n\t Min Length cutoff missed. Check $option_file.\n\n" unless ($MinSeqLength ne "");

my $MaxSeqLength = $parse[5];
die "\n\n\t Max Length cutoff missed. Check $option_file.\n\n" unless ($MaxSeqLength ne "");

# my $Nbase = $parse[8];
# die "\n\n\t Ambiguous base value missed. Check $option_file.\n\n" unless ($Nbase ne "");
# die "\n\n\t Ambiguous base value should be yes or no. Check $option_file.\n\n" unless (($Nbase eq "yes") or ($Nbase eq "no")) ;
my $Nbase="yes";

#Fichiers R1 R2
my $inputSeqNameF = $parse[2];

die "\n\n\tCannot find sequence file. Check $option_file.\n\n" unless (($inputSeqNameF ne "") and (-e $inputSeqNameF));
die "\n\n\tSequence file, $inputSeqNameF, appears to be empty!\n\n" if (-z $inputSeqNameF);
chomp ($inputSeqNameF);

my $inputSeqNameR = $parse[10];

die "\n\n\tCannot find sequence file. Check $option_file.\n\n" unless (($inputSeqNameR ne "") and (-e $inputSeqNameR));
die "\n\n\tSequence file, $inputSeqNameR, appears to be empty!\n\n" if (-z $inputSeqNameR);
chomp ($inputSeqNameR);

my $kept_primers=$parse[31];

my $forward = $parse[15];
my $reverse = $parse[16];


####### barcode file
my $barIn= $parse[13];

my $barcode;
my $nbrBar;
my $ficBar;

if (!(defined $barIn ) or ($barIn eq "")) {
	$barcode =0;
	$nbrBar = 0;
	$ficBar = "no barcode file" ;
	$barIn="NULL";
}

elsif ((defined $barIn) and ($barIn ne "")){
	$barcode = 1;
	$ficBar = $barIn;

	# See if the barcode file exist and contains something
	die "\n\n\tCannot find barcode file: $barIn. Check $option_file.\n\n" unless (-e $barIn);
	die "\n\n\tBarcode file, $barIn, appears to be empty!\n\n" if (-z $barIn);

	open (BB, $barIn);
	@barIN= <BB>;
	close BB;

	foreach my $e (@barIN) {
		if ($e !~ /\t/) {die "\n\n\tYour bar code file $barIn does not seem to be in the right format. Check $option_file.\n\n";}
		my @ty = split(/\t/, $e);
		if (($ty[0] =~ /^\d/) or ($ty[0]=~ /_/) or ($ty[0]=~ /\//)) { die "\n\n\tThe bar code Ids in $barIn should start with no number and contain no underscorses. Check $option_file.\n\n";}
		if ($ty[1] !~ /[ATCGU]/i) {die "\n\n\tYour bar code file $barIn does not seem to be in the right format. Check $option_file.\n\n";}
		my $nbr = `wc -l $barIn`;
		if ($nbr =~ /(\d.?)\s/) { $nbrBar = $1} 
	}
}


#### input files and parameters

my $MinOverlap = $parse[6]; 
my $MismatchOverlap = $parse[7]; 
my $merge = $parse[11];
my $score = $parse[12];

##### ajout 21/8/2014 $CheckChimeras et $CheckHomopol
my $CheckChimeras = $parse[27];
die "\n\n\t Check chimeras value missed. Check panam.ini.\n\n" unless (($CheckChimeras ne ""));
die "\n\n\t Invalid check chimeras value!\n\n" unless (($CheckChimeras eq "yes") or ($CheckChimeras eq "no"));	

# my $CheckHomopol = $parse[28]; # utile pour MISEQ ???
# $CheckHomopol = "no";
# die "\n\n\t Check homopolymers value missed. Check panam.ini.\n\n" unless (($CheckHomopol ne ""));
# die "\n\n\t Invalid check homopolymers value!\n\n" unless (($CheckHomopol eq "yes") or ($CheckHomopol eq "no"));

#### 20/05/2015 le merger : pandaseq, usearch, pear
my $MergePairsCom=$parse[30];

# print "$NGS_id_Results $parse[13] $forward $reverse  $MinSeqLength $MaxSeqLength $CheckChimeras\n";

my $quality_output;

if ($inputSeqNameR !~ /.*.fastq.\w+/){ # 11/05/2016 Fichier non compressé
my $ligneF1;my $ligneF2;
my $ligneR1;my $ligneR2;

# inputSeqNameF format and integrity # 11/05/2016 suprresion de ces tests pour fichier compressé
# 	my $ligneF1= `sed -n '1p' $inputSeqNameF`; #7/9/2016 Inutile
# 	my $ligneF2= `sed -n '2p' $inputSeqNameF`;
	if (($ligneF1 !~ /^@/) or ($ligneF2 !~ /^[ATCGUNatcgun]/)) {
		die "\n\n\tYour forward fastq file $inputSeqNameF does not seem to contain sequences in fastq format.\n\n";
	}

	## inputSeqNameR format and integrity
	die "\n\n\tCannot find reverse sequence file. Check $option_file.\n\n" unless (($inputSeqNameR ne "") and (-e $inputSeqNameR));
	die "\n\n\tSequence file, $inputSeqNameF, appears to be empty!\n\n" if (-z $inputSeqNameF);
	chomp ($inputSeqNameR);
# 	my $ligneR1= `sed -n '1p' $inputSeqNameR`;
# 	my $ligneR2= `sed -n '2p' $inputSeqNameR`;
	if (($ligneR1 !~ /^@/) or ($ligneR2 !~ /^[ATCGUNatcgun]/)) {
		die "\n\n\tYour reverse fastq file $inputSeqNameR does not seem to contain sequences in fastq format.\n\n";
	}



  $command1=qx(grep -c "^+$"  $inputSeqNameF); # Pas de sens si les fichiers sont compressés
  $command2=qx(grep -c "^+$"  $inputSeqNameR);


  if ($command1 >$command2){
  print "Sequences to merge: $command2\n";
  }
  else
  {
  print "Sequences to merge: $command1\n";
  }



}

if ($dataType eq "IlluminaMiSeq") {

# 
# my $path=abs_path($0);
# $path =~ /(.*)\/quality.*/; # déterminer le path du bin
my $path_bin=$path_panam."/bin/"; # uvsearch serainstalla dans le bin (les autres mergers devront être accessibles)


# 	
# 	
# 	## parameters for merging paired-end reads
# 	die "\n\n\t Invalid value for the minimum overlap between forward and reverse reads. Check $option_file.\n\n" unless (($MinOverlap ne "") and ($MinOverlap > 0));
# 	die "\n\n\t Merging option missed. Check panam.ini.\n\n" unless ($merge ne "");
# 	die "\n\n\t Invalid value for merging option. Check panam.ini.\n\n" unless (($merge eq "fix_sequences") or ($merge eq "perfect_match"));



	if($MergePairsCom eq "pandaseq") # 20/05/2015 DD
	{

            if(`which pandaseq | wc -l` > 0) # check the installation of the software
            {

            if ($merge eq "perfect_match") {
		die "\n\n\t Invalid value for the mismatch on the overlap between forward and reverse reads. Check $option_file.\n\n" unless (($MismatchOverlap ne "") and ($MismatchOverlap >=0)) ;
            }
	
            if ($merge eq "fix_sequences") {
		die "\n\n\t Invalid value for the score that a sequence must meet to be kept in the output. Thse score must be between 0 and 1. Check $option_file.\n\n" unless (($score >0) and ($score < 1)) ;
            }


            print "$inputSeqNameF\n$inputSeqNameR\n$MinSeqLength\n$MaxSeqLength\n$Nbase\n$MinOverlap\n$MismatchOverlap\n$merge\n$score\n";
            my @merge_results = &merge_pairs_pandaseq($inputSeqNameF, $inputSeqNameR, $MinSeqLength, $MaxSeqLength, $Nbase, $MinOverlap, $MismatchOverlap, $NGS_id_Results, $merge, $score);
            $quality_output = $merge_results[0];
            }
            else { print "\n Pandaseq software is not installed, paman process aborted ... \n\n"; die;}  
        }
        

#VSEARCH
	 elsif($MergePairsCom eq "usearch")# 20/05/2015 DD
	 {

	  my @merge_results = &merge_pairs_vsearch($inputSeqNameF, $inputSeqNameR, $MinSeqLength, $MaxSeqLength, $Nbase, $MinOverlap, $MismatchOverlap, $NGS_id_Results, $merge, $score, $path_bin);
	  $quality_output = $merge_results[0];
	 }

#PEAR

	  elsif($MergePairsCom eq "pear")# 21/05/2015 DD
	 {
            if(`which pear | wc -l` > 0) # check the installation of the software
            {
            my @merge_results = &merge_pairs_pear($inputSeqNameF, $inputSeqNameR, $MinSeqLength, $MaxSeqLength, $Nbase, $MinOverlap, $MismatchOverlap, $NGS_id_Results, $merge, $score);
            $quality_output = $merge_results[0];
            }
            else { print "\n Pear software is not installed, panam process aborted ... \n\n"; die;}  
	 }


$command1=qx(grep -c ">"  $quality_output/merged.fasta);chomp($command1);
print "Merging generated $command1 pairs\n";


#DEMULTIPLEXING

demultiplex_miseq_exact ($NGS_id_Results, $barIn, $forward,$reverse, $MinSeqLength, $MaxSeqLength, $CheckChimeras, $kept_primers, $path_bin); 
	
# Séquences par échantillon

open (LOG, ">>".$NGS_id_Results."/panam.log");
$command1="$NGS_id_Results/quality_output/";
print LOG "\n";
print LOG "\t\t- Quality output statistics:\n:\n";

print "\n";
print "\t\t- Quality output statistics:\n";

  foreach my $sample (@barIN){
  chomp($sample);
  @ligne=split("\t", $sample);
    my  $fichier=$NGS_id_Results."/quality_output/seqAll_".$ligne[0].".fasta";
    if(-e $fichier) {
    $nb_sequences=qx(grep -c ">"  $fichier);
    chomp($nb_sequences);
    }
    else{
    $nb_sequences=0;
  }
  print LOG "Sample $ligne[0]: $nb_sequences clean sequences\n";
  print "Sample $ligne[0]: $nb_sequences clean sequences\n";
  }

my $bad_directory="$NGS_id_Results/quality_output/BAD";
print LOG "\nSequences removed:\n"; print "\nSequences removed:\n";
    if (-e $bad_directory."/seqAll_bad_ATGC.fasta") {
    $nb_sequences=qx(wc -l   $bad_directory/seqAll_bad_ATGC.fasta);chomp($nb_sequences);
    print LOG "Sequences with Ns: $nb_sequences\n"; print "Sequences with Ns: $nb_sequences\n";
    }
    if (-e $bad_directory."/seqAll_bad_length.fasta") {
    $nb_sequences=qx(wc -l  $bad_directory/seqAll_bad_length.fasta);chomp($nb_sequences);
    print LOG "Sequences shorter or longer: $nb_sequences\n"; print "Sequences shorter or longer: $nb_sequences\n";
    }
    if (-e $bad_directory."/seqAll_bad_primers.fasta") {
    $nb_sequences=qx(wc -l   $bad_directory/seqAll_bad_primers.fasta);chomp($nb_sequences);
    print LOG "Sequences without primers (not mismatch allowed): $nb_sequences\n";print "Sequences without primers (not mismatch allowed): $nb_sequences\n";
    }
    if (-e $bad_directory."/seqAll_withoutTags.fasta") {
    $nb_sequences=qx(wc -l   $bad_directory/seqAll_withoutTags.fasta);chomp($nb_sequences);
    print LOG "Sequences without tags (not mismatch allowed): $nb_sequences\n";print "Sequences without tags (not mismatch allowed): $nb_sequences\n";
    }

close LOG;
}

}


sub quality_withoutdemultiplexing
{

if (-e "$NGS_id_Results/quality_output/") {qx(rm -R "$NGS_id_Results/quality_output/");}# Sinon ajout de séquences à l'existant
system("mkdir $NGS_id_Results/quality_output/");



$log_file=">".$NGS_id_Results."/panam.log";# Initialisation du fichier/début d'analyse

open (LOG, $log_file);
print LOG "Script Quality for demultiplexed files PANAM2\n";
open (TAG, ">".$NGS_id_Results."/.barcode.tag");

my @demultiplexed_files=<"$demul_folder/*_R1*.fastq*">;


foreach my $R1 (@demultiplexed_files)
{
    if ($R1 =~ /(^[^_]*)(_.*)_R1.*fastq(.*)/)
    {
    my $path_demul=$1;
    $path_demul=~ /$demul_folder\/(.*)/;
    my $tag_demul=$1;
    print "\t- Sample: $tag_demul\n";
    print TAG $tag_demul."\n";
    
    my $R2=$R1;
    $R2 =~ s/_R1/_R2/g;
    print "Quality step for the folowing files :\n$R1\n$R2\n";
    
        my $path_new_results= $NGS_id_Results."/quality_output/".$tag_demul;
        system("mkdir $path_new_results");
        # Fichier INI
        open (INI,  $NGS_id_Results."/.first.ini");
        open (NEW, ">".$path_new_results."/".$tag_demul.".ini");
        while( my $li=<INI>)
        {
            if ( $li =~ /PATH_RESULTS/)
            {
            print NEW "PATH_RESULTS\t$path_new_results\n";
            }
            elsif($li =~ /INPUT_FILE_FORWARD/)
            {
            print NEW "INPUT_FILE_FORWARD\t$R1\n";
            }
            elsif($li =~ /INPUT_FILE_REVERSE/)
            {
            print NEW "INPUT_FILE_REVERSE\t$R2\n";
            }
            elsif($li =~ /DEMUL_FOLDER/)
            {
            print NEW "\n";
            }
            else
            {
            print NEW $li;
            }
        }
    print NEW "DEMUL\tY\n"; # variable tmp pour affichage des étapes
    close NEW;
    close INI;
    $command1="perl ".$path_panam."/quality_panam.pl ".$path_new_results."/".$tag_demul.".ini";
    print LOG $command1."\n";
    system($command1);
    $command1="sed s/sample/".$tag_demul."/g -i ".$path_new_results."/quality_output/seqAll_sample.fasta";
    print LOG  $command1."\n";
    system($command1);
    $command1="cat ".$path_new_results."/quality_output/seqAll_sample.fasta >> ".$NGS_id_Results."/quality_output/seqAll.fasta";
    print LOG  $command1."\n";
    system($command1);
    }
}


}







