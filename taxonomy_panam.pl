#!/usr/bin/perl


use Cwd 'abs_path';
use Cwd;
use FindBin;

use lib "$FindBin::RealBin/modules";
use lib "$FindBin::RealBin/bin/bioperl-1.5.2_102";

use phylogeny ;
use Bio::TreeIO;

## VERSION 3.03 ## le fichier PANAM_Affiliation est écrasé en cas de relance / les fichiers PANAM_Affliation_<group> sont effacés / ce fichiers ouverts sont fermés
## VERSION 3.02 ## ajout d'une fonction de tri des barcodes ## JCC avril 2013 ##
## VERSION 3.01 ## construction des taxonomic_distribution par échantillons/méthodes/normalisation ## JCC février 2013 ##

my ($USAGE) =  "\n\n\t****** USAGE of panam.pl PROGRAM ******\n\n\n\n\tUSAGE: perl $0 <panam.ini file> \n\n\n\n";
die "$USAGE" if ((scalar(@ARGV))< 1);
my $option_file = $ARGV[0];
chomp($option_file);

die "\n\n\t=> Cannot find configuration file: $option_file.\n\n" unless (-e $option_file);
die "\n\n\t=> Configuration file, $option_file, appears to be empty!\n\n" if (-z $option_file);
open(OFILE, "<$option_file") || die "Cannot open input file : $option_file\n";


############################################################27/04/2016

my $path_abs=abs_path($0);
$path_abs =~ /(.*)\/taxonomy_panam.*/; # déterminer le path du bin
my $path_panam="$FindBin::RealBin/";

my $usearch;
my $NGS_id_Results;
my $panam_output = "panam_output";
# my $parsing;
my $query_seq; 
my $user_file;
my $preprocess_output = "preprocess_output";
my $seq_F;
my $seq_R;
my $primF;
my $primR;
my $trim =0;
# my $path = "/usr/local/panam/panam_v3/bin/uclust3.0.617";
# my $us_version = "3.0.617";
my $path = $path_panam."bin/vsearch";
my $us_version ="";
####
my $fast = $path_panam."bin/FastTree";
my $align = $path_panam."bin/hmmer-2.3.2/src/hmmalign";
#my $Reference = "/home/panam/Reference" ;
my $Reference;
my $build = $path_panam."bin/hmmer-2.3.2/src/hmmbuild";

############################################################27/04/2016
####
my $errorlog; # Fichier error log
my $log_file;
my %concat;
my $dom; my $d;
my %dom; # domaines
my $aff_norm; 
## VERSION 3.02 ##
my $nbrBar; my $ficBar;
my $barIn; my $barcode;
#### V3.02

my @best;
my %best_lca_similarity; # Taxo lca sur identités identiques
my $path_tmpdiv; # Dossier résultats indice de diversité
my @clean; # Retour du module alignement

################################################ module parse_ini DD 13/04/2016 ##################


use parse_ini ;
my @parse = &parse_ini($option_file);
# 
#Run Identifier
$NGS_id_Results = $parse[0];
if (-d $NGS_id_Results."/".$panam_output) { `rm -r $NGS_id_Results/$panam_output` } 
if (!(-d $NGS_id_Results."/".$panam_output)) 
    { 
    mkdir $NGS_id_Results."/".$panam_output || die "\n\tCould not create Directory $NGS_id_Results/$panam_output $!\n"; 
    mkdir $NGS_id_Results."/".$panam_output."/"."tmp"; # Dossier temporaire des calculs de diversité
    $path_tmpdiv=$NGS_id_Results."/".$panam_output."/"."tmp";
    }

$errorlog= $NGS_id_Results."/.error.log";
$log_file=">>".$NGS_id_Results."/panam.log";

open (LOG, $log_file);
print "\n\n* STEP 3 - Phylogenetic affiliation by PANAM2\n\n";
print LOG "\n\n* STEP 3 - Phylogenetic affiliation by PANAM2\n\n";

#Séquences clustérisées 
$query_seq=$NGS_id_Results."/preprocess_output/pooled_sample/pooled_sample_SEQ_OTU"; # PATH fixe
$user_file = "pooledSamples";

#Primers
$seq_F = $parse[15];
$seq_R = $parse[16];

$primF = $parse[15];
$primR = $parse[16];
$aff_norm = $parse[21]; # Nombre pour la normalisation

####### Domaines/Groupes

print "Domains studied:";
foreach $d (keys %{$parse[32]}){
$dom{$d}=1;
print " $d ";
}
print "\n";

foreach my $d  (keys %dom) {
die "\n\n\t Domain name is not correct. Check panam.ini.\n\n" unless (($d eq "bacteria") or ($d eq "eukaryota") or ($d eq "archaea") or ($d eq "enterovirus"));
}

####### Base de reference

$Reference=$parse[26];
if (! defined ($Reference)){

  if (exists ($dom{'enterovirus'})){
  $Reference=$path_panam."/bd_entero";
  }
  else {

  $Reference=$path_panam."/bd_ssrna";
  }

}
print "Reference database: $Reference\n";

####### barcode file
$barIn= $parse[13];

if ($parse[34] ne ""){ $barIn=$NGS_id_Results."/.barcode.tag" ; } # Données démultiplexées



if (!(defined $barIn ) or ($barIn eq "")) {
	$barcode =0;
	$nbrBar = 0;
	$ficBar = "no barcode file" ;

}
elsif ((defined $barIn) and ($barIn ne "")){
	$barcode = 1;
	$ficBar = $barIn;

	# See if the barcode file exist and contains something
	die "\n\n\tCannot find barcode file: $barIn. Check $option_file.\n\n" unless (-e $barIn);
	die "\n\n\tBarcode file, $barIn, appears to be empty!\n\n" if (-z $barIn);

# 	open (BB, $barIn); # Effacé le 26/06/2017 -> uniquemen les labes sont nécessaires
# 	my @barIN= <BB>;
# 	close BB;
# 
# 	foreach my $e (@barIN) {
# 		if ($e !~ /\t/) {die "\n\n\tYour bar code file $barIn does not seem to be in the right format. Check $option_file.\n\n";}
# 		my @ty = split(/\t/, $e);
# 		if (($ty[0] =~ /^\d/) or ($ty[0]=~ /_/) or ($ty[0]=~ /\//)) { die "\n\n\tThe bar code Ids in $barIn should start with no number and contain no underscorses. Check $option_file.\n\n";}
# 		if ($ty[1] !~ /[ATCGU]/i) {die "\n\n\tYour bar code file $barIn does not seem to be in the right format. Check $option_file.\n\n";}
# 		my $nbr = `wc -l $barIn`;
# 		if ($nbr =~ /(\d.?)\s/) { $nbrBar = $1} 
# 	}
}


## VERSION 3.02 ##
##################### Ajout JCC avril 2013 : by sample order

# Le 17/5/2017 -> absence de codes bars (1 seul échantillon pas de tags)

# if ($barcode==1)
# {
    my $socat = `cat -n $barIn | awk -F"\t" '{ print\$1"\t"\$2 }'`;
    my @socat = split("\n", $socat);
    my %so = ();

    foreach my $line (@socat) {
    my ($val, $key) = split("\t", $line);
    $so{$key} = $val;
    }

# }
# else
# {
# $so{"sample"}=1;
# }

########################### ajout position

my $forw_regexp; my $reve_regexp;
my @tab_forw=split(//,$seq_F);
$forw_regexp=regexp(@tab_forw);

my $rev=reverse($seq_R);
my @tab_reve=split(//,$rev);
$reve_regexp=regexp(@tab_reve);
$reve_regexp=~tr/ACGURMBD/UGCAYKVH/;

sub regexp{
	my $regexp;
	foreach (@_){
		if ($_ eq "T"){
			$regexp.="U";
		}
		elsif ($_ eq "R"){
			$regexp.="[AGR]";
		}
		elsif ($_ eq "Y"){
			$regexp.="[CUY]";
		}
		elsif ($_ eq "M"){
			$regexp.="[CAM]";
		}
		elsif ($_ eq "K"){
			$regexp.="[UGK]";
		}
		elsif ($_ eq "W"){
			$regexp.="[UAW]";
		}
		elsif ($_ eq "S"){
			$regexp.="[CGS]";
		}
		elsif ($_ eq "B"){
			$regexp.="[CUGB]";
		}
		elsif ($_ eq "D"){
			$regexp.="[AUGD]";
		}
		elsif ($_ eq "H"){
			$regexp.="[AUCH]";
		}
		elsif ($_ eq "V"){
			$regexp.="[ACGV]";
		}
		elsif ($_ eq "N"){
			$regexp.="[ACGUN]";
		}
		else{
			$regexp.=$_;
		}
	}
	return $regexp;
}

############################################# trimming profiles ######################################################

if (!(defined $primF) and (defined $primR))  {
	die "\n\n\t Primer forward not defined. Check panam.ini.\n\n";
}

if (defined $primF) {
	die "\n\n\t Primer reverse not defined. Check panam.ini.\n\n" unless (defined $primR);
	die "\n\n\t Primers' sequences missed. Check panam.ini.\n\n" unless ((defined $seq_F) and (defined $seq_R));
	die "\n\n\t Primers' sequences seem to be in the wrong format. Check panam.ini.\n\n" unless (($seq_F =~ /[ATCGUatcgu]/) and ($seq_R =~ /[ATCGUatcgu]/));
	$trim = 1
}

my @tax=();
my $Profiles;
if ($trim == 1) {
	$Profiles = $Reference."/Reference_all/Profiles_".$primF."_".$primR;
	if ((-d $Profiles)) {
		print "Trimmed profiles exist in $Profiles/\n";

	}
	
	else {
		print "\n\tTrimming profiles ...\n"; 

		############### Average positions of the primers
		my %list_profils;
		my $acces; 
		my $taxo;
# 		my @tax = ();
		my $rep = $Reference."/Reference_all/Profiles/fasta/";
		my @listefic = <$rep*.fasta>;
			
		my $pos_forward; my $pos_reverse;
			
		foreach my $nn (@listefic) { 
			#foreach my $domain (@dom){
			foreach my $domain (keys %dom) { # All sequences dor the domain(s) studied
				if ($nn =~ /$domain/) { 
					if ($nn=~ /$Reference\/Reference_all\/Profiles\/fasta\/(.*?).fasta/) { $taxo = $1 ; push (@tax, $taxo); }	
					open (my $profils,"<".$Reference."/Reference_all/Profiles/fasta/$taxo.fasta") || die "can not open file";
					while (<$profils>){
						chomp $_;
						if ($_=~/^>(.*)$/){
							$acces=$1;
						}
						else{
							$_=~s/ //g;
							$list_profils{$taxo}{$acces}.=$_;
						}
					}
				close $profils;
				}
			}
		}
			
		my $pos_reverse_som=0; my $pos_forward_som=0; my $if =0; my $ir=0; 
		foreach my $taxo (@tax){
			my $pos_forward;my $pos_reverse;
			foreach my $acces (keys %{$list_profils{$taxo}}){
				my $sequence=$list_profils{$taxo}{$acces};
		
				$sequence=~s/-//g;
				$sequence=~s/\.//g;
		
				if ($sequence=~/$forw_regexp/){ # casse ?
					$if++;
					$pos_forward=length($`) ;
					$pos_forward_som+=$pos_forward;
				}
				if ($sequence=~/$reve_regexp/){
					$ir++;
					$pos_reverse=length($`) + length($&);
					$pos_reverse_som+= $pos_reverse;
				}
						
				if (defined($pos_reverse) && defined($pos_forward)){
					last;	
				}
			}
		#print LOG "\t\t- Primers positions in $taxo-> $forw_regexp : $if $reve_regexp $ir\n";
		}
		
		

			if(	$if > 0 && $ir > 0){ #  10/01/2017 -> condition possible si un primer trouvé et longueur du primer
			my $posF = int ($pos_forward_som/$if); my $posR = int ($pos_reverse_som/$ir);
			print LOG "\t\t- Average primers positions in  $forw_regexp : $posF $reve_regexp $posR\n";
			#############fin recuperation des positions	moyennes sur les séquences
			trim ($primF, $primR, $posF, $posR);
			print "\n\tTrimmed profiles have been generated in ".$Profiles."\n\n";
			}
			else
			{
			die "Primers not found in reference database, PANAM cancelled...\n\n"; 
			}
	}
}

sub trim {
	my $posF; my $posR; my $primF; my $primR; my $fasta_file; my $hmm_file;
	($primF, $primR, $posF, $posR) = @_;

	my $nom = $primF."_".$primR ;
	
	`mkdir $Reference"/Reference_all/Profiles_"$nom/`;
	`mkdir $Reference"/Reference_all/Profiles_"$nom/hmmprofil/`;
	`mkdir $Reference"/Reference_all/Profiles_"$nom/fasta/`;
	
	open (C, $Reference."/seq_model_trim") || die "cannot open file";
	my %corresp;
	while (<C>) {
		my @tab = split (/\t/, $_);
		chomp ($tab[0]); chomp ($tab[1]);
		$corresp{$tab[0]} = $tab[1];
	}
	close C;
	
	foreach my $n (@tax) {
		open (F, $Reference."/Reference_all/Profiles/fasta/".$n.".fasta") || die "can not open file";
		my %seq ; my $a;
		while (<F>) {
			if ($_=~ />(.*?)\s+/ ) {
				$a = $1;
			}
			else {
				chomp($_);
				$seq{$a}.= $_;
			}
		}
		close F;
		
		my $seq = $corresp{$n}; # Comment est choisi cette séquence de référence ?
		
		my @tab = split (//, $seq{$seq}) ;
		
		my $pos_reel =1;
		my $pos_nuc =1;
		my $pos_nuc_fin;
		
		my $posF_gaped;
		my $posR_gaped;
		
		foreach my $e (@tab) {
		
			if (($e ne "-") and ($e ne ".")) { 
				$pos_nuc ++;
			}
			if ($pos_nuc == $posF) {
				$posF_gaped = $pos_reel;	
			}
			
			elsif ($pos_nuc == $posR) {
				$posR_gaped = $pos_reel;	
				last;
			}
	
			$pos_reel++;
		}
			
		if (!(defined $posR_gaped)) {$posR_gaped = $pos_reel}
			
		my $length = $posR_gaped - $posF_gaped;
		
		open (A, ">".$Reference."/Reference_all/Profiles_".$nom."/fasta/".$n.".fasta") || die "can not open file"; 
		foreach my $k (keys %seq) {
			my $sub = substr ($seq{$k}, $posF_gaped, $length);
			if ($sub=~ /[ATCGUatcgu]/) {
				print A ">$k\n$sub\n"
			}
		}
		
		close A;
			
		$fasta_file = $Reference."/Reference_all/Profiles_".$nom."/fasta/".$n.".fasta" ;
		`$build -g $Reference"/Reference_all/Profiles_"$nom"/hmmprofil/$n".hmmprofil $fasta_file` ;
			
		$hmm_file = $Reference."/Reference_all/Profiles_".$nom."/hmmprofil/".$n.".hmmprofil"
	}
	#return ($fasta_file , $hmm_file)
}
 
if ($trim == 0) {
	$Profiles= $Reference."/Reference_all/Profiles"
}

############################################# trimming profiles END ######################################################

########### hash query sequences 
my %h; my $c;
open (F, $query_seq);
while (<F>){
	if ($_ =~ /^>(.*?)\n/) {
		$c=$1;
		$c =~ s/\s+$//; 
		
	}
	else {
		chomp ($c);
		$h{$c}.=$_ ;
	}
}
close F;

########### hash tax
my %taxeuk;
open (T, $Reference."/Taxonomy") || die "can not open file";
while (<T>) {
	my @e = split (/\t/, $_);
	chomp($e[0]); chomp ($e[1]); chomp ($e[2]);
	$taxeuk{$e[0]}=$e[1]."\t".$e[2];
}

#################################################
########### hash reference sequences + uclust
`mkdir $NGS_id_Results/$panam_output"/Similarity_annotation"` ;

if (($us_version eq "1.1.579q") or ($us_version eq "3.0.617")) {# DD Modif 23/12
`$path --sort $query_seq --output $NGS_id_Results/$panam_output/Similarity_annotation/seq_sorted_seed 2>> $errorlog` ; 
}# Modif DD 23/12 usearch -sortbylength seqs.fasta -output seqs_sorted.fasta
else
{
`$path -sortbylength $query_seq -output $NGS_id_Results/$panam_output/Similarity_annotation/seq_sorted_seed 2>> $errorlog` ;
}


print "Similarity assignment ...\n";

my $repp = $Reference."/Reference_all/Reference_bases/";
my $e; my %seq;
# 
######################################"""""""
#####comparer aux bases seed

my $rep_seed = $Reference."/Reference_all/Reference_bases/seed_bases/";
#my @listeseed = <$rep_seed*_seed> ; ### modif 28/11/2014 seed concat
my @listeseed = <$rep_seed*_concat> ;  ### modif 28/11/2014 seed concat
my %seed_affil;
foreach (@listeseed) {
	my $seed_file = $_;

	#if ($seed_file =~ /\/Reference_all\/Reference_bases\/seed_bases\/(.*?)_seed/ ){    ### modif 28/11/2014 seed concat
	if ($seed_file =~ /\/Reference_all\/Reference_bases\/seed_bases\/(.*?)_concat/ ){   ### modif 28/11/2014 seed concat
		my $seed = $1;
		print "Sort by domain  --- $seed\n"; # DD
		
		if (($us_version eq "1.1.579q") or ($us_version eq "3.0.617")) {
			`$path --input $NGS_id_Results/$panam_output/Similarity_annotation/seq_sorted_seed --lib $seed_file --uc $NGS_id_Results/$panam_output/Similarity_annotation/seq_results_seed_"$seed".uc --id 0.10 --libonly --allhits --maxaccepts 1 --rev 2>> $errorlog`;
		}


		else {# Version 7.0 et plus DD 23/12
			`$path -usearch_global $NGS_id_Results/$panam_output/Similarity_annotation/seq_sorted_seed -db $seed_file -uc $NGS_id_Results/$panam_output/Similarity_annotation/seq_results_seed_"$seed".uc  -id 0.10 -maxaccepts 1 -strand both 2>> $errorlog`;

# print "$path -usearch_global $NGS_id_Results/$panam_output/Similarity_annotation/seq_sorted_seed -db $seed_file -uc $NGS_id_Results/$panam_output/Similarity_annotation/seq_results_seed_"$seed".uc  -id 0.10 --maxaccepts 1 \n";
		}	

		open (KK, $NGS_id_Results."/".$panam_output."/Similarity_annotation/seq_results_seed_".$seed.".uc") || die "can not open file" ;
# L	258	1736	*	*	*	*	*	eukaryota_AM501959	*
# H	258	444	44.6	+	0	0	529I75M3D28M7D113M16D47M28D127M817I	ARNAydat_OTU1924	eukaryota_AM501959
# H	258	442	45.5	+	0	0	529I85M5D78M2D53M17D47M28D127M817I	ARNFades_OTU7287	eukaryota_AM501959
# L	278	1721	*	*	*	*	*	eukaryota_AB546139	*
# H	278	442	48.9	+	0	0	527I82M6D22M6D13M32D30MD73M3D67MD106M801I	ADNVichy_OTU10812	eukaryota_
# AB5
		while (<KK>) {
			if (!($_=~ /\#/)) {
				my @t = split (/\t/, $_);
				if ($t[0] eq "H") {
					chomp($t[8]); chomp($t[3]); 
					$seed_affil{$t[8]}{$seed} = $t[3];

					#########  ajout 28/11/2014 seed concat
					my @l = split(/_/, $t[9]);
					$concat{$t[8]}{'domain'} = $l[0]; # DD 0.1 pour le domaine ?
					########

				}
			}
		} 
		close KK;		
	}
}
 ### modif 28/11/2014 seed concat

=pod
my %prim_affil ;
foreach my $u (keys %seed_affil) {
	my $max = 0;
	foreach my $uu (keys %{$seed_affil{$u}}) {
		if ($seed_affil{$u}{$uu} > $max) { 
			$max = $seed_affil{$u}{$uu}	; 		
			$prim_affil{$u} = $uu ;
		}
	}
}			

#open (DISC, ">".$NGS_id_Results."/".$panam_output."/Similarity_annotation/discarded.fasta") || die "can not open file" ;
foreach my $k (keys %prim_affil) {
	if (!(exists $dom{$prim_affil{$k}})) {
		#print DISC ">$k\n$h{$k}";
		delete ($h{$k});
	}
}
#close DISC;

open (QUERY, ">".$NGS_id_Results."/".$panam_output."/Similarity_annotation/query_seq") || die "can not open file" ;
foreach my $mm (keys %h) {
	print QUERY ">$mm\n$h{$mm}"
}		
close QUERY;
=cut

############## concat -> les séquences qui n'appartiennent pas au domaine sont enlevées

open (DISC, ">".$NGS_id_Results."/".$panam_output."/Similarity_annotation/discarded.fasta") || die "can not open file" ;
foreach my $k (keys %concat) {
	if (!(exists $dom{$concat{$k}{'domain'}})) {
		print DISC ">$k\n$h{$k}";
		delete ($h{$k});
	}
}
close DISC;	



open (QUERY, ">".$NGS_id_Results."/".$panam_output."/Similarity_annotation/query_seq") || die "can not open file" ;
foreach my $mm (keys %h) {
	print QUERY ">$mm\n$h{$mm}"
}
close QUERY;

### fin modif 28/11/2014 seed concat



if (($us_version eq "1.1.579q") or ($us_version eq "3.0.617")) 
{# DD Modif 23/12
`$path --sort $NGS_id_Results/$panam_output/Similarity_annotation/query_seq --output $NGS_id_Results/$panam_output/Similarity_annotation/seq_sorted 2>> $errorlog` ; 
}
else #DD 23/12
{
`$path -sortbylength  $NGS_id_Results/$panam_output/Similarity_annotation/query_seq -output $NGS_id_Results/$panam_output/Similarity_annotation/seq_sorted 2>> $errorlog` ;
}
###################################################"" fin comparaison aux seed bases

print "Similarity search in progress...\n"; # DD 23/12

my @listeficc = <$repp*.fasta>;
my $num_uc=1;
foreach(@listeficc) { 
	my $ficTemp = $_; 
	#foreach my $domain (@dom){	
	foreach my $domain (keys %dom) {
		if ($ficTemp =~ /$domain/) {
		print "Domain $ficTemp\n"; # DD 23/12

		#################################################
			if (($us_version eq "1.1.579q") or ($us_version eq "3.0.617")) {
				`$path --input $NGS_id_Results/$panam_output/Similarity_annotation/seq_sorted --lib $ficTemp --uc $NGS_id_Results/$panam_output/Similarity_annotation/seq_results_"$num_uc".uc --id 0.10 --libonly --allhits --maxaccepts 5 --rev 2>> $errorlog`;
			}

#  usearch -usearch_global  preprocess_output/pooled_sample/pooled_sample_SEQ_OTU  -db /usr/local/panam/panam_v4/Reference_nvxProfils_bact115_euk108_PlusSeqEnv/Reference_all/Reference_bases/usearch_Reference_base_totale.db -id 0.2 -maxaccepts 5 -maxhits 5 -strand plus -uc_allhits -uc results.uc 
			else { # DD 23/12
# 				`$path --query $NGS_id_Results/$panam_output/Similarity_annotation/seq_sorted --db $ficTemp --uc $NGS_id_Results/$panam_output/Similarity_annotation/seq_results_"$num_uc".uc  --id 0.10 --maxaccepts 5 --rev --local 2>&1`
			      `$path -usearch_global $NGS_id_Results/$panam_output/Similarity_annotation/seq_sorted -db $ficTemp -uc $NGS_id_Results/$panam_output/Similarity_annotation/seq_results_"$num_uc".uc  -id 0.10 -maxaccepts 5 -maxhits 5 -strand both -uc_allhits 2>> $errorlog`;

			}
		##################################################


			`cat $NGS_id_Results/$panam_output/Similarity_annotation/seq_results_"$num_uc".uc >> $NGS_id_Results/$panam_output/Similarity_annotation/seq_results.uc` ;

			$num_uc++;
			open(Q, "$ficTemp");
			while (<Q>){
				if ($_ =~ /^>(.*?)\s+/) {
					$e=$1;
				}
				else {
					chomp ($_);
					$seq{$e}.=$_ ; # DD séquences de référence non alignées ????
				}
			}
			close Q ;
		}
	}
}


####################### sorting sequences in profiles defined by sortingtax
 
if (!(-z $NGS_id_Results."/".$panam_output."/Similarity_annotation/seq_results.uc")) {
	print  "vsearch done\n"; # DD sinon quoi ??? 11/05/2016
}

##le premier tableau @tax n'existe pas toujours, si les profiles existent déjà, @tax n'est pas créé. 
print "Sorting sequences ...\n"; 
my @tax = ();
open (IL, $Reference."/sortingtax");
@tax=<IL>;
close IL;


open (F, $NGS_id_Results."/".$panam_output."/Similarity_annotation/seq_results.uc") || die "can not open file $NGS_id_Results/$panam_output/Similarity_annotation/seq_results.uc" ;
# H	59066	435	82.2	+	0	0	551I183M42D39MI59MI112M804I	ADNAydat_OTU16	DQ514895
# H	63117	435	82.2	+	0	0	513I183M42D39MI59MI112M804I	ADNAydat_OTU16	DQ514911
# H	62135	435	82.2	+	0	0	555I183M42D39MI59MI112M848I	ADNAydat_OTU16	GQ844874

open (B, ">".$NGS_id_Results."/".$panam_output."/Similarity_annotation/best_hit_uc");
# ADNRocheMoines_OTU2778	FJ559386	Eukaryota;Viridiplantae;Chlorophyta;Chlorophyceae;Chlorococcales;Chlorococcaceae;Chlorococ
# cum;	Chlorophyta
# ADNRocheMoines_OTU2778	GQ122363	Eukaryota;Viridiplantae;Chlorophyta;Chlorophyceae;Chlorococcales;Chlorococcaceae;Chlorococ
# cum;	Chlorophyta

my %hits; my %target_hits;
while (<F>) {
	if (!($_=~ /\#/)) {
		my @t = split (/\t/, $_);
		if ($t[0] eq "H") {
			chomp ($t[8]) ,chomp($t[9]); chomp ($t[3]) ;
			$t[8] =~ s/\s+$//; $t[9] =~ s/\s+$//; $t[3] =~ s/\s+$//;
			$target_hits{$t[8]}{$t[9]} = $t[3];

			if ($t[4] eq "-") { # DD Ne devrait pas exister 23/12/2015 -> pb de NN dans la base ou de séquences à éliminer....
			print LOG "OTUs strand - $t[8]\n";
# 				my $rev = reverse ($h{$t[8]});
# 				$rev =~ tr/ACGTUacgtu/TGCAAtgcaa/ ; ##### modif 24/5
# 				chomp ($rev);
# 				$h{$t[8]} =$rev;
			}
		}
	}
}
 
my %tax_hits; my %seq_hits;
foreach my $query (keys %target_hits) {
	my $nb_target=0;
	foreach my $target (sort {$target_hits{$query}{$b} <=> $target_hits{$query}{$a}} keys %{$target_hits{$query}}) { # DD Classement en fonction de l'identité
		my $ligne_taxo = $taxeuk{$target};
		my @t = split (/\t/, $ligne_taxo);
		print B "$query\t$target\t$t[0]\t$t[1]\t$target_hits{$query}{$target}\n"; # DD 11/1/2016 Ajout identité ici si $target_hits{$query}{$target} version > 7 
		$tax_hits{$query}{$nb_target} = $t[1];
		$seq_hits{$query}{$nb_target} = $target;
		$nb_target++;
		if ($nb_target == 5) {last}
	}
}

close B;




`mkdir $NGS_id_Results/$panam_output"/Similarity_annotation/Sorted_Sequences"`;
foreach my $tax (@tax) { # DD Tous les domaines/profiles  sont parcourus ????
	chomp ($tax);
	my %tax_bis;
	if (-e $Profiles."/fasta/".$tax.".fasta") 
	{ 
		open (T, $Profiles."/fasta/".$tax.".fasta") || die "can not open file"; 
		my $aa; my %seq_ref;
		while (<T>) {
			if ($_=~ />(.*?)\s+/) {
				$aa = $1;
				chomp ($aa);
				$seq_ref{$aa}=1;
			}
		}
		close T;
		
		if (!(-d $NGS_id_Results."/".$panam_output."/Similarity_annotation/Sorted_Sequences/files_".$tax)) {
			`mkdir $NGS_id_Results/$panam_output"/Similarity_annotation/Sorted_Sequences/files_"$tax`;
		}

		open (TAX, ">".$NGS_id_Results."/".$panam_output."/Similarity_annotation/Sorted_Sequences/files_".$tax."/seq_".$tax);
		
	
		foreach my $seq_c (keys %seq_hits) { # DD : parcours du tableau best hit $seq_hits{ARNAydat_OTU1924}{$nb_target} = $target =AY485487
# 			if ($tax_hits{$seq_c}{'0'} =~ /$tax/) {
			if ($tax =~ /$tax_hits{$seq_c}{'0'}/) { # DD $tax_hits{ARNAydat_OTU1924}{$nb_target} = profil dans la taxonomie
				print TAX ">$seq_c\n$h{$seq_c}\n"; # DD impression de l'OTU et de sa séquence dans le profil 
		
				my %aze;
				foreach my $e (keys %{$seq_hits{$seq_c}}) { # DD Parcours de 1 à 5 pour l'OTU analysé
					$aze{$seq_hits{$seq_c}{$e}}= 1
					
				}
	
				foreach my $uu (keys %seq_ref) { # DD Les séquences du profils sont parcourues 
					delete ($aze{$uu}) # DD si la séquence fait partie du profil elle est éliminée sinon elle est gardée même si elle n'a pas la bonne taxonomie ??? -> elle ne fait pas partie du profil car similarité très éloignée !!!!!! -> on introduit du bruit !!!???
				}
	
				foreach my $qs (keys %aze) {
					$tax_bis{$qs} = 1;
				}
			}
		}
	
		foreach my $ww (keys %tax_bis) {
			print TAX ">$ww\n$seq{$ww}\n"; # DD Ecriture de la séquence non alignéee
		}
		close TAX; # DD OTU + séquences proches complètes   ??? Pour utiliser SINA il faut couper les séquences ici avant d'introduire les séquences de ref coupées ou alors aligner sur le profil complet comme pou Arb!!! Les séquencees de ref sont à la fin donc si fichier coupé en 2 -> les OTUs au début ????
	}
}
# 	
###################################### cutting ...

my $max=0;
foreach my $e (@tax) { 
	my $aa; my %hh; my $i= 0; my $num=1; 
	if (-e $NGS_id_Results."/".$panam_output."/Similarity_annotation/Sorted_Sequences/files_".$e."/seq_".$e) {

		open (AA,  $NGS_id_Results."/".$panam_output."/Similarity_annotation/Sorted_Sequences/files_".$e."/seq_".$e);
		while (<AA>) {
			if ($_=~/^(>.*?)\s/) {
				$aa = $1;
			}
			else {
				$hh{$aa}.=$_;
			}	
		}
		
		open (P, ">". $NGS_id_Results."/".$panam_output."/Similarity_annotation/Sorted_Sequences/files_".$e."/file_cut_".$num);

		foreach my $ee (keys %hh) {
			$i++;
			if ($i != 27001) { # DD 27000 séquences introduites ??
 				print P "$ee\n$hh{$ee}\n";
			}
			else  {
				close P;
				$num++;
				open (P, ">". $NGS_id_Results."/".$panam_output."/Similarity_annotation/Sorted_Sequences/files_".$e."/file_cut_".$num);
				print P "$ee\n$hh{$ee}\n";
				$i= 1;
				
			}
		}
		if ($num > $max) { $max = $num; }
		close P;
	}		
}

my @num = (1 .. $max);

# #################################""
 
`mkdir $NGS_id_Results/$panam_output"/Phylogeny"` ;
`mkdir $NGS_id_Results/$panam_output"/Alignment"`;
`mkdir $NGS_id_Results/$panam_output"/Phylogeny/Node_files"`;


################## changé /!\ 28/3/2012
open (AFFIL, ">".$NGS_id_Results."/PANAM_Affiliation.txt");		## VERSION 3.03 ## ">>" changed to ">" 

###############"ajout /!\ fichiers clades 28/3/2012
open (CLADENN, ">".$NGS_id_Results."/PANAM_Clades_NN.txt") || die "can not open file PANAM_Clades_NN.txt";
open (CLADELCA, ">".$NGS_id_Results."/PANAM_Clades_LCA.txt") || die "can not open file PANAM_Clades_LCA.txt";
 my %taxOTU; my %leveltax; my %sample;
################" fin ajout fichiers clades

##########ajout 24/5
my %seq_affil;
######################" fin ajout 24/5

foreach my $kk (@tax) {
	foreach my $fff (@num) {
#########Alignment
		if (-e $NGS_id_Results."/".$panam_output."/Similarity_annotation/Sorted_Sequences/files_".$kk."/file_cut_".$fff && !(-z $NGS_id_Results."/".$panam_output."/Similarity_annotation/Sorted_Sequences/files_".$kk."/file_cut_".$fff)) { # DD 18/01/2016 Fichier non vide et existe
			my $ca = $align." -o ".$NGS_id_Results."/".$panam_output."/Alignment/".$kk.".".$fff.".fasta --outformat a2m --mapali ".$Profiles."/fasta/".$kk.".fasta ".$Profiles."/hmmprofil/".$kk.".hmmprofil ". $NGS_id_Results."/".$panam_output."/Similarity_annotation/Sorted_Sequences/files_".$kk."/file_cut_".$fff;	
			system "$ca  2>> $errorlog"  ;
			print "Alignment files_$kk.$fff done\n";
			print LOG "$ca\nAlignment files_$kk.$fff done\n";
			@clean=clean_alignments($NGS_id_Results."/".$panam_output."/Alignment/", $kk.".".$fff); #DD 11/1/2016
			print "Length of filtered alignment for building phylogeny: $clean[0]\n";# DD 26/08/2016
                        print LOG "Length of filtered alignment for building the phyleny: $clean[0]\n";
    
		}
	
#######Phylogeny	
		if (-e $NGS_id_Results."/".$panam_output."/Alignment/".$kk.".".$fff.".fasta") {
			my $ct = $fast." -nt -boot 100 ".$NGS_id_Results."/".$panam_output."/Alignment/".$kk.".".$fff.".fasta > ".$NGS_id_Results."/".$panam_output."/Phylogeny/".$kk.".".$fff.".newick";
			system "$ct  2>> $errorlog" ;
			print "Phylogeny $kk.$fff done\n";
			print LOG "Phylogeny $kk.$fff done\n";
		}

######## Rooting tree
		if ( -e $NGS_id_Results."/".$panam_output."/Phylogeny/".$kk.".".$fff.".newick"){
#18/05/2016
		  if (exists ($dom{'enterovirus'})){
		  root_enterovirus($NGS_id_Results, $panam_output, $kk, $fff);
		  }
		  else{
		  root_life($NGS_id_Results, $panam_output, $kk, $fff);
		  }
		}
#18/05/2016
	
# 		`rm $NGS_id_Results/$panam_output/n_tree`;

#############Nodes_file
		if ( -e $NGS_id_Results."/".$panam_output."/Phylogeny/".$kk.".".$fff."_rooted.newick"){
	
			my $nom_fichier = $NGS_id_Results."/".$panam_output."/Phylogeny/".$kk.".".$fff."_rooted.newick"  ;
			open (E, ">".$NGS_id_Results."/".$panam_output."/Phylogeny/Node_files/file_get_descendent_seq_".$kk.".".$fff);
	
			my $input = new Bio::TreeIO -> new(-file   => $nom_fichier,
				-format => "newick");
			
			while (my $tree = $input->next_tree) {
				for my $node ( grep { ! $_->is_Leaf } $tree->get_nodes) {
					print E" Node : ", $node->id ;
					for my $child ( $node -> get_Descendents ) {
						print E" child : ", $child->id ;
					}
					print E"\n";
				}
			}
			close E;
		}
		########################### parsing
		if ( -e $NGS_id_Results."/".$panam_output."/Phylogeny/Node_files/file_get_descendent_seq_".$kk.".".$fff){
			open (W, $NGS_id_Results."/".$panam_output."/Phylogeny/Node_files/file_get_descendent_seq_".$kk.".".$fff) ;
			my $i=0;
			my %hh; my %child; my %node_seq; my %boot;
			while (<W>) {
				if ($_=~ /Node :\s(.*?)\schild/) { 
					my $a = "Node_".$i; 
					my $boot_value = $1;
					my $child = ($_=~ tr/child//);
					$child{$a}= $child; 
					$child{$a}{'bootstrap'} = $boot_value;
	
					$boot{$a}= $boot_value;
					my @tab = split (/:/, $_);
					foreach my $e (@tab) { 
						my @t = split (/ /,$e);
						foreach my $p (@t) { 
							my $pp; 
							if (($p=~ /^\w/) and ($p !~ /\./)) { 
								chomp ($p);
								if (($p !~ /Node/) and ($p !~ /child/)){
									$pp = $p;
								}
								if (defined $pp ) { 
									$node_seq{$pp} .="\t".$a; 
									my @tax_seq = split(/\t/, $taxeuk{$pp});
									my $tax_seq = $tax_seq[0];
									if (defined ($tax_seq)) {
										chomp ($tax_seq); 
										$hh{$a}{'tax'}.="\t".$tax_seq;
										$hh{$a}{'seq'}.="\t".$pp;		
									}	
								}
							}
						}
					}
				}
				$i++;
			}
	
	##############"NN	
			my %nn_node;
			foreach my $nn (keys %hh) {
				my @mmm = split (/\t/, $hh{$nn}{'tax'});
				if ($mmm[1] ne "") {
					$nn_node{$nn} = $mmm[1]
				}
			}
	
	###################"LCA
			my %taxo;
			foreach my $p (keys %hh) {	
				my @ttax=();
				my @tab = split (/\t/,$hh{$p}{'tax'}); 
				my $pp = $#tab;
				
			##### toutes les taxonomies reliées à un noeud sont stockées dans @t; 
		
				for (my $i=0; $i<=$pp; $i++) {
					if ($tab[$i] ne "") {
						push(@ttax, $tab[$i]);
					}
				}
				$taxo{$p}{'tax'} = \@ttax;
			}
	
			##########################

			my %taille_taxo;
			foreach my $bb (keys %taxo) {
				my $lmin_silva=100000; 
		
				my $mm = $#{$taxo{$bb}{'tax'}};
				my @r; my @rr; my $nnn ; my $nn;
				
				for(my $i=0; $i<=$mm; $i++) {
					my @r = split(/;/, ${$taxo{$bb}{'tax'}}[$i]);
					my $nn = $#r;
					if (($nn >=0) and ($nn < $lmin_silva)) {
						$lmin_silva = $nn;
					}
				}
					
				$taille_taxo{$bb}{'tax'} = $lmin_silva;
			}

			my %term; my %taxo_com; my $qq; my @tt; my %rang_taxo;
				
			foreach my $p (keys %taxo) {
				my $qs = $#{$taxo{$p}{'tax'}}; 
				for (my $oo =0; $oo<= $qs; $oo++) { 
					my @t = split (/;/, ${$taxo{$p}{'tax'}}[$oo]); 
					for (my $i = 0; $i <= $taille_taxo{$p}{'tax'}; $i++) {
						if (exists $t[$i]) { #print "ti --- $t[$i]\n";
							$term{$p}{'tax'}{$i}.="\t".$t[$i];
						}
					}	
				}
			}
	
			my %tax_node;
			foreach my $tt (keys %term) { 
				foreach my $cc (keys %{$term{$tt}}) { 
					$tax_node{$tt}{$cc}=""; 
					my %temp=%{$term{$tt}{$cc}};
					my @liste_level=sort keys(%temp);
					foreach(@liste_level){
						my $ss=$_; 
						my %rang_tax;		
						if (defined $term{$tt}{$cc}{$ss}) {
							my @ff = split(/\t/, $term{$tt}{$cc}{$ss});
							my $t_max = $#ff; 
							foreach my $ll (@ff) {
								$rang_tax{$ll}++;
							} 
				
							foreach my $aa (keys %rang_tax) {
								if ($aa ne "") { 
									if ($rang_tax{$aa} == $t_max) { 
										$tax_node{$tt}{$cc}.=$aa.";";
									}
								}
							}
						}
					}
				}
			} 
	
			################ printing results 
			open (R, ">".$NGS_id_Results."/PANAM_Affiliation".$kk."_".$fff.".txt");
			if (-e $NGS_id_Results."/".$panam_output."/Similarity_annotation/Sorted_Sequences/files_".$kk."/file_cut_".$fff) {
				open (G, $NGS_id_Results."/".$panam_output."/Similarity_annotation/Sorted_Sequences/files_".$kk."/file_cut_".$fff)  ;
				my %nom_seq;
				while (<G>) {	
					if ($_=~ /^>(.*?)\s+/) {
						my $aa = $1; 
						$nom_seq{$aa} = 1;
					}
				}

				foreach my $e (keys %nom_seq) { 
				##############"ajout 24/5	
#					$seq_affil{$e}=1;
				###########fin ajout 24/5  
					if (exists $h{$e}) {
						my $min = 1000000;
						my $res_node;
						my @tab = split (/\t/, $node_seq{$e});
						foreach my $a (@tab) {
							if ($a ne "") {
								if (($child{$a} < $min) and defined ($nn_node{$a}) ) {
									$min = $child{$a}; 
									$res_node = $a; 
								}
							}
						}
	
						my @mmm = split (/\t/,$hh{$res_node}{'seq'});
						my $nearest = $mmm[1];
						my @ttax_nearest = split(/\t/,$taxeuk{$nearest});
						my $tax_nearest	= $ttax_nearest[0];
						
						my $tax_node;
						if ($tax_node{$res_node}{'tax'} eq "") { #print "res_node --- $res_node\n"; 
							$tax_node = "Other;"}	
						else { $tax_node = $tax_node{$res_node}{'tax'}	 }
				
						########## print dans 2 fichiers : le final AFFIL et le provisoir $k

						my @ppp = split (/;/, $tax_nearest);
						my $ppp = $#ppp;
						if ($ppp == 1) {
							$tax_nearest .="unclassified ".$ppp[1].";"
						}
						
						my @ooo = split (/;/, $tax_node);
						my $ooo = $#ooo;
						if ($ooo == 1) { 
							$tax_node .="unclassified ".$ooo[1].";" ; 
						}

						print R ">Read: $e\tBootstrap:$boot{$res_node}\nNearest neighbor: $nearest\t$tax_nearest\nLowest node: $res_node\t$tax_node\n\n";				
						print AFFIL ">Read: $e\tBootstrap:$boot{$res_node}\nNearest neighbor: $nearest\t$tax_nearest\nLowest node: $res_node\t$tax_node\n\n";

					}
				}
			}
			close R;

#####################################################################

			if (-e $NGS_id_Results."/PANAM_Affiliation".$kk."_".$fff.".txt") { 
				open (AA, $NGS_id_Results."/PANAM_Affiliation".$kk."_".$fff.".txt");
	
	 			my $c; my %hhh; my $bb;  my %bootstrap; 
				while (<AA>) {	
					if ($_=~ /^>Read: (.*?)_/) {    
						my $aa = $1;             
						$sample{$aa}=1;          
					}                               
					if ($_=~ /^>Read: (.*?)\t(.*?)\n/) {
						$c = $1; 
						$bootstrap{$c} = $2;
					}
					else {
						chomp($_);
						$hhh{$c}.= $_
					}
				}
				close AA; 
# 			}

				my %clade;my %nbr; my %tax; my %boot;

				foreach my $f (keys %hhh) { 
					if (($hhh{$f} !~ /Metazoa/) and ($hhh{$f} !~ /Other/)) { # DD à quoi ça sert ???
						$bb = $bootstrap{$f};
						my @tt = split (/Lowest node:/, $hhh{$f}) ;
						chomp ($tt[0]); chomp ($tt[1]);
			
						################"" %taxOTU pour user_file = pooled et %leveltax pour user_file = eachSample
			
						my @uu = split (/\t/, $tt[1]);
						my $node = $uu[0]; 
		
						$boot{'node'}{$node} = $bb;
						$clade{'node'}{$node}.=$f.", ";
						$nbr{'node'}{$node}++;
						$tax{'node'}{$node} =  $uu[1];
			
						if ((defined $uu[1])  and ($uu[1] ne "")) {
							$taxOTU{$f}{'LCA'} = $uu[1]; 
							chomp($uu[1]); 

							$uu[1]=~s/\t//gi;
							$uu[1]=~s/^\s+//gi;
							$uu[1]=~s/\s+$//gi;

							my @aqs = split (/;/, $uu[1]);
							my $aqs = $#aqs; 
	
							if ($aqs == 0) { }
				
							else {
								$leveltax{'LCA'}{$aqs[1]}{$aqs[2]}++
							}
						}
			
						my @u = split (/\t/, $tt[0]);
						if ((defined $u[1])  and ($u[1] ne "")) {
							$taxOTU{$f}{'NN'} = $u[1];
							chomp($u[1]); 
					
							my @bqs = split (/;/, $u[1]); 
							my $bqs = $#bqs;

							if ($bqs == 0) { }
							else {
								$leveltax{'NN'}{$bqs[1]}{$bqs[2]}++
							}
						}

		###########################################"

						my @t = split (/Nearest neighbor:/, $tt[0]);
				
						my @proc = split (/\t/, $t[1]); 
						my $voisin = $proc[0]; 
			
						$boot{'voisin'}{$voisin} = $bb;
						$tax{'voisin'}{$voisin}=$proc[1]."\n\t\t";
						$clade{'voisin'}{$voisin}.=$f.", ";
						$nbr{'voisin'}{$voisin}++;	
					}
				}

				foreach my $p (keys %{$clade{'voisin'}}) {
					chomp ($tax{'voisin'}{$p}); 
					if ($nbr{'voisin'}{$p} == 1 ) {
						print  CLADENN ">$p\t$tax{'voisin'}{$p}Phylogeny:  ".$kk.".".$fff."_rooted.newick\t".$boot{'voisin'}{$p}."\n\n\t".$nbr{'voisin'}{$p}." neighbor: ".$clade{'voisin'}{$p}."\n\n";					
					}
					else {
						print  CLADENN ">$p\t$tax{'voisin'}{$p}Phylogeny: ".$kk.".".$fff."_rooted.newick\t".$boot{'voisin'}{$p}."\n\n\t".$nbr{'voisin'}{$p}." neighbors: ".$clade{'voisin'}{$p}."\n\n";
					}
				}

				foreach my $p (keys %{$clade{'node'}}) { 
					chomp ($tax{'node'}{$p}); 
					if ($nbr{'node'}{$p} == 1 ) {
						print  CLADELCA ">$p\t$tax{'node'}{$p}\n\t\tPhylogeny:  ".$kk.".".$fff."_rooted.newick\t".$boot{'node'}{$p}."\n\n\t".$nbr{'node'}{$p}." Descendent: ".$clade{'node'}{$p}."\n\n";
					}
					else {
						print  CLADELCA ">$p\t$tax{'node'}{$p}\n\t\tPhylogeny:  ".$kk.".".$fff."_rooted.newick\t".$boot{'node'}{$p}."\n\n\t".$nbr{'node'}{$p}." Descendents : ".$clade{'node'}{$p}."\n\n";
					}
				}

			}
			`rm $NGS_id_Results/PANAM_Affiliation$kk"_"$fff.txt`;			## VERSION 3.03 ##
		}
	}
}
close AFFIL; 			## VERSION 3.03 ## 
close CLADENN;			## VERSION 3.03 ## 
close CLADELCA;			## VERSION 3.03 ## 

################################ eliminer les seq qui n'appartiennent pas aux groupes d'intérêt

my @dom_enlv = ("eukaryota", "bacteria" , "archaea" , "Metazoa", "Other"); # DD Les Embriophyta ??? A vérifier, que signifie other ?

foreach my $k (keys %taxOTU) { #print "k -- $k\n";
	foreach my $aff (keys %{$taxOTU{$k}}) { 
		foreach my $d (@dom_enlv) { 
			if (exists $dom{$d}) {} # DD : Cette variable est égale à 1 si le domaine est déclaré
			else {
				my @tab = split(/;/, $taxOTU{$k}{$aff}) ;
				if ($tab[1] eq $d ) {
					delete ($taxOTU{$k}{$aff})
				}

				#### ajout 16/9/2014
				if ($tab[2] eq "Chloroplast" ) {
					delete ($taxOTU{$k}{$aff})
				}
				#### fin ajout 16/9/2014
			}
		}
	}
}
#################################
##### récupérer les OTUs avec seed, children et nombre ; selon le type d'analyse (Pooled ou eachSample)
##### récupérer $nbrSeq pour calcul des indices

my %nbrSeq;
my %seed; my %child;  
my %barIn;
my %NbrSeqClean;

if ($user_file eq "pooledSamples") { 
	open (FF, $NGS_id_Results."/".$preprocess_output."/pooled_sample/pooled_sample_OTU") || die "can not open file";
	#my %barIn;   ### modif 16/9/2014
	my %hps; 
	my $l;
	while (<FF>) {
		#if ($_=~ /^OTU/) {} #### modif 20/8/2014
		if ($_=~ /seed Sequence/) {} #### modif 20/8/2014
		else {
			chomp ($_);
			my @tab = split (/\t/, $_);

			############################################ !!!!!!!!!!!
			#if (exists $taxOTU{$tab[1]}{'LCA'}) {  ##### modif 20/8/2014
			if (exists $taxOTU{$tab[0]}{'LCA'}) {  ##### modif 20/8/2014

				my $id = $tab[0]."\t".$tab[1];
				#$hps{$id}{'id'} = $tab[1]; ##### modif 20/8/2014
				$hps{$id}{'id'} = $tab[0]; ##### modif 20/8/2014
				my @tt = split (/, /, $tab[3]);			
				
				foreach my $e (@tt) {
					my $d;
					if ($e=~ /(\w.*?)_/) { 
						$d = $1; 
						########
					 	$NbrSeqClean{$d}++;
						#####
						#$barIn{$d} = 1; #### modif 16/9/2014
						$hps{$id}{$d}++;
						$nbrSeq{$hps{$id}{'id'}}{$d}{'pooled'}++ ;
					}
				}	
		
				my $u = scalar(@tt);
				$u-= 1;
				my @t = @tt[0..$u];

				#$child{$tab[1]} = \@t; ##### modif 20/8/2014
				$child{$tab[0]} = \@t; ##### modif 20/8/2014
				foreach my $e (@t) {	
					if (defined $e) {
						#$seed{$e} = $tab[1];   ##### modif 20/8/2014
						$seed{$e} = $tab[0];   ##### modif 20/8/2014
					}
				}
			} 
			#else { delete $taxOTU{$tab[1]} }	### éliminer les seq qui n'ont pas été affiliées phylogénétiquement   ### modif 20/8/2014
			else { delete $taxOTU{$tab[0]} }	### éliminer les seq qui n'ont pas été affiliées phylogénétiquement   ### modif 20/8/2014
			
		}
	}
	close FF;

#####################################" OTU_distribution_taxo ###################################################################################

		
	open (REPS, ">".$NGS_id_Results."/OTU_distribution_tax.txt");
 	print REPS "OTU_Id\tOTU_Seed\t";
	#foreach my $pp ( (sort {$so{$a} <=> $so{$b}} keys %barIn)) { 	#### modif 16/9/2014
	foreach my $pp ( (sort {$so{$a} <=> $so{$b}} keys %so)) {								
		print REPS "$pp\t";
	}
	print REPS "NN taxonomy\tLCA taxonomy\n"; 
		
	foreach my $p (keys %hps) {
		print REPS "$p\t";
		#foreach my $pp ( (sort {$so{$a} <=> $so{$b}} keys %barIn)) {   #### modif 16/9/2014
		foreach my $pp ( (sort {$so{$a} <=> $so{$b}} keys %so)) {								
	 		if (defined $hps{$p}{$pp}) {
				print REPS"$hps{$p}{$pp}\t";
			}
			else {
				print REPS"0\t"
			}
		}
		print REPS "$taxOTU{$hps{$p}{'id'}}{'NN'}\t$taxOTU{$hps{$p}{'id'}}{'LCA'}\n" ; #20/05/2016
	}

	close REPS;
}


if ($user_file eq "eachSample") {
	foreach my $s (keys %sample) { #print "s --- $s\n";
		open (FFF, $NGS_id_Results."/".$preprocess_output."/".$s."_OTU") || die "can not open file" ;
		while (<FFF>) {
			if ($_=~ /^OTU/) {}
			else {
				my @tab = split (/\t/, $_);
	
				############################################ !!!!!!!!!!!  /!\ NN ou LCA  dans le test d'existence ????
				#if (exists $taxOTU{$tab[1]}{'LCA'}) {    ##### modif 20/8/2014
				if (exists $taxOTU{$tab[0]}{'LCA'}) {    ##### modif 20/8/2014
					###### $nbrSeq
					#$nbrSeq{$tab[1]}{$s}{'each'} = $tab[2];    ##### modif 20/8/2014
					$nbrSeq{$tab[0]}{$s}{'each'} = $tab[2];    ##### modif 20/8/2014
					#########			
				
					my @tt = split (/, /, $tab[3]);
					my $u = scalar(@tt);
					$u-= 1;
					my @t = @tt[0..$u];
					#$child{$tab[1]} = \@t;   ##### modif 20/8/2014
					$child{$tab[0]} = \@t;   ##### modif 20/8/2014
					foreach my $e (@t) {
						if (defined $e) {
							#$seed{$e} = $tab[1]   ##### modif 20/8/2014
							$seed{$e} = $tab[0]   ##### modif 20/8/2014
						}
					}
				}
				#else { delete $taxOTU{$tab[1]} }	### éliminer les seq qui n'ont pas été affiliées phylogénétiquement   ##### modif 20/8/2014
				else { delete $taxOTU{$tab[0]} }	### éliminer les seq qui n'ont pas été affiliées phylogénétiquement   ##### modif 20/8/2014
			}
		}
		close FFF;
	}
}


##################################################################
	
my %seq_clean;					
foreach my $k (keys %taxOTU) { 
	foreach my $aff (keys %{$taxOTU{$k}}) { 	
		foreach my $e (@{$child{$k}}) {		
			if ((defined $e) and ($e ne "")) { 
				my @oo = split (/_/,$e);			
# 				foreach my $s (%sample) { 20/04/2016 $sample ne contien plus les environnements qui sont sans %so
				foreach my $s (%so){	#Version 5 20/04/2016 DD
					if ($oo[0] eq $s){ 										
						push (@{$seq_clean{$aff}{$s}}, $e) ;
					}
				}	
			}	
		}
	}
}

%child = (); 	

													## VERSION 3.01 ##

############# Version 5 - 20/04/2016-  Normalisation automatique

my $abundance_norm=1000000; # valeur maximale de la normalisation
my $abundance_min=1000; #Valeur minimale
my $abundance; #Abondance courante
my %min_norm;

if ($aff_norm==0) # Pas de valeurs définies par l'utilisateur
{
  foreach my $aff (keys %seq_clean) {
	foreach my $ss(keys %{$seq_clean{$aff}}) {
	 $abundance = scalar (@{$seq_clean{$aff}{$ss}});
		if ($abundance < $abundance_norm && $abundance >=  $abundance_min) {
			$abundance_norm=$abundance;
		}
    }
  }

$min_norm{'NN'} = $abundance_norm;
$min_norm{'LCA'} = $abundance_norm;

}
else
{
$min_norm{'NN'} = $aff_norm;
$min_norm{'LCA'} = $aff_norm;

}


##################################################################



my %hash; 
# my %min_norm; 20/4/2016
my %min_norm_env;

# $min_norm{'NN'} = $aff_norm; 20/04/2016
# $min_norm{'LCA'} = $aff_norm; 20/04/2016


foreach my $aff (keys %seq_clean) {
	foreach my $ss(keys %{$seq_clean{$aff}}) {											# $ss = env 
		my $m = scalar (@{$seq_clean{$aff}{$ss}});
		if ($m >= $aff_norm) {
			$min_norm_env{$ss}=1;
		}
		foreach my $g (@{$seq_clean{$aff}{$ss}}) { 
			$hash{$ss}{$aff}{$g} = 1;
		}
	}
}

%seq_clean = ();													## VERSION 3.01 ##


##### sélectionner $min_norm séq dans chaque sample ####
my %subsample;
foreach my $k (keys %hash) { #print "hash-k --- $k\n";										# $k = $env
	foreach my $aff (keys %{$hash{$k}}) { #print "hash-aff ---- $aff\n"; 
		my $sub=0;

		my $keys = keys(%{$hash{$k}{$aff}});		
		if ($keys >= $min_norm{$aff}){ #### éliminer les ech avec un nbr de seq < le seuil de normalisation
			foreach my $kk (keys %{$hash{$k}{$aff}}) { #print "hash-kk ---- $kk\n";
				$subsample{$k}{$aff}{$kk}=1;
				$sub++;
				if ($sub == $min_norm{$aff}) { last }
			}
		}
	}
}
#######

%hash = (); %min_norm = ();												## VERSION 3.01 ##		

my %d;
my %synt;
foreach my $k (keys %subsample) { #print "synt-k ---- $k\n";											# $k = $env
	foreach my $aff (keys %{$subsample{$k}}) { #print "synt-aff ---- $aff\n"; 
		foreach my $kk (keys %{$subsample{$k}{$aff}}) { #print "synt-kk ---- $kk\n";
			my $aa = $seed{$kk};
			$synt{$k}{$aff}{$aa}{'reads'} .= $kk.",";
			$synt{$k}{$aff}{$aa}{'som'}++ ;
			$d{$aff}{$aa}= 1;
		}
	}
}	

%subsample = (); %seed = ();

my %distr; 
foreach my $v (keys %synt) {											# $v = $env
	foreach my $aff (keys %{$synt{$v}}) {
		foreach my $k (keys %{$d{$aff}} ) { 
			if (defined $synt{$v}{$aff}{$k}{'som'}) { 
				$distr{$k}{$aff}{$v} = $synt{$v}{$aff}{$k}{'som'};
			}
			else {
				$distr{$k}{$aff}{$v} = 0;
			}
		}
	}
}

%d = (); 														## VERSION 3.01 ##		

################### OTU_distribution_tax normalisé
if ($user_file eq "pooledSamples") { 
	open (ODNN, ">".$NGS_id_Results."/OTU_distribution_tax_normalized_NN.txt");
	open (ODNL, ">".$NGS_id_Results."/OTU_distribution_tax_normalized_LCA.txt");
	print ODNN "OTU_Seed\t"; print ODNL "OTU_Seed\t";
 
	foreach my $pp ( (sort {$so{$a} <=> $so{$b}} keys %synt)) { 								## VERSION 3.02 ##
		print ODNN "$pp\t" ;
		print ODNL "$pp\t" ;
	}
	print ODNN "NN taxonomy\n"; 
	print ODNL "LCA taxonomy\n"; 

	#### Print data
	foreach my $k (keys %distr) {
		print ODNL "$k\t"; print ODNN "$k\t";
		foreach my $p ( (sort {$so{$a} <=> $so{$b}} keys %synt)) { 								
			if (exists ($synt{$p})) {
			print ODNL "$distr{$k}{'LCA'}{$p}\t";
			print ODNN "$distr{$k}{'NN'}{$p}\t";
			}
			else {
			print ODNL "NA\t" ;
			print ODNN "NA\t";
			}
		}
		print ODNL "$taxOTU{$k}{'LCA'}\n";
		print ODNN "$taxOTU{$k}{'NN'}\n";
	}

	close ODNN;
	close ODNL;
}

%distr = ();  																


################################################ hashs pour calcul des indices et printing

print "\nComputing Alpha- and Beta- diversities:\n";
print LOG "\nComputing Alpha- and Beta- diversities:\n";

###################################################################################### BOUCLER sur les $env (%sample)


foreach my $env ( (sort {$so{$a} <=> $so{$b}} keys %so)) { #### modif 16/9/2014

print "Sample: $env\n";
print LOG "Sample: $env\n";

#foreach my $env ( (sort {$barIn{$a} <=> $barIn{$b}} keys %barIn)) { #### modif 16/9/2014								
	chomp($env);															
	
	my %res;
	my %res_norm;															
	#foreach my $env (keys %synt) {													
		foreach my $aff (%{$synt{$env}}) {
			foreach my $Id_OTU (%{$synt{$env}{$aff}}) {

				my @bqs = split (/;/, $taxOTU{$Id_OTU}{$aff});

				if ($#bqs == 0) {
					$res_norm{$aff}{$env}{$bqs[0]}{'nbrSEQ'}+= $synt{$env}{$aff}{$Id_OTU}{'som'};
					$res_norm{$aff}{$env}{$bqs[0]}{'OTU'}{$Id_OTU} = $synt{$env}{$aff}{$Id_OTU}{'som'}	;
				
					if ($synt{$env}{$aff}{$Id_OTU}{'som'} != 0) {
						$res_norm{$aff}{$env}{$bqs[0]}{'nbrOTU'}++;
					}
					if ($synt{$env}{$aff}{$Id_OTU}{'som'} == 1) {
						$res_norm{$aff}{$env}{$bqs[0]}{'1'}++;
					}
					if ($synt{$env}{$aff}{$Id_OTU}{'som'} == 2) {
						$res_norm{$aff}{$env}{$bqs[0]}{'2'}++;
					}
				}	

				else { 
					$res_norm{$aff}{$env}{$bqs[1]}{$bqs[2]}{'nbrSEQ'}+= $synt{$env}{$aff}{$Id_OTU}{'som'};
					$res_norm{$aff}{$env}{$bqs[1]}{$bqs[2]}{'OTU'}{$Id_OTU} =$synt{$env}{$aff}{$Id_OTU}{'som'};

					$res_norm{$aff}{$env}{$bqs[1]}{'nbrSEQ'}+= $synt{$env}{$aff}{$Id_OTU}{'som'};
					$res_norm{$aff}{$env}{$bqs[1]}{'OTU'}{$Id_OTU}= $synt{$env}{$aff}{$Id_OTU}{'som'};

					###
					if ($synt{$env}{$aff}{$Id_OTU}{'som'} != 0) {						
						$res_norm{$aff}{$env}{$bqs[1]}{$bqs[2]}{'nbrOTU'}++;
						$res_norm{$aff}{$env}{$bqs[1]}{'nbrOTU'}++;
					}

					if ($synt{$env}{$aff}{$Id_OTU}{'som'} == 1) {					
						$res_norm{$aff}{$env}{$bqs[1]}{$bqs[2]}{'1'}++;
						$res_norm{$aff}{$env}{$bqs[1]}{'1'}++				
					}
					if ($synt{$env}{$aff}{$Id_OTU}{'som'} == 2) {
						$res_norm{$aff}{$env}{$bqs[1]}{$bqs[2]}{'2'}++;
						$res_norm{$aff}{$env}{$bqs[1]}{'2'}++
					}
				}
			}
		}
	#}

########################### no normalizing

	foreach my $l (keys %taxOTU) { 
		my $nbrSeq; 
		if ($user_file eq "eachSample") {
			$nbrSeq = $nbrSeq{$l}{$env}{'each'};
		}
		if ($user_file eq "pooledSamples") { 
			$nbrSeq = $nbrSeq{$l}{$env}{'pooled'};
		}

		#print "$l --- nbrSeq --- $nbrSeq\n";
		################## parsing method ############## 

		foreach my $aff (keys %{$taxOTU{$l}} ) { 

			my @aqs = split (/;/, $taxOTU{$l}{$aff});
	
			if ($#aqs == 0 ) {
				$res{$aff}{$env}{$aqs[0]}{'nbrSEQ'}+= $nbrSeq;
				$res{$aff}{$env}{$aqs[0]}{'OTU'}{$l} = $nbrSeq	;
				if ($nbrSeq != 0) {
					$res{$aff}{$env}{$aqs[0]}{'nbrOTU'}++;
				}
				if ($nbrSeq == 1) {
					$res{$aff}{$env}{$aqs[0]}{'1'}++;
				}
				if ($nbrSeq == 2) {
					$res{$aff}{$env}{$aqs[0]}{'2'}++;
				}
			}	
			else { 
				$res{$aff}{$env}{$aqs[1]}{$aqs[2]}{'nbrSEQ'}+= $nbrSeq;
				$res{$aff}{$env}{$aqs[1]}{$aqs[2]}{'OTU'}{$l} = $nbrSeq;
				###
				$res{$aff}{$env}{$aqs[1]}{'nbrSEQ'}+= $nbrSeq;
				$res{$aff}{$env}{$aqs[1]}{'OTU'}{$l}= $nbrSeq;
				###
				if ($nbrSeq != 0) {
					$res{$aff}{$env}{$aqs[1]}{$aqs[2]}{'nbrOTU'}++;
					$res{$aff}{$env}{$aqs[1]}{'nbrOTU'}++;
				}
				if ($nbrSeq == 1) {
					$res{$aff}{$env}{$aqs[1]}{$aqs[2]}{'1'}++;
					$res{$aff}{$env}{$aqs[1]}{'1'}++;
				}
				if ($nbrSeq == 2) {
					$res{$aff}{$env}{$aqs[1]}{$aqs[2]}{'2'}++;
					$res{$aff}{$env}{$aqs[1]}{'2'}++;
				}
			}
		}
	}


	######################################################### pas de normalisation ##################
	###################################### calcul 
	foreach my $parsing (keys %res) {
		foreach my $ss (keys %{$res{$parsing}}) { 
			foreach my $a (keys %{$res{$parsing}{$ss}}) { 
		
				my $nv2_n1; my $nv2_n2;

				################## calcul de schao
				
				$nv2_n1 = $res{$parsing}{$ss}{$a}{'1'};  $nv2_n2 = $res{$parsing}{$ss}{$a}{'2'};
			
				if ((($nv2_n1 >0) and ($nv2_n2 >= 0)) or (($nv2_n1==0 ) and ($nv2_n2 == 0))) {
					$res{$parsing}{$ss}{$a}{'chao'} =  $res{$parsing}{$ss}{$a}{'nbrOTU'} + (($nv2_n1*($nv2_n1 - 1))/(2*($nv2_n2+1)));
					$res{$parsing}{$ss}{$a}{'chao'} = sprintf("%.2f", $res{$parsing}{$ss}{$a}{'chao'});
				}
				else {
					$res{$parsing}{$ss}{$a}{'chao'} = $res{$parsing}{$ss}{$a}{'nbrOTU'} + (($nv2_n1*$nv2_n1)/(2*$nv2_n2));
					$res{$parsing}{$ss}{$a}{'chao'} = sprintf("%.2f", $res{$parsing}{$ss}{$a}{'chao'});
				}
		
				############### calcul de Coverage
				my $cov;
				if ($res{$parsing}{$ss}{$a}{'nbrSEQ'} != 0 ) {	
					$cov = 1- ($nv2_n1/$res{$parsing}{$ss}{$a}{'nbrSEQ'});
					$cov = $cov * 100;
					$res{$parsing}{$ss}{$a}{'cov'} =  sprintf("%.2f", $cov);
				}
				else {
					$res{$parsing}{$ss}{$a}{'cov'} =0
				}
		
			########################################"nv3

				foreach my $aa (keys %{$leveltax{$parsing}{$a}}) { 
					my $nv3_n1; my $nv3_n2;
		
				    ##################### calcul de schao
					$nv3_n1 = $res{$parsing}{$ss}{$a}{$aa}{'1'};  
					$nv3_n2 = $res{$parsing}{$ss}{$a}{$aa}{'2'};
									
					if ((($nv3_n1 >0) and ($nv3_n2 >= 0)) or (($nv3_n1==0 ) and ($nv3_n2 == 0))) {
						$res{$parsing}{$ss}{$a}{$aa}{'chao'} =  $res{$parsing}{$ss}{$a}{$aa}{'nbrOTU'} + (($nv3_n1*($nv3_n1 - 1))/(2*($nv3_n2+1)));
						$res{$parsing}{$ss}{$a}{$aa}{'chao'} = sprintf("%.2f", $res{$parsing}{$ss}{$a}{$aa}{'chao'});
					}
					else {
						$res{$parsing}{$ss}{$a}{$aa}{'chao'} = $res{$parsing}{$ss}{$a}{$aa}{'nbrOTU'} + (($nv3_n1*$nv3_n1)/(2*$nv3_n2));
						$res{$parsing}{$ss}{$a}{$aa}{'chao'} = sprintf("%.2f", $res{$parsing}{$ss}{$a}{$aa}{'chao'});
					}
			
				    ################## calcul de Coverage
					my $cov;
					if ($res{$parsing}{$ss}{$a}{$aa}{'nbrSEQ'} != 0 ) {	
						$cov = 1- ($nv3_n1/$res{$parsing}{$ss}{$a}{$aa}{'nbrSEQ'});
						$cov = $cov * 100;
						$res{$parsing}{$ss}{$a}{$aa}{'cov'} =  sprintf("%.2f", $cov);
					}
					else {
						$res{$parsing}{$ss}{$a}{$aa}{'cov'} =0;
					}
				}
			}
		}
	}

	###################################### printing # 26/08/2016 -> dossier $path_tmpdiv
	foreach my $parsing (sort keys %leveltax) { 												# $parsing = LCA ou NN
		
		open (TAX, ">".$path_tmpdiv."/grp_tax".$parsing."_".$env."_tmp_all");
		print TAX "\n\t";
		
		foreach my $pp (sort keys %{$leveltax{$parsing}}) {
			if (($pp ne "") and ($pp ne ";")) {	
				print TAX "\n$pp\t";
				foreach my $ppp (sort keys %{$leveltax{$parsing}{$pp}}) {
					if (($ppp ne "") and ($ppp ne ";") and ($pp ne "Metazoa") and ($pp ne "Other")) {
						print TAX "\n\t$ppp";
					}
				}
			}
		}
		
		close TAX;

		my $o = $env;
				
		open (OO, ">".$path_tmpdiv."/res_".$o."_".$parsing."_tmp_all");
		print OO "\t\t\t\t\t$o\n";
		print OO "\t#seq\t#OTUs\tSchao1\tShannon\tCoverage\n";
						
		foreach my $j (sort keys %{$leveltax{$parsing}}) { 
			if ((defined $j) and ($j ne "") and ($j ne ";")) {
				my $t2; my $sh=0; 
				if (exists $res{$parsing}{$o}{$j}{'nbrOTU'}) {
					foreach my $e (keys %{$res{$parsing}{$o}{$j}{'OTU'}}) {
						if ($res{$parsing}{$o}{$j}{'OTU'}{$e} != 0) {
							$t2=$res{$parsing}{$o}{$j}{'OTU'}{$e}/$res{$parsing}{$o}{$j}{'nbrSEQ'};
							$t2 = $t2 * log($t2);
							$sh += $t2;
						}
					}
				
					$sh = $sh*(-1);
					$res{$parsing}{$o}{$j}{'shannon'} = sprintf("%.2f", $sh);
					print OO "\t$res{$parsing}{$o}{$j}{'nbrSEQ'}\t$res{$parsing}{$o}{$j}{'nbrOTU'}\t$res{$parsing}{$o}{$j}{'chao'}\t$res{$parsing}{$o}{$j}{'shannon'}\t$res{$parsing}{$o}{$j}{'cov'}\n";
				}
							
				else {
					print OO "\t0\t0\t0\t0\t0\n";
				}
			}
			foreach my $jj (sort keys %{$leveltax{$parsing}{$j}}) { 
				my $sh=0; my $t2;
				if (exists $res{$parsing}{$o}{$j}{$jj}{'nbrOTU'}) {
					foreach my $e (keys %{$res{$parsing}{$o}{$j}{$jj}{'OTU'}}) {
						if ($res{$parsing}{$o}{$j}{$jj}{'OTU'}{$e} != 0) {
							$t2=$res{$parsing}{$o}{$j}{$jj}{'OTU'}{$e}/$res{$parsing}{$o}{$j}{$jj}{'nbrSEQ'};
							$t2 = $t2 * log($t2);
							$sh += $t2;
						}
					}
									
					$sh = $sh*(-1);
					$res{$parsing}{$o}{$j}{$jj}{'shannon'} = sprintf("%.2f", $sh);
		
					print OO "\t$res{$parsing}{$o}{$j}{$jj}{'nbrSEQ'}\t$res{$parsing}{$o}{$j}{$jj}{'nbrOTU'}\t$res{$parsing}{$o}{$j}{$jj}{'chao'}\t$res{$parsing}{$o}{$j}{$jj}{'shannon'}\t$res{$parsing}{$o}{$j}{$jj}{'cov'}\n";
				}
				else {
					print OO "\t0\t0\t0\t0\t0\n";
				}
			}	
		}
		close OO;
		`paste $path_tmpdiv"/grp_tax"$parsing"_"$o"_tmp_all" $path_tmpdiv"/res_"$o"_"$parsing"_tmp_all" > $path_tmpdiv/taxonomic_distribution_"$o"_"$parsing".txt`;

# 		`rm $path_tmpdiv/*_tmp_*`;
	}
	
	%res = (); 																	


	######################################################### normalisation ##################
	###################################### calcul 
	foreach my $parsing (keys %res_norm) { 
		foreach my $ss (keys %{$res_norm{$parsing}}) { 
			foreach my $a (keys %{$res_norm{$parsing}{$ss}}) { 
		
				my $nv2_n1; my $nv2_n2;

				################## calcul de schao
				
				$nv2_n1 = $res_norm{$parsing}{$ss}{$a}{'1'};  $nv2_n2 = $res_norm{$parsing}{$ss}{$a}{'2'};
			
				if ((($nv2_n1 >0) and ($nv2_n2 >= 0)) or (($nv2_n1==0 ) and ($nv2_n2 == 0))) {
					$res_norm{$parsing}{$ss}{$a}{'chao'} =  $res_norm{$parsing}{$ss}{$a}{'nbrOTU'} + (($nv2_n1*($nv2_n1 - 1))/(2*($nv2_n2+1)));
					$res_norm{$parsing}{$ss}{$a}{'chao'} = sprintf("%.2f", $res_norm{$parsing}{$ss}{$a}{'chao'});
				}
				else {
					$res_norm{$parsing}{$ss}{$a}{'chao'} = $res_norm{$parsing}{$ss}{$a}{'nbrOTU'} + (($nv2_n1*$nv2_n1)/(2*$nv2_n2));
					$res_norm{$parsing}{$ss}{$a}{'chao'} = sprintf("%.2f", $res_norm{$parsing}{$ss}{$a}{'chao'});
				}
		
				############### calcul de Coverage
				my $cov;
				if ($res_norm{$parsing}{$ss}{$a}{'nbrSEQ'} != 0 ) {	
					$cov = 1- ($nv2_n1/$res_norm{$parsing}{$ss}{$a}{'nbrSEQ'});
					$cov = $cov * 100;
					$res_norm{$parsing}{$ss}{$a}{'cov'} =  sprintf("%.2f", $cov);
				}
				else {
					$res_norm{$parsing}{$ss}{$a}{'cov'} =0
				}
		
			########################################"nv3

				foreach my $aa (keys %{$leveltax{$parsing}{$a}}) { 
					my $nv3_n1; my $nv3_n2;
		
				    ##################### calcul de schao
					$nv3_n1 = $res_norm{$parsing}{$ss}{$a}{$aa}{'1'};  
					$nv3_n2 = $res_norm{$parsing}{$ss}{$a}{$aa}{'2'};
									
					if ((($nv3_n1 >0) and ($nv3_n2 >= 0)) or (($nv3_n1==0 ) and ($nv3_n2 == 0))) {
						$res_norm{$parsing}{$ss}{$a}{$aa}{'chao'} =  $res_norm{$parsing}{$ss}{$a}{$aa}{'nbrOTU'} + (($nv3_n1*($nv3_n1 - 1))/(2*($nv3_n2+1)));
						$res_norm{$parsing}{$ss}{$a}{$aa}{'chao'} = sprintf("%.2f", $res_norm{$parsing}{$ss}{$a}{$aa}{'chao'});
					}
					else {
						$res_norm{$parsing}{$ss}{$a}{$aa}{'chao'} = $res_norm{$parsing}{$ss}{$a}{$aa}{'nbrOTU'} + (($nv3_n1*$nv3_n1)/(2*$nv3_n2));
						$res_norm{$parsing}{$ss}{$a}{$aa}{'chao'} = sprintf("%.2f", $res_norm{$parsing}{$ss}{$a}{$aa}{'chao'});
					}
			
				    ################## calcul de Coverage
					my $cov;
					if ($res_norm{$parsing}{$ss}{$a}{$aa}{'nbrSEQ'} != 0 ) {	
						$cov = 1- ($nv3_n1/$res_norm{$parsing}{$ss}{$a}{$aa}{'nbrSEQ'});
						$cov = $cov * 100;
						$res_norm{$parsing}{$ss}{$a}{$aa}{'cov'} =  sprintf("%.2f", $cov);
					}
					else {
						$res_norm{$parsing}{$ss}{$a}{$aa}{'cov'} =0
					}
				}
			}
		}
	}

	###################################### printing
	foreach my $parsing (sort keys %leveltax) {
		open (TAX, ">".$path_tmpdiv."/grp_tax".$parsing."_".$env."_norm_tmp_all");
		print TAX"\n\t";
		
		foreach my $pp (sort keys %{$leveltax{$parsing}}) {
			if (($pp ne "") and ($pp ne ";")) {	
				print TAX "\n$pp\t";
				foreach my $ppp (sort keys %{$leveltax{$parsing}{$pp}}) {
					if (($ppp ne "") and ($ppp ne ";") and ($pp ne "Metazoa") and ($pp ne "Other")) {
						print TAX "\n\t$ppp";
					}
				}
			}
		}
		
		close TAX;
		
		my $o = $env;
	
		open (OO, ">".$path_tmpdiv."/res_".$env."_".$parsing."_norm_tmp_all");
		print OO "\t\t\t\t\t$o\n";
		print OO "\t#seq\t#OTUs\tSchao1\tShannon\tCoverage\n";
						
		foreach my $j (sort keys %{$leveltax{$parsing}}) { 
			if ((defined $j) and ($j ne "") and ($j ne ";")) {
				my $t2; my $sh=0; 
				######
				if (exists $res_norm{$parsing}{$o}) {
				#####
					if (exists $res_norm{$parsing}{$o}{$j}{'nbrOTU'}) {
						foreach my $e (keys %{$res_norm{$parsing}{$o}{$j}{'OTU'}}) {
							if ($res_norm{$parsing}{$o}{$j}{'OTU'}{$e} != 0) {
								$t2=$res_norm{$parsing}{$o}{$j}{'OTU'}{$e}/$res_norm{$parsing}{$o}{$j}{'nbrSEQ'};
								$t2 = $t2 * log($t2);
								$sh += $t2;
							}
						}
	
						$sh = $sh*(-1);
						$res_norm{$parsing}{$o}{$j}{'shannon'} = sprintf("%.2f", $sh);
						print OO "\t$res_norm{$parsing}{$o}{$j}{'nbrSEQ'}\t$res_norm{$parsing}{$o}{$j}{'nbrOTU'}\t$res_norm{$parsing}{$o}{$j}{'chao'}\t$res_norm{$parsing}{$o}{$j}{'shannon'}\t$res_norm{$parsing}{$o}{$j}{'cov'}\n";
					}
						
					else {
						print OO "\t0\t0\t0\t0\t0\n";
					}
				}
				else {
					print OO "\tNA\tNA\tNA\tNA\tNA\n";
				}


			######
			}
			foreach my $jj (sort keys %{$leveltax{$parsing}{$j}}) { 
				my $sh=0; my $t2;
				######
				if (exists $res_norm{$parsing}{$o}) {
				#####
					if (exists $res_norm{$parsing}{$o}{$j}{$jj}{'nbrOTU'}) {
						foreach my $e (keys %{$res_norm{$parsing}{$o}{$j}{$jj}{'OTU'}}) {
							if ($res_norm{$parsing}{$o}{$j}{$jj}{'OTU'}{$e} != 0) {
								$t2=$res_norm{$parsing}{$o}{$j}{$jj}{'OTU'}{$e}/$res_norm{$parsing}{$o}{$j}{$jj}{'nbrSEQ'};
								$t2 = $t2 * log($t2);
								$sh += $t2;
							}
						}
										
						$sh = $sh*(-1);
						$res_norm{$parsing}{$o}{$j}{$jj}{'shannon'} = sprintf("%.2f", $sh);
		
						print OO "\t$res_norm{$parsing}{$o}{$j}{$jj}{'nbrSEQ'}\t$res_norm{$parsing}{$o}{$j}{$jj}{'nbrOTU'}\t$res_norm{$parsing}{$o}{$j}{$jj}{'chao'}\t$res_norm{$parsing}{$o}{$j}{$jj}{'shannon'}\t$res_norm{$parsing}{$o}{$j}{$jj}{'cov'}\n";
					}
					else {
						print OO "\t0\t0\t0\t0\t0\n";
					}
				}
				else {
					print OO "\tNA\tNA\tNA\tNA\tNA\n";
				}
			}	
		}

		close OO;
		`paste $path_tmpdiv"/grp_tax"$parsing"_"$o"_norm_tmp_all" $path_tmpdiv"/res_"$o"_"$parsing"_norm_tmp_all" > $path_tmpdiv/taxonomic_distribution_"$o"_"$parsing"_norm.txt`;

# 		`rm $path_tmpdiv/*_tmp_*`;
		
	}
	 
	###############################################################

	#foreach my $k (keys %NbrSeqClean) {
	#	print "$k --- $NbrSeqClean{$k}\n"
	#}

}

# Fin calculs par environnement

my @normtypes = ("", "_norm");
my @parsingmeth = ("LCA", "NN");

my $head = `head -n 1 $NGS_id_Results/OTU_distribution_tax.txt`;
my @samples = split("\t", $head);
shift(@samples); shift(@samples);
pop(@samples); pop(@samples); 
chomp(@samples);

foreach my $parsing (@parsingmeth) {
	foreach my $ntype (@normtypes) {

		my $cpt = 1;
		
		foreach my $env (@samples) {	
			my $coden = "";
			if ($cpt<10) {
				$coden = "0";
			} 
			#print $env."\n";
			if (!( -e "tax$parsing$ntype"."_tmp0" ) )  {
				`awk -F\"\t\" '{print \$1"\t"\$2;}' $path_tmpdiv/taxonomic_distribution_$env"_"$parsing$ntype.txt > $path_tmpdiv/tax$parsing$ntype"_"tmp0`;
			}
			`awk -F\"\t\" '{print \$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8;}' $path_tmpdiv/taxonomic_distribution_$env"_"$parsing$ntype.txt > $path_tmpdiv/tax$parsing$ntype"_tmp"$coden$cpt`;
# 			`rm $NGS_id_Results/taxonomic_distribution_$env"_"$parsing$ntype.txt`; 
			$cpt++;
		}
		`paste $path_tmpdiv/tax$parsing$ntype"_"tmp* > $NGS_id_Results/taxonomic_distribution_$parsing$ntype.txt`;
# 		`rm $path_tmpdiv/tax$parsing$ntype"_"tmp*`;
	}
}




#################################################################### Annotation similitude ################## 20/05/2016
#20/5/2016

@best=&best_hit_similarity($NGS_id_Results, $NGS_id_Results."/".$panam_output."/Similarity_annotation/best_hit_uc"); # Renvoi un hash avec la taxo LCA si identité identique pour ajouter au fichier OTU_distribution_tax.txt

while( my ($d,$v) = each($best[0])){

chomp($d); chomp($v);
$best_lca_similarity{$d}=$v;

}
#20/5/2016
##########################################################################################################################################


