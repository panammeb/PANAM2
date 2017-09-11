#!/usr/bin/perl
use warnings;
use strict;
use FindBin;

use lib "$FindBin::RealBin/modules";
use lib "$FindBin::RealBin/bin/Text-CSV_XS-1.04/lib/x86_64-linux-gnu/perl/5.22.1/";

use Cwd 'abs_path';
use Text::CSV_XS;






my $R_scripts="$FindBin::RealBin/R";
my $path_panam="$FindBin::RealBin";



## Modification du 21/07/2015 - DD
# Récuperation des mnnd calculé avec les abondances et presence/absence - Variables renomées
# Fichiers .fig -> .svg



my ($USAGE) =  "\n\n\t****** USAGE of $0 PROGRAM ******\n\n\n\n\tUSAGE: perl $0 <panam.ini file> \n\n\n\n";
die "$USAGE" if ((scalar(@ARGV))< 1);
my $option_file = $ARGV[0];
chomp($option_file);

die "\n\n\t=> Cannot find configuration file: $option_file.\n\n" unless (-e $option_file);
die "\n\n\t=> Configuration file, $option_file, appears to be empty!\n\n" if (-z $option_file);


##### recuperation des variables

use parse_ini ;
my @parse = &parse_ini($option_file);

# my $otu_distrib = $parse[23];
# die "\n\n\tCannot find OTU_distribution_tax file: $otu_distrib. Check $option_file.\n\n" unless (-e $otu_distrib);
# die "\n\n\tOTU_distribution_tax file, $otu_distrib, appears to be empty!\n\n" if (-z $otu_distrib);
# 
# my $clades_file = $parse[24];
# die "\n\n\tCannot find PANAM_Clades_NN file: $clades_file. Check $option_file.\n\n" unless (-e $clades_file);
# die "\n\n\tPANAM_Clades_NN file, $clades_file, appears to be empty!\n\n" if (-z $clades_file);
# 
# my $phylo_dir = $parse[25];  #Fichier clade
# die "\n\n\tCannot find the phylogenies' folder: $phylo_dir. Check $option_file.\n\n" unless (-d $phylo_dir);
# die "\n\n\tThe phylogenies folder, $phylo_dir, appears to be empty!\n\n" unless ( my @files = glob("$phylo_dir/*"));


####### Domaines/Groupes
my %dom; # domaines

print "Domains studied:";
foreach my $d (keys %{$parse[32]}){
$dom{$d}=1;
print " $d ";
}
print "\n";

foreach my $d  (keys %dom) {
die "\n\n\t Domain name is not correct. Check panam.ini.\n\n" unless (($d eq "bacteria") or ($d eq "eukaryota") or ($d eq "archaea") or ($d eq "enterovirus"));
}



my $refBase = $parse[26];


if (! defined ($refBase)){

  if (exists ($dom{'enterovirus'})){
  $refBase="$FindBin::RealBin/bd_entero";
  }
  else {

  $refBase="$FindBin::RealBin/bd_ssrna";
  }

}

die "\n\n\tCannot find the reference base path at: $refBase. Check $option_file.\n\n" unless (-d $refBase);
die "\n\n\tThe outgroup file is missing in $refBase!\n\n" unless (( -e "$refBase/outgrp_profiles") and !(-z "$refBase/outgrp_profiles"));  ### modif 27/8/2014

# working folder
my $NGS_id_Results = $parse[0];
die "\n\n\tOutput directory for NGS analyses is not defined. Check $option_file.\n\n" unless ( defined $NGS_id_Results);
unless (-d $NGS_id_Results) { mkdir $NGS_id_Results || die "Could not create Directory $NGS_id_Results !\n"; }

# 15/01/2016
#Fichier clade
my $otu_distrib = $NGS_id_Results."/OTU_distribution_tax.txt";
my $clades_file = $NGS_id_Results."/PANAM_Clades_NN.txt";
my $phylo_dir = $NGS_id_Results."/panam_output/Phylogeny";

my $dir_results = $NGS_id_Results."/PhyloDiv_output";
if (-e "$NGS_id_Results/PhyloDiv_output/") {qx(rm -R "$NGS_id_Results/PhyloDiv_output/");}
unless (-d $dir_results) { mkdir $dir_results || die "Could not create Directory $dir_results !\n"; }

my $otu_distrib_tmp = $dir_results."/OTU_distribution_tmp";
my $errorlog= $NGS_id_Results."/.error.log";



my $log_file=">>". $NGS_id_Results."/panam.log";
open (LOG, $log_file);
#print LOG "Script phylodiv\n\n";

my $system;

print "\n\n* STEP 4 - Phylogenetic indices by PANAM2\n\n";
print LOG "\n\n* STEP 4 - Phylogenetic indices by PANAM2\n\n";


################## regénérer le tableau otu_distrib pour enlever tous les espaces

my $csv = Text::CSV_XS->new( { sep_char => "\t" } );
my %abund;

open(F, $otu_distrib) or die "Could not open $otu_distrib\n";
open(FF, ">".$otu_distrib_tmp);
while (my $line = <F>) {
	chomp $line;
	if ($csv->parse($line)) {
		my @fields = $csv->fields();

		##### ajout 28/8/2014
		print FF "$fields[0]\t";
		########
		
		#for (my $i=1; $i<$#fields-1; $i++) {
		for (my $i=2; $i<$#fields-1; $i++) {   ### modif 28/8/2014

			$fields[$i]=~ s/\s//g ;
			print FF"$fields[$i]\t";

			#### récupérer l'abondance de chaque otu
			if ($fields[$i] =~  m/^\d+$/){
				#$abund{$fields[1]} += $fields[$i]; 
				$abund{$fields[0]} += $fields[$i];  ### modif 28/8/2014
			}
			##############
		}
	}
	print FF"\n"
}

close F;
close FF;

###################################
### grp_ext

my $file_grpext = $refBase."/outgrp_profiles" ;	#### modif 27/8/2014
my %grp_ext;
open (GE, $file_grpext) ;
my @grp_ext = <GE>;
foreach my $e (@grp_ext) {
	my @tt = split (/\t/, $e) ;
	$tt[0]=~s/\s+//gi; $tt[2]=~s/\s+//gi;	
	$grp_ext{$tt[0]} = $tt[2]; 
}
close GE;


###################################

open (CLADE, $clades_file) or die "can not open file $clades_file";
my %clade; my $clade_id; my %profil; my $profile_id; my $bootstrap; my $tax;

while (<CLADE>) {
	if ($_=~ /> (.*?)[\s+\t](.*?)\n/) {
		$clade_id = $1; 
		$tax = $2;
	}
	else {
		#if ($_=~ /.*Phylogeny:\s*(\w*_rooted.newick).*Bootstrap:(.*)\n/) { 
		if ($_=~ /Phylogeny:\s(.*?)\tBootstrap:(.*?)\n/) {
			$profile_id = $1; 
			$profile_id=~s/^\s+//gi;
			$profile_id=~s/\s+$//gi;
			
#                             my $folder=""; # 16/07/2017 Changemment du nom du profil
#         
#                             if(($profile_id=~ /^(.*)_cp_.*\.(\d)_rooted\.newick/) || ($profile_id=~ /^(.*)_.*\.(\d)_rooted\.newick/)){
#                             $folder="$1_$2";
#                             }
#         
#                             $profile_id=$folder;# 16/07/2017 Changemment du nom du profil
			
			

			$profil{$profile_id}{'clade'}.= $clade_id.", " ;

			$bootstrap = $2;
			$clade{$profile_id}{$clade_id}{'bootstrap'} = $bootstrap;
		}

		elsif ($_=~ /.*neighbor[s]{0,1}: (.*?)\n/) {
			my $otu = $1;
			$otu=~s/\s+//gi;
			$profil{$profile_id}{'OTU'} .= $otu ;

			$clade{$profile_id}{$clade_id}{'OTUs'} = $otu;
			$clade{$profile_id}{$clade_id}{'tax'} = $tax;		
		}
	}
}
close CLADE;

#### stocker les profils avec leurs OTUs dans un/des fichiers tmp ####

#open (RES, ">".$dir_results."/Profile_Results");   ### modif 3/11/2014
open (RES, ">".$NGS_id_Results."/Profile_Results.txt");
print RES "Taxonomic Profile\tOTU number\tmean PD\tmean MPD\tmean MNND\tmean NRI\tmean NTI\tUnweighted Unifrac\n" ;

print "\nPhylogenetic indices for each phyletic groups in process... \n"; 
print LOG "\n\nPhylogenetic indices for each phyletic groups in process... \n";

foreach my $k (keys %profil) { 

	`mkdir $dir_results/$k`;
	my $k_results = $dir_results."/".$k;

	open (PR, ">".$k_results."/".$k."OTUs_tmp") ;
	my @t = split (/,/, $profil{$k}{'OTU'});
	foreach my $e (@t) {
		print PR "$e\n";
	}
	close PR;
        

        print  LOG "\t$k\n";
        my $folder="";
        
        if(($k=~ /^(.*)_cp_.*\.(\d)_rooted\.newick/) || ($k=~ /^(.*)_.*\.(\d)_rooted\.newick/)){
            $folder="$1_$2";
        }
        
        print "\t$folder\n";

	$system="R --vanilla --args ".$k." ".$otu_distrib_tmp." ".$phylo_dir."/".$k." ".$k_results."/".$k."OTUs_tmp ".$k_results."/".$k.".svg ".$k_results."/".$k.".shared_branches ".$k_results."/".$k.".unifrac ".$k_results."/".$k.".pd ".$k_results."/".$k.".mpd ".$k_results."/".$k.".mntd ".$k_results."/".$k.".nri ".$k_results."/".$k.".nti ".$k_results."/".$k.".tab_tmp ".$k_results."/".$k."_unifrac_cluster.svg ".$k_results."/".$k."_unifrac_pcoa.svg < ".$R_scripts."/phylo_01.R 2>> ".$errorlog;### 21/07/2015
	print  LOG $system."\n";
	qx($system);
	
        if (-e $k_results."/".$k.".tab_tmp"){
	`sed -i 's/-Inf/NA/g' $k_results/$k".tab_tmp" ` ;
        `sed -i 's/_rooted.newick//g' $k_results/$k".tab_tmp" ` ;
	#`cat $k_results/$k".tab_tmp" >> $dir_results"/Profile_Results"`   #### modif 3/11/2014
	
	`cat $k_results/$k".tab_tmp" >> $NGS_id_Results"/Profile_Results.txt"`    #### modif 3/11/2014
	}
	else {
	print LOG $k_results."/".$k.".tab_tmp does not exist\n";
	}
	
}

####################################"

####### Résultats par clade	##############

foreach my $k (keys %clade) {	
	my $nom_profil;
	if ($k =~ /(\w*)_(.*?)_rooted.newick/){
		$nom_profil = $1; 
	}

	foreach my $kk (keys %{$clade{$k}}) {
		######### test si la seq ref du clade appartient au grp externe
		### si la séq ref appartient au grp externe mais qu'on est dans le même profil, on garde le clade, sinon on drop
		if (exists $grp_ext{$kk}) {	
			if ($grp_ext{$kk} eq $nom_profil) {	
				my @tab = split (/,/, $clade{$k}{$kk}{'OTUs'});
				$clade{$k}{$kk}{'nbrOTUs'} = scalar(@tab) ; 

				foreach my $e (@tab) {
					$e=~s/^\s+//gi;
					$e=~s/\s+$//gi;
					$clade{$k}{$kk}{'seq'}+=$abund{$e};
					if ($abund{$e} == 1 ) {
						$clade{$k}{$kk}{'singl'}++
					}				
				}
			}
			else {
				delete $clade{$k}{$kk} ;
			}
		}

		####
		else {
			my @tab = split (/,/, $clade{$k}{$kk}{'OTUs'});
			$clade{$k}{$kk}{'nbrOTUs'} = scalar(@tab);

			foreach my $e (@tab) {
				$e=~s/^\s+//gi;
				$e=~s/\s+$//gi;
				if (exists $abund{$e}) {# 15/06/2017
                                    $clade{$k}{$kk}{'seq'}+=$abund{$e};
                                    if ($abund{$e} == 1 ) {
					$clade{$k}{$kk}{'singl'}++;
                                    }
                                    
                                }
                                else { $clade{$k}{$kk}{'seq'}=0; 
#                                 print "$k\t$kk\n"; # Absence d'abondances dans la table OTU -> Chloroplastes par exemple
                                }# 15/06/2017
			}
		}
		####
	}
}

######### deuxième boucle sur %clade pour ne pas lancer R pour tous les clades, mais uniquement sur ceux avec au moins 3 (?) otus et avec des bootstraps > 0.6 (?)

print "\nPhylogenetic indices for each clades in process... \n"; 
print LOG "\n\nPhylogenetic indices for each clades in process... \n"; 

my %print;
foreach my $k (keys %clade) {

#	my $nom_profil; my $domain ; 
#	if ($k =~ /(\w*)_(.*?)_1_rooted.newick/){
#		$nom_profil = $1; 
#		$domain = $2;
#	}

	foreach my $kk (keys %{$clade{$k}}) {
	
		if (defined $clade{$k}{$kk}{'singl'}) {}
		else { $clade{$k}{$kk}{'singl'} = 0 }

		if (($clade{$k}{$kk}{'nbrOTUs'} > 2) and $clade{$k}{$kk}{'seq'} > 0 and ($clade{$k}{$kk}{'bootstrap'} > 0)) {	# Modif du 15/07/2017 Ajour du nb des séquences > 0
			my $mntd; my $deepest;
			
# 			my $folder;
# 			if(($k=~ /^(.*)_cp_.*\.(\d)_rooted\.newick/) || ($k=~ /^(.*)_.*\.(\d)_rooted\.newick/)){
#                         $folder="$1_$2";
#                         }
                        
			
                    
		
			open (PR, ">".$dir_results."/".$k."/".$kk."_OTUs_tmp") ;
			print PR "$kk\n";
			my @t = split (/,/, $clade{$k}{$kk}{'OTUs'});
			foreach my $e (@t) {
				print PR "$e\n";
			}
			close PR;
			
			
			
			print "\t $kk\n";
			print LOG "\t $kk\n";
			$system="R --vanilla --args ".$kk." ".$phylo_dir."/".$k." ".$dir_results."/".$k."/".$kk."_OTUs_tmp ".$dir_results."/".$k."/".$kk."_clade.svg ".$otu_distrib_tmp." ".$dir_results."/".$k."/".$kk."_tmp ".$dir_results."/".$k."/".$kk.".nri < ".$R_scripts."/clade_02.R  2>> $errorlog";     ### 21/07/2015
                        qx($system);
                        print LOG $system."\n";

			open (H, $dir_results."/".$k."/".$kk."_tmp");
			my @tt = <H>;
			close H;

			my @nb_samples = split(/\s/, $tt[0]); ### 21/07/2015
			my @mnnd_ab = split(/\s/, $tt[1]); ### 21/07/2015
			my @mnnd_pa = split(/\s/, $tt[2]);  ### 21/07/2015
			my @depth_deepest = split(/\s/, $tt[3]);  ### 21/07/2015
		

			$print{$k}{$kk} = $clade{$k}{$kk}{'bootstrap'}."\t".$clade{$k}{$kk}{'nbrOTUs'}."\t".$clade{$k}{$kk}{'seq'}."\t".$clade{$k}{$kk}{'singl'}."\t".$nb_samples[1]."\t".$mnnd_ab[1]."\t".$mnnd_pa[1]."\t".$depth_deepest[1]."\t".$k."\t".$kk."\t".$clade{$k}{$kk}{'tax'};   ### 21/07/2015

		}
	}
}

#open (RR, ">".$dir_results."/Clade_Results");   ### modif 3/11/2014
open (RR, ">".$NGS_id_Results."/Clade_Results.txt");   ### modif 3/11/2014
# print RR "Clade\tBootstrap\t#OTUs\t#Sequences\t#Singletons\tMNND\tdepth_deepest\tProfile\tNN\tTaxonomy\n" ;   ### modif 3/11/2014
print RR "Clade\tBootstrap\t#OTUs\t#Sequences\t#Singletons\t#Samples\tMNND with abundances\tMNND presence/absence\tdepth_deepest\tTree\tNN\tTaxonomy\n" ;    ### 21/07/2015


open (NB, ">".$NGS_id_Results."/PhyloDiv_output/PANAM_Clades_Phylo.txt");   ##### 21/07/2015

my $i=1;
foreach my $k (keys %print) {
	foreach my $kk (keys %{$print{$k}}) {
		my $cl = "clade_".$i;
		print RR "$cl\t$print{$k}{$kk}\n";
		##### ajout 3/11/2014
		print NB ">$cl\t$clade{$k}{$kk}{'tax'}\n\t\tNearest Neighbor: $kk\tPhylogeny: $k\tBootstrap: $clade{$k}{$kk}{'bootstrap'}\n\t\t$clade{$k}{$kk}{'nbrOTUs'} neighbors: $clade{$k}{$kk}{'OTUs'}\n\n"; 
		$i++;
	}
}

close RR;
close NB; ##### ajout 3/11/2014



# `rm $otu_distrib_tmp`;

