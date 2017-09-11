#!/usr/bin/perl
use strict;
# use warnings;


sub clean_alignments 
{


#############################
# Le 10/01/2016 DD
# Mothur n'accepte pas des noms de directory avec carectères spéciaux (-)->  Script plus lent
# 26/08/2016 -> modification du calcul du nb de positions Rque : le saut de ligne pour certains "-" sous kate et gedit n'est pas visualisé avec un more (?)
#############################



my $path;
my $fichier;
my $align;
my $clean;

($path, $fichier)=@_;

# $align="$path$fichier".".fasta";
# $clean="$path"."tmp.fasta";


# Variables ############################

my $alignement="$path$fichier".".fasta";
my $resultat="$path$fichier".".filter.fasta";
my $filter="$path$fichier".".filter";

my %position; #caractère pour chaque position 
my %filter; # 1 ou 0
my %sequence;
my @sequence; # séquence parcourue
my $soft=0.5;
my $nb_position=0;

my $i=0; my $l; my $name; my $v;
my $ct=0; # Nombre de séquences

######################################### 

# `cp $alignement "$path$fichier"."old.fasta"`; # pour garder l'ancien résultat

open (AL, $alignement);
open (R, ">$resultat");
open (F, ">$filter");

while($l=<AL>)
{
chomp($l);
  if($l =~ />(.*)/)
  {
  $name=$1;
  $ct++;
  }
  else
  {
  $sequence{$name}.=$l;
  }
}


foreach $name (keys %sequence)
{
@sequence = split("", $sequence{$name});

  for($i=0; $i<=$#sequence; $i++)
  {
    if(uc($sequence[$i]) eq "T" or $sequence{$i} eq "U"){
    $position{$i}{'T'}++; 
    }
    elsif(uc($sequence[$i]) eq "A"){
    $position{$i}{'A'}++;
    }    
    elsif(uc($sequence[$i]) eq "G"){
    $position{$i}{'G'}++;
    }
    elsif(uc($sequence[$i]) eq "C"){
    $position{$i}{'C'}++;
    }
    else{
    $position{$i}{'other'}++;
    }

  }


}


foreach $v (keys %position){

if($position{$v}{'other'} == $ct){
    $filter{$v}=0;
}
else{

  if ($position{$v}{'A'}/$ct >= $soft)
  {
   $filter{$v}=1; 
  }
  elsif ($position{$v}{'T'}/$ct >= $soft)
  {
   $filter{$v}=1; 
  }
  elsif ($position{$v}{'G'}/$ct >= $soft)
  {
   $filter{$v}=1; 
  }
 elsif ($position{$v}{'C'}/$ct >= $soft)
  {
   $filter{$v}=1; 
  }
  else # ($position{$v}{'other'}/$ct >= $soft)
  {
  $filter{$v}=0;
  }
  
}

print F $filter{$v};
if ($filter{$v}==1){  $nb_position++;}


}


foreach $name (keys %sequence)
{
print R ">$name\n";
@sequence = split("", $sequence{$name});


  for($i=0; $i<=$#sequence; $i++)
  {
    if ($filter{$i}==1){
    print R $sequence[$i];

    }
  }
print R "\n";
}


`mv "$path$fichier"."filter.fasta" $alignement`; # 

return($nb_position);
}



sub clean_alignments_mothur 
{

my $path;
my $fichier;
my $align;
my $clean;

($path, $fichier)=@_;

$align="$path$fichier".".fasta";
$clean="$path"."tmp.fasta";

# `sed "s/\./-/g" $align > $clean`;

`mothur "#filter.seqs(fasta=$align, soft=50, vertical=T)"`;


`mv "$path$fichier"."filter.fasta" $align`; # Ou dans un dossier spécial ?
# `rm $clean`;

}
1;


########################################################################################

sub root_life { # prok et eukaryotes

my $NGS_id_Results;
my $panam_output;
my $kk;
my $fff;

($NGS_id_Results, $panam_output, $kk, $fff)=@_;

open (FF , $NGS_id_Results."/".$panam_output."/Phylogeny/".$kk.".".$fff.".newick") ;

print "\nRooting  ".$kk.".".$fff.".newick\n";
	
my $tree="";
my $id_externe = "9999";
while (<FF>) {
	$tree.=$_;
}
		
my $n_tree;
if ($tree=~ /(^\(.*\(AX409584:.*?,DD141196:.*?\))(.*?)(:.*)/) {
				$n_tree = "$1".$id_externe."$3";
}

elsif ($tree=~ /(^\(.*\(DD141196:.*?,\(AX409584:.*?\).*?\))(.*?)(:.*)/){
				$n_tree = "$1".$id_externe."$3";
}

elsif ($tree=~ /(^\(.*\(DD141196:.*?,AX409584:.*?\))(.*?)(:.*)/) {
				$n_tree = "$1".$id_externe."$3";
}

elsif ($tree=~ /(^\(.*\(AX409584:.*?,\(DD141196:.*?\).*?\))(.*?)(:.*)/) {
				$n_tree = "$1".$id_externe."$3";
}	

else { 
				$n_tree = $tree;
				$id_externe = "DD141196"
}

open (FFF,">".$NGS_id_Results."/".$panam_output."/n_tree");
print FFF "$n_tree";
close FFF;
	
my $data = $NGS_id_Results."/".$panam_output."/n_tree";
my $res = $NGS_id_Results."/".$panam_output."/Phylogeny/".$kk.".".$fff."_rooted.newick";
	
	
my $in = Bio::TreeIO->new(-format => 'newick', -file => $data);
	
my $out = Bio::TreeIO->new(-format => 'newick', -file => ">$res");
	
while( my $t = $in->next_tree ){
  my ($a) = $t->find_node(-id =>"9999");
  $t->reroot($a);
  $out->write_tree($t);
}


}

################################################################### Root enterovirus

sub root_enterovirus { 

my $NGS_id_Results;
my $panam_output;
my $kk;
my $fff;

($NGS_id_Results, $panam_output, $kk, $fff)=@_;

open (FF , $NGS_id_Results."/".$panam_output."/Phylogeny/".$kk.".".$fff.".newick") ;

print "\nRooting enterovirus ".$kk.".".$fff.".newick with AF326765/AF406813\n";
	
my $tree="";
my $id_externe = "9999";
while (<FF>) {
	$tree.=$_;
}
		
my $n_tree;

if ($tree=~ /(^\(.*\(AF326765:.*?,AF406813:.*?\))(.*?)(:.*)/) {
				$n_tree = "$1".$id_externe."$3";
			}

			elsif ($tree=~ /(^\(.*\(AF406813:.*?,\(AF326765:.*?\).*?\))(.*?)(:.*)/){
				$n_tree = "$1".$id_externe."$3";
			}

			elsif ($tree=~ /(^\(.*\(AF406813:.*?,AF326765:.*?\))(.*?)(:.*)/) {
				$n_tree = "$1".$id_externe."$3";
			}

			elsif ($tree=~ /(^\(.*\(AF326765:.*?,\(AF406813:.*?\).*?\))(.*?)(:.*)/) {
				$n_tree = "$1".$id_externe."$3";
			}	

			else { 
				$n_tree = $tree;
				$id_externe = "AF406813"
			}


open (FFF,">".$NGS_id_Results."/".$panam_output."/n_tree");
print FFF "$n_tree";
close FFF;
	
my $data = $NGS_id_Results."/".$panam_output."/n_tree";
my $res = $NGS_id_Results."/".$panam_output."/Phylogeny/".$kk.".".$fff."_rooted.newick";
	
	
my $in = Bio::TreeIO->new(-format => 'newick', -file => $data);
	
my $out = Bio::TreeIO->new(-format => 'newick', -file => ">$res");
	
while( my $t = $in->next_tree ){
  my ($a) = $t->find_node(-id =>"9999");
  $t->reroot($a);
  $out->write_tree($t);
}


}

########################################### Similitude
# TAXONOMIE sur la base du meilleur hit - LCA sur meilleur hit si identité égale - Table m8 du usearch - alignement global
# 20/05/2016 vsearch remplace usearch -> pas de base en formatage "db" disponible


sub best_hit_similarity {

my $NGS;
my $blast; # Patt complet vers le fichier best_hit_uc;

($NGS, $blast)=@_;

# my $tab_OTU="/home/didier/DATA/Recherche/BoiteAoutils/panamv5quick/OTU_distribution_tax.txt";



open (BLAST,$blast ) || "Fichier non trouvé";;
# open (TABLEAU, ">OTU_distribution_identite.txt");

my $tab_OTU=$NGS."/OTU_distribution_tax.txt";
my $ct_otu=0;
my @usearch=<BLAST>;
my @ligne;
my $best_hit_courant;
my $identite_courant;
my $identite;
my %best_hit_identite;

while ($ct_otu < $#usearch)
{

my @hits_5;
# $hits_5 = ();
my $hit=0;
my $trouve="N";
my $best_hit;

  while ($trouve eq "N") # Initilalisation du tableau des meilleurs hits de même identité
  {
  chomp ($usearch[$ct_otu]);
  @ligne=split("\t", $usearch[$ct_otu]);
#     print "$ct_otu  ligne = @ligne\n";
    $best_hit_courant=$ligne[0];
    $identite_courant=$ligne[4];
# print "Courant -> $best_hit_courant - $identite_courant\n";

    if($hit==0)
    {
    $best_hit=$best_hit_courant;
    $identite= $identite_courant;
    $hits_5[$hit]=$ligne[2];
#     print "hit=$hit $best_hit\t$identite\n";
    $hit=$hit+1;
    $ct_otu=$ct_otu+1;


    }
    else
    {
      if (($best_hit eq  $best_hit_courant) and ($identite eq $identite_courant))
      {
      $hits_5[$hit]=$ligne[2];
#           print " 2 hit=$hit :$best_hit\t$identite\n";
      $hit=$hit+1;
      $ct_otu=$ct_otu+1;

      }
      elsif (($best_hit eq  $best_hit_courant) and ($identite ne $identite_courant))
      {
      $ct_otu=$ct_otu+1;

      }
      else
      {
      $trouve ="Y";
      }
   
   }
      
  } # Fin de recherche meilleur hit avec même identité


#############################################
# print $hits_5[$hit];

	


	my $ct=0;
	my $recherche=1;
	my $lca="";
	my $ct_lca=0;
# 	print "bloc : @hits_5\n\n";

	    #Recherche LCA
	    while($recherche==1)
	    {

		my $i=0;
		my $chaine_courante="";
		while($i<$hit && $recherche==1) # Parcours pour le niveau x =$ct_cla
		{
# 		print "$i-$hits_5[$i]\n";
		my @chaine=split(";", $hits_5[$i]);
# 		print "Niveau taxo = $ct_lca  Taxo = $chaine[$ct_lca]\n";
		  if(exists($chaine[$ct_lca]))
		  {
		      if($i==0)
		      {
		      $chaine_courante=$chaine[$ct_lca]; # le première chaine doit être identique aux 5 autres
		      }
		      elsif ($i > 0 && $chaine_courante ne $chaine[$ct_lca])
		      {
		       $recherche=0; # pour un niveau donné plus de congruence
		      }
		      

		  
		    }
		    else
		    {
		    $recherche=0; # fin car plus de niveau
		    }

	      

		$i=$i+1;
		}
	    if($recherche==1)
	    {
	      if($lca eq "")
	      {
	      $lca=$chaine_courante;
	      }
	      else
	      {
	      $lca=$lca.";".$chaine_courante;
	      }
	    }
	     $ct_lca=$ct_lca+1;
	   } # Fin de recherche


$best_hit_identite{$best_hit}="$lca\t$identite\n";
  	
# print "$best_hit\t$lca\t$identite\n";
}
close BLAST;

my $label;
my %otu_distribution;
my $id;

open (F1,$tab_OTU) || die "$tab_OTU not found \n" ;
my %otu_distribution;
while (my $l = <F1>)
{
	chomp ($l);
	if ($l !~ /OTU_Id.*/){
	@ligne=split("\t", $l);

	
		$otu_distribution{$ligne[0]} = $l; 
	}
	else
	{
	$label="$l\tBest_hit_identity\tIdentity (%)\n";
	}
}
close F1;

#Ancien tableau écrasé
open (F1,">".$tab_OTU);
print F1 $label;

foreach $id (keys %otu_distribution){
print F1 "$otu_distribution{$id}\t$best_hit_identite{$id}";

}

return(\%best_hit_identite);
}
