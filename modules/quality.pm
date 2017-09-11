#!/usr/bin/perl -w
use strict;
use warnings;

####################################################################################
#DEMULTIPLEXAGE

# Script qui à partir d'un assemblage : sequences en sens + et demultiplexage
# /i pour prendre en compte les tags en minuscule
# 1/06/2015 Traitement de la dégénrescence des primers  - BUG ce n'est pas l'amplicon qui est écrit mais la séquence totale !!! tmp remplacé par BAD - tri des séquences non ATGC -> elimination des Ns (ici après assemblage)
# 2/6/2015 Introduction du traitement des lignes de commande 
# 2/6/2015 V2 : ajout des chimères et homopolymères
# 2/6/2015 V3 : bugs sur primers en reverse
# 07 5/06/2015 -> déreplication pour le traitement des chimères
# 27/07/2015 -> Elimination des polytags (JCC)
# 27/07/2015_1 : imbrication des boucles : sens + et recherche des tags. - Traitement des chimères - 
# Numérotation remplacée par un sed : effet sur PANAM ?
#Elimination des lignes vides dans le fichier tag
#29/07/2015 -> renomage des séquences dans un dossier tmp sinon les identifiants ne sont pas uniques
#10/01/2016 -> taille de l'amplicon-> module
#26/01/2016 Traitement des fichiers vsearch -> 80 colonnes
#15/2/2016 données sans tags
#3/05/2016 Garde ou non les primers (important pour pour script Enteros)
#12/05/2016 Test si fichier existe pour chimeres
#17/05/2016 use strict, kept_primers -> yes/no
#18/05/2016 -> validation sur données
# 1/06/2016 dossier BAD avec uniquement l'identifiant
# 25/08/2016 traitement des primers multiples
#12/07/2017 éimination des séquences des non mergées

sub demultiplex_miseq_exact{

my $path; #path du dossier résultat
my $tags;
my $seq_F;
my $seq_R;
my $min; 
my $max;
my $CheckChimeras;
my $kept_primers; # yes valeur par défaut
my $path_panam; # <panam>/bin/

($path, $tags,  $seq_F,  $seq_R, $min, $max, $CheckChimeras, $kept_primers,$path_panam)= @_ ;

####################
my $uclust;
my $command;
my $log_file; # Fichier log;
my $path_quality_output; my $quality_output;
my $fasta_merged;
my $fichiertag; # tag présent ou non
my $forw_regexp; my $reve_regexp;my @tab_forw; my @tab_rev;
my $comp_rev_regexp;
my $rev_regexp; my $comp_forw_regexp; 
my $revcomp;
my $amplicon;
my $size;
my $ct; my $i; my $id;
my $n; # presence de Ns
my %sequence; # hash des séquences
my $sens; # strand +/-

my $trouvechant;
my $trouvfasta;
my $nom;
my $numtrouve;
my $tag_amont; my $tag_aval; my $echantillon; my $nb_tags;
my @ligne_tag; my @tag;
####################

$uclust=$path_panam."vsearch"; #path vers l'exécutable


$log_file=">>".$path."/panam.log"; 

open (LOG, $log_file);
print LOG "Demultiplexing and cleaning sequences
Primers -> Forward: $seq_F  Reverse: $seq_R length bewteen $min and $max
Check chimeras: $CheckChimeras \n";

print  "\t\tDemultiplexing and cleaning sequences
Primers -> Forward: $seq_F  Reverse: $seq_R length bewteen $min and $max
Check chimeras: $CheckChimeras \n";



####################

#print "$min $max";

$path_quality_output=$path."/quality_output";


# if (-e $path_quality_output) {qx(rm -R $path_quality_output);}# Sinon ajout de séquences à l'existant
# # system("mkdir $path_quality_output");
# system("mkdir $path_quality_output/BAD");
# system("mkdir $path_quality_output/tmp");

open(EQ,">".$path_quality_output."/tmp/eq_seq.csv") ||  die "Non créé" ;
open (BAD, ">".$path."/quality_output/BAD/seqAll_bad_primers.fasta") ||  die " seqAll_bad_primers.fasta Non créé" ;
open (ATGC, ">".$path."/quality_output/BAD/seqAll_bad_ATGC.fasta") ||  die " seqAll_bad_ATGC.fasta Non créé" ;
# open (LENGTH, ">".$path."/quality_output/BAD/seqAll_bad_length.fasta") ||  die " seqAll_bad_length.fasta Non créé" ; # 23/06/2017 à vérifier

# open (MULTPRIMER, ">".$path."/quality_output/BAD/seqAll_multi_primers.fasta") ||  die " seqAll_multi_primers.fasta Non créé" ; 

$fasta_merged=$path."/quality_output/pairs/merged.fasta";
open(MERGE, $fasta_merged) || die "$fasta_merged not found";
# open(TAG, $tags)  || die "$tags not found";

if  ($tags ne "NULL")#15/2/2016
{
qx(sed -i '/./!d' $tags);  # Elimination des lignes vides dans le fichier tag
open(TAG, $tags);
@tag=<TAG>;
$fichiertag=1;
}
else
{
$fichiertag=0;
}

####################



@tab_forw=split(//,$seq_F);
$forw_regexp=regexp(@tab_forw);

$comp_forw_regexp=reverse($forw_regexp);
$comp_forw_regexp=~tr/ACGT[]/TGCA][/;

@tab_rev=split(//,$seq_R);
$rev_regexp=regexp(@tab_rev);

$comp_rev_regexp=reverse($rev_regexp);
$comp_rev_regexp=~tr/ACGT[]/TGCA][/;


$amplicon="$forw_regexp.*$comp_rev_regexp";
print LOG "Amplicon strand + : $amplicon\n";

# $reve_regexp=~tr/ACGURMBD/UGCAYKVH/;

sub regexp{
	my $regexp;
	foreach (@_){
		
		if ($_ eq "R"){
			$regexp.="[AG]";
		}
		elsif ($_ eq "Y"){
			$regexp.="[CT]";
		}
		elsif ($_ eq "M"){
			$regexp.="[CA]";
		}
		elsif ($_ eq "K"){
			$regexp.="[TG]";
		}
		elsif ($_ eq "W"){
			$regexp.="[TA]";
		}
		elsif ($_ eq "S"){
			$regexp.="[CG]";
		}
		elsif ($_ eq "B"){
			$regexp.="[CTG]";
		}
		elsif ($_ eq "D"){
			$regexp.="[ATG]";
		}
		elsif ($_ eq "H"){
			$regexp.="[ATC]";
		}
		elsif ($_ eq "V"){
			$regexp.="[ACG]";
		}
		elsif ($_ eq "N"){
			$regexp.="[ACGT]";
		}
		else{
			$regexp.=$_;
		}
	}
	return $regexp;
}



# Création du hash
$ct=1;
$n=0;

while(<MERGE>) # Le 26/01/2016 Traitement des fichiers vsearch -> 80 colonnes
{
chomp($_);

  if($_ =~ />(.*)/)
  {

	if($n==1)
	{
	print ATGC $nom."\n"; #12/04/2016
	delete($sequence{$nom});# Supprime les Ns...
  	}
	
	$n=0;
  	$nom="Seq$ct";
	print EQ "$1\t$nom\n"; # 12/04/2016 ajout de l'écriture dans le fichier équivalent (du 26/01 au 12/04 le fichier est vide
	$ct++;
  }
  else
  {
    if( $_ =~ /[^ATGC\s]+/i) # Supprime les Ns...
    {
   
	$n=1;
    }
    else
    {
    $sequence{$nom}.=$_;
    }

  }

}




# Dans le bon sens-> +

foreach $id (keys(%sequence))
{

$sens=0; # Primers trouvés

# Le 25/08/2016

#  if ($sequence{$id} =~ /.*$forw_regexp{1}.*$comp_rev_regexp{1}.*/i ) # sens +
#   { 
#   $sens=1;
#   }
#  elsif ($sequence{$id} =~ /.*$rev_regexp{1}.*$comp_forw_regexp{1}.*/i ) #sens -
#   {
#   $sens=1;
#   chomp($sequence{$id});
#     $revcomp = reverse($sequence{$id});
# 
# 	# complement the reversed DNA sequence
#         $revcomp =~ tr/ACGTacgt/TGCAtgca/;
#         $sequence{$id}= "$revcomp\n";
#         
#   }
#   else
#   {
# #   print BAD ">".$id."\n".$sequence{$id}."\n";
#   print BAD $id."\n";
#   delete ($sequence{$id}); # Rajouté le 27/05
#   # Absence des 2 primers
#   }

#Modif

my @countf=($sequence{$id} =~ /$forw_regexp/gi);
my @countr=($sequence{$id} =~ /$comp_rev_regexp/gi);

 if (@countf ==1 && @countr ==1) # sens +
  { 
  $sens=1;
  }
 elsif (@countf ==0 && @countr ==0) #pas de multiprimers
  {
    my @countrevf=($sequence{$id} =~ /$rev_regexp/gi);
    my @countrevr=($sequence{$id} =~ /$comp_forw_regexp/gi);
    
    if (@countrevf ==1 && @countrevr ==1) # Séquence en sens -
    {
    $sens=1;
    chomp($sequence{$id});
    $revcomp = reverse($sequence{$id});

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        $sequence{$id}= "$revcomp\n";
    }
    else
    {
    print BAD $id."\n";
    delete ($sequence{$id});
    }
        
  }
  else # Multiprimers et/ou primers absents
  {
#   print BAD ">".$id."\n".$sequence{$id}."\n";
  print BAD $id."\n";
  delete ($sequence{$id}); # Rajouté le 27/05
  # Absence des 2 primers
  }


# Fin du 25/08




#Les primers ont été trouvés

  if ($sens==1 && $fichiertag==1) # Presence du fichier tag -> TODO est-ce indispensable -> 1 échantillon ?

  {

	$trouvechant="";
	$trouvfasta="";
	$numtrouve=0;
	# Parcours de tous les échantillons
	for ($i=0; $i<=$#tag; $i++)
	{
		chomp($tag[$i]);
		@ligne_tag=split("\t", $tag[$i]);
		$echantillon=$ligne_tag[0];
		$tag_amont=$ligne_tag[1];
				
		if (exists $ligne_tag[2]) # nb_tags = 2
	  	{
	  		$nb_tags=2;
			$tag_aval=reverse($ligne_tag[2]);
			$tag_aval =~ tr/ACGTacgt/TGCAtgca/;
	
			if( $sequence{$id} =~ /.*($tag_amont)($amplicon)($tag_aval).*/i) 
			{
	    			$trouvechant=$echantillon;
				$size=length($2);
				$trouvfasta=$2;
				$numtrouve++;
				if ($numtrouve==2) # multitag
				{
					last;
				}
			}
	  	} 
		else  # nb_tags = 1
		{	
			if( $sequence{$id} =~ /.*($tag_amont)($amplicon).*/i) 
			{
	    			$trouvechant=$echantillon;
                                $size=length($2);
				$trouvfasta=$2;
				$numtrouve++;
				if ($numtrouve==2) # multitag
				{
					last;
				}
	   		 }		
		}		
 	}

# Il faut ajouter les conditions sur la taille ici
	
	if ($numtrouve==2) # multitag
	{
		open (RES, ">>".$path."/quality_output/BAD/seqAll_multitag.fasta");
		print RES $id."\n";
# 		print RES "$trouvfasta\n";
	}
	elsif ($numtrouve==1) # monotag donc OK
	{

		if($size >= $min && $size <= $max)
		{
		open(RES,">>".$path."/quality_output/seqAll_".$trouvechant.".fasta" );
		print RES ">".$trouvechant."_".$id."\n";
# 		print RES "$trouvfasta\n";
		  # On enregistre uniquement la séquances sans les primser 3/05/2016
		    if ( $kept_primers eq "no") { # Bug 17/05/2016
		    
		      if ($trouvfasta=~ /$forw_regexp(.*)$comp_rev_regexp/ ){ # Bug 17/05/2016
		      print RES "$1\n";
		      }
		      else {
		      print "Problème : ".$trouvechant."_".$id."\n";
		      }
		    }
		    else{
		    print RES "$trouvfasta\n";
		    }
		     
		# 3/05/2016
		}
		else
		{
		open (RES, ">>".$path."/quality_output/BAD/seqAll_bad_length.fasta");
                print RES $id."\n";
#                 print RES "$trouvfasta\n";
		}
	}
	else
	{
		open(RES,">>".$path."/quality_output/BAD/seqAll_withoutTags.fasta" );
		print RES $id."\n";	
# 		print RES $sequence{$id};	
	}
	
	close RES;    





  }
  elsif ($sens==1 && $fichiertag==0) # Absence du fichier tag 15/02/2016
  {
  $sequence{$id} =~ /.*($amplicon).*/i;
  open(RES,">>".$path."/quality_output/seqAll_sample.fasta" );
  print RES ">sample_".$id."\n";
  print RES "$1\n";
  }
			
  
}


close(BAD);


# Traitement des chimères


if ($CheckChimeras eq "yes") 
{
print LOG "Chimeras checking in progess\n";print "Chimeras checking in progess\n";
system("mkdir ".$path."/quality_output/chimeras");
system("mkdir ".$path."/quality_output/demultiplexing");
system("mv ".$path."/quality_output/seqAll_* ".$path."/quality_output/demultiplexing/"); # les séquencees démultipléxées totales sont stockées dans un dossier
$quality_output=$path."/quality_output/"; 
    
    if($fichiertag==1)
    {
        for ($i=0; $i<=$#tag; $i++)
        {
	chomp($tag[$i]);
	@ligne_tag=split("\t", $tag[$i]);
	$echantillon=$ligne_tag[0];

# Test si l'échantillon existe 12/05/2016
	if(-e $quality_output."demultiplexing/seqAll_".$echantillon.".fasta"){
	$command=$uclust." -derep_fulllength ".$quality_output."demultiplexing/seqAll_".$echantillon.".fasta  -output ".$quality_output."derep_".$echantillon.".fasta -sizeout >> ".$path."/panam.log";
	system($command);
	print LOG "Dereplication sample : $echantillon\n";

	check_chimeras_64 ($echantillon, $path, $uclust);
	}
	else {
	print LOG "Sample $echantillon empty\n";
	}
# 12/05/2016
        }
    }
    else
    {
        if(-e $quality_output."demultiplexing/seqAll_sample.fasta"){
	$command=$uclust." -derep_fulllength ".$quality_output."demultiplexing/seqAll_sample.fasta  -output ".$quality_output."derep_sample.fasta -sizeout >> ".$path."/panam.log";
	system($command);
	print LOG "Dereplication sample\n";

	check_chimeras_64 ("sample", $path, $uclust);
	}
	else {
	print LOG "Sample empty\n";
	}
    }

  		
}


close LOG;


}
######################################################################## Fin demultiplexage



######################################################################## check chimeras 
# check chimeras  # 5/06/2015
# A faire enlever le size des séquences chimera_free_
# 27/07/2015 : mapping
# 12/04/2016 : changement pour vsearch
# 12/05/2016 interface utilisateur : error.log
#############################



sub check_chimeras_64 {

	my $sample_name;
	my $fasta_dir;
	my $path;
	my $uclust;
	my $errorlog;
	my $command;my $nom;
	my %fasta_total; # Fasta

	($sample_name, $path, $uclust) = @_ ; 

$fasta_dir=$path."/quality_output/";
$errorlog=$path."/.error.log";

open (LOG, ">>".$path."/panam.log");
print LOG "Chimeras $sample_name in process...\n" ;

$command=$uclust." -uchime_denovo ".$fasta_dir."derep_".$sample_name.".fasta  -uchimeout ".$fasta_dir."chimeras/".$sample_name.".uchime -nonchimeras ".$fasta_dir."chimeras/nonchimeras_".$sample_name.".fasta 2>> $errorlog";
print LOG $command."\n";
system($command);

#mapping exact contre les séquences non chimériques
$command=$uclust." -search_exact $fasta_dir/demultiplexing/seqAll_".$sample_name.".fasta -db $fasta_dir/chimeras/nonchimeras_".$sample_name.".fasta  -uc $fasta_dir/chimeras/".$sample_name."_mapping.uc 2>> $errorlog";
print LOG $command."\n";
system($command);

open(CHIM, "$fasta_dir/chimeras/".$sample_name."_mapping.uc") || die ("UC not found");
open(FASTA,  "$fasta_dir/demultiplexing/seqAll_".$sample_name.".fasta") || die ("Fasta not found");
open(GOOD,   ">".$fasta_dir."seqAll_".$sample_name.".fasta");

while (my $l = <FASTA>){
	chomp ($l);
	if ($l =~ />(.*)/)
	{
	
	$nom=$1;
	}
	else	
	{
	$fasta_total{$nom}.=$l;

	}
}


while (my $l = <CHIM>){
	chomp ($l);
	
	my @tab = split ("\t", $l);

	print GOOD ">".$tab[8]."\n".$fasta_total{$tab[8]}."\n";
	
}

system("rm ".$fasta_dir."derep_".$sample_name.".fasta");

}

















#############################
# USEARCH  # le 20/05 DD
#10/01/2016 Ajout fichier log - Quality : -fastq_maxee 1
#9/05/2016 ajout de fasta not merged
#11/05/2016 redirection vers .error.log
#############################


sub merge_pairs_vsearch {
	## inputs
	my $MinSeqLength;
	my $MaxSeqLength;
	my $Nbase; 
	my $MinOverlap;
	my $MismatchOverlap;
	my $inputSeqNameF; 	 	
	my $inputSeqNameR;
	my $NGS_id_Results;
	my $merge;
	my $score;
	my $path_panam;

	($inputSeqNameF, $inputSeqNameR, $MinSeqLength, $MaxSeqLength, $Nbase, $MinOverlap, $MismatchOverlap, $NGS_id_Results, $merge, $score, $path_panam) = @_ ;

my $uclust = $path_panam."/vsearch" ;
my $command;
my $errorlog; #Redirection des erreurs


my $log_dir = $NGS_id_Results."/quality_output/pairs/";
$errorlog= $NGS_id_Results."/.error.log";


my $log_file;
$log_file=">>".$NGS_id_Results."/panam.log";
my $date=`date`;
open (LOG, $log_file);
print LOG "\n\t\t- Merging pairs with $uclust\n";
print LOG "Analysis starting: $date\n";

print  "\n\t\t- Merging pairs with $uclust\n";
print  "\t\t- Analysis starting: $date\n";


system("mkdir $log_dir");


# my $MaxDiffs= sprintf("%.0f",$MismatchOverlap*$MinOverlap); # A vérifier....
#my $MaxDiffs=5; #default values version 8
my $MaxDiffs=$MismatchOverlap;

# $command = $uclust." -fastq_mergepairs ".$inputSeqNameF." -reverse ".$inputSeqNameR." -fastq_minovlen ". $MinOverlap. " -fastq_minmergelen ".$MinSeqLength." -fastqout ".$log_dir."assembled.fastq"." -fastq_allowmergestagger -fastq_maxdiffs ".$MaxDiffs." -fastaout_notmerged_fwd ".$log_dir."not_merged_fwd.fasta -fastaout_notmerged_rev ".$log_dir."not_merged_rev.fasta 2>> $errorlog";

$command = $uclust." -fastq_mergepairs ".$inputSeqNameF." -reverse ".$inputSeqNameR." -fastq_minovlen ". $MinOverlap. " -fastq_minmergelen ".$MinSeqLength." -fastq_maxee 1  -fastaout ".$log_dir."merged.fasta"." -fastq_allowmergestagger -fastq_maxdiffs ".$MaxDiffs." -fastaout_notmerged_fwd ".$NGS_id_Results."/quality_output/BAD/not_merged_fwd.tmp -fastaout_notmerged_rev ".$NGS_id_Results."/quality_output/BAD/not_merged_rev.tmp 2>> $errorlog";

print LOG "$command\n";
system($command) ;
# $command = $uclust." -fastq_filter ".$log_dir."assembled.fastq"." -fastq_maxee 1  -fastaout ".$log_dir."merged.fasta  2>> $errorlog";
# print LOG "$command\n";
# system($command) ;

$command="rm ".$log_dir."tmp.fastq";
# print LOG "$command\n";
# system($command) ;

$date=`date`;
print LOG "Analysis ending $date\n";

$command="grep '>' ".$NGS_id_Results."/quality_output/BAD/not_merged_fwd.tmp > ".$NGS_id_Results."/quality_output/BAD/not_merged_fwd.fasta";
qx($command);
$command="rm ".$NGS_id_Results."/quality_output/BAD/not_merged_fwd.tmp ";
qx($command);

$command="grep '>' ".$NGS_id_Results."/quality_output/BAD/not_merged_rev.tmp > ".$NGS_id_Results."/quality_output/BAD/not_merged_rev.fasta";
qx($command);
$command="rm ".$NGS_id_Results."/quality_output/BAD/not_merged_rev.tmp ";
qx($command);




	return ($log_dir);
close LOG;

}







#############################
# PEAR  # le 21/05 DD
#
#############################


sub merge_pairs_pear {
	## inputs
	my $MinSeqLength;
	my $MaxSeqLength;
	my $Nbase; 
	my $MinOverlap;
	my $MismatchOverlap;
	my $inputSeqNameF; 	 	
	my $inputSeqNameR;
	my $NGS_id_Results;
	my $merge;
	my $score;
	my $path_panam;
	my $soft;

	($inputSeqNameF, $inputSeqNameR, $MinSeqLength, $MaxSeqLength, $Nbase, $MinOverlap, $MismatchOverlap, $NGS_id_Results, $merge, $score, $path_panam) = @_ ;



	#my $fasta_dir =  $NGS_id_Results."/fasta_dir/";
	my $log_dir = $NGS_id_Results."/quality_output/pairs/";



	if (-d $log_dir) {
		`rm  -r $log_dir`;
		`mkdir $log_dir`
	}
	else {
		`mkdir $log_dir`
	}
my $log_file;
$log_file=">>".$NGS_id_Results."/panam.log";

open (LOG, $log_file);
print LOG "Module merge_pairs_pear 10/01/2016 \n";

	
	my $command;

$command = "pear -f ".$inputSeqNameF." -r ".$inputSeqNameR." --min-overlap ". $MinOverlap. " --threads 16 -o ".$log_dir."tmp";

print LOG "Merged by PEAR\n";
print "$command\n\n";

system($command) ;
# $command = "usearch -fastq_filter ".$fasta_dir."tmp.assembled.fastq"." -fastaout ".$fasta_dir."seqAll_checkPrimers.fasta"; 
# system($command) ;

my $tag=0;
my %corresp ;

open(FQR,$log_dir."tmp.assembled.fastq") || die $log_dir."tmp.assembled.fastq not found \n";
open(FAR,">".$log_dir."merged.fasta") || die $log_dir."merged.fasta not found\n";
while(<FQR>){
	chomp($_);
	if (($_=~/^@/) && ($tag==0)){
		my $or = $_;
		my $id=$_ ; 
# 		$id =~ s/://g ;
# 		$id =~ s/-//g ;
# 		$id =~ s/\s//g ;
		$id =~ s/\@//g ; #print "$_ ---- $id\n";
		print FAR ">$id\n";
		$corresp{$or} = $id;
		

	}
	elsif($tag==1){
		print FAR "$_\n";
# 		print "$_\n";
		$tag=-3;
	}
	$tag++;
}
close FQR;
close FAR;



$command="rm ".$log_dir."tmp.assembled.fastq";
system($command) ;


#	`rm -r $fastq_split_dir`;
	return ($log_dir);

}






#############################
		# pandaseq  #
#############################


sub merge_pairs_pandaseq {
	## inputs
	my $MinSeqLength;
	my $MaxSeqLength;
	my $Nbase; 
	my $MinOverlap;
	my $MismatchOverlap;
	my $inputSeqNameF; 	 	
	my $inputSeqNameR;
	my $NGS_id_Results;
	my $merge;
	my $score;
	my $path_panam;

	($inputSeqNameF, $inputSeqNameR, $MinSeqLength, $MaxSeqLength, $Nbase, $MinOverlap, $MismatchOverlap, $NGS_id_Results, $merge, $score,$path_panam ) = @_ ;


	my $log_dir = $NGS_id_Results."/quality_output/pairs/";



	if (-d $log_dir) {
		`rm  -r $log_dir`;
		`mkdir $log_dir`
	}
	else {
		`mkdir $log_dir`
	}
my $log_file;
$log_file=">>".$NGS_id_Results."/panam.log";

open (LOG, $log_file);
print LOG "Merge pairs with pandaseq 10/01/2016 \n";
	## pandaseq


	my $panda_dir = "pandaseq" ;
	my $command;

	if ($merge eq "perfect_match") {
		#foreach my $k (keys %h) { 
			if ($Nbase eq "yes" ) {
				#$command = $panda_dir."./pandaseq -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -N -o ".$MinOverlap." -C completely_miss_the_point:".$MismatchOverlap." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";
				$command = $panda_dir." -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -N -o ".$MinOverlap." -C completely_miss_the_point:".$MismatchOverlap." -g ".$log_dir."seqAll_log.txt -w ".$log_dir."merged.fasta";
			 
			}

			if ($Nbase eq "no" ) {
				#$command = $panda_dir."./pandaseq -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -o ".$MinOverlap." -C completely_miss_the_point:".$MismatchOverlap." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";
				$command = $panda_dir." -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -o ".$MinOverlap." -C completely_miss_the_point:".$MismatchOverlap." -g ".$log_dir."seqAll_log.txt -w ".$log_dir."merged.fasta";

			}

			system($command) ;
		#}
	}

	elsif ($merge eq "fix_sequences") {
			#foreach my $k (keys %h) { 
				if ($Nbase eq "yes" ) {
					#$command = $panda_dir."./pandaseq -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -N -o ".$MinOverlap." -t ".$score." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";
					$command = $panda_dir." -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -N -o ".$MinOverlap." -t ".$score." -g ".$log_dir."seqAll_log.txt -w ".$log_dir."merged.fasta";

				}

				if ($Nbase eq "no" ) {
					#$command = $panda_dir."./pandaseq -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -o ".$MinOverlap." -t ".$score." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";
					$command = $panda_dir." -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -o ".$MinOverlap." -t ".$score." -g ".$log_dir."seqAll_log.txt -w ".$log_dir."merged.fasta";

				}

				#print "$command\n";
				system($command) ;
			#}
		}



print LOG "$command\n";

# 	my @fasta_file = <$fasta_dir*>;
# 	foreach my $fasta_file (@fasta_file) { 
# 		`sed -i s/:/-/g $fasta_file 2> /dev/null`;
# 		`sed -i s/ /-/g $fasta_file 2> /dev/null`;
# 	}
# 
# 	#`rm -r $fastq_split_dir`;
	return ( $log_dir);

}

1;
