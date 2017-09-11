#!/usr/bin/perl


use Cwd 'abs_path';
use FindBin;

use lib "$FindBin::RealBin/modules";


##################################################
my $path; my $path_panam; my $path_intallation; my $path_results; # path d'installation

my $id; # colonne similarité


##################################################

$path=abs_path($0);
# $path =~ /(.*)\/postprocess.*/; 
$path_panam="$FindBin::RealBin/";

$path_intallation=$path_panam."/R"; #


my ($USAGE) =  "\n\n\t****** USAGE of $0 PROGRAM ******\n\n\n\n\tUSAGE: perl $0 <panam.ini file> \n\n\n\n";
die "$USAGE" if ((scalar(@ARGV))< 1);
my $option_file = $ARGV[0];
chomp($option_file);

die "\n\n\t=> Cannot find configuration file: $option_file.\n\n" unless (-e $option_file);
die "\n\n\t=> Configuration file, $option_file, appears to be empty!\n\n" if (-z $option_file);

use parse_ini ;
my @parse = &parse_ini($option_file);

$path_results = $parse[0];
my $errorlog= $path_results."/.error.log";
die "\n\n\tOutput directory for NGS analyses is not defined. Check $option_file.\n\n" unless ( defined $path_results);



#####################Fin des modifs du 26 mai 2015

my $log_file;
$log_file=">>".$path_results."/panam.log";# Initialisation du fichier/début d'analyse

open (LOG, $log_file);
print LOG "\n\n* LAST STEP Processing final results and displaying it with HTML page PANAM2 \n\n";
print "\n\n* LAST STEP Processing final results and displaying it with HTML page PANAM2\n\n";





### Krona process #####

my $ktImportText= "$path_panam/bin/KronaTools-2.7/scripts/ImportText.pl";

`mkdir $path_results/samples 2>/dev/null`;

# From OTU_distribution_tax.txt 
open (TAG,"$path_results/OTU_distribution_tax.txt") or die "Couldn't open file OTU_distribution_tax.txt in the folder $path_results !\n";
@tag=<TAG>;
@ligne_label=split('\t',$tag[0]);


$nb_labels=$#ligne_label-5; # 2 premières colonnes + 4 affiliations taxonomiques (NN, LCA, identié, %)
$nn= $nb_labels+3;
$lca=$nn+1;
$id=$nn+2;

print "\nSamples processed: \n"; 

$ct=2;

while($ligne_label[$ct] !~ /NN.*/) # NN, premier label de la taxonomie
{
print "$ligne_label[$ct]\n";
$fichier="$ligne_label[$ct]_NN.txt";
$fichier_csv="$ligne_label[$ct]_NN.csv";
$colonne=$ct+1;

`cut  "\t" -f$colonne,$nn  $path_results/OTU_distribution_tax.txt >  $path_results/samples/$fichier 2>/dev/null`;
`sed "s/;/\t/g" $path_results/samples/$fichier | sed "1d" > $path_results/samples/$fichier_csv`;

$fichier="$ligne_label[$ct]_LCA.txt";
$fichier_csv="$ligne_label[$ct]_LCA.csv";
`cut  "\t" -f$colonne,$lca  $path_results/OTU_distribution_tax.txt > $path_results/samples/$fichier 2>/dev/null`;
`sed "s/;/\t/g" $path_results/samples/$fichier  | sed "1d" > $path_results/samples/$fichier_csv`;

$fichier="$ligne_label[$ct]_SIM.txt";
$fichier_csv="$ligne_label[$ct]_SIM.csv";
`cut  "\t" -f$colonne,$id  $path_results/OTU_distribution_tax.txt > $path_results/samples/$fichier 2>/dev/null`;
`sed "s/;/\t/g" $path_results/samples/$fichier  | sed "1d" > $path_results/samples/$fichier_csv`;



$ct++;

}

print  "\nTaxonomic distribution processed by Krona...\n";

`mkdir $path_results/Figures 2>/dev/null`;
`$ktImportText $path_results/samples/*_NN.csv -o $path_results/Figures/Fig_taxo_distribution_NN.html `;
`$ktImportText $path_results/samples/*_LCA.csv -o $path_results/Figures/Fig_taxo_distribution_LCA.html`;
`$ktImportText $path_results/samples/*_SIM.csv -o $path_results/Figures/Fig_taxo_distribution_SIM.html`;


# From OTU_distribution_tax_normalized_XX.txt (n° of samples could be different from non normalized data)

open (NN,"$path_results/OTU_distribution_tax_normalized_NN.txt");
@NN=<NN>;
@ligne_tag=split('\t',$NN[0]);


$nb_labels=$#ligne_tag-1;
$taxonomy= $#ligne_tag+1; # indice du cut

$ct=1;

	while($ct <= $nb_labels)
	{

	$fichier="$ligne_tag[$ct]_normalized_NN.txt";
	$fichier_csv="$ligne_tag[$ct]_normalized_NN.csv";
	$colonne=$ct+1;

	

	`cut  "\t" -f$colonne,$taxonomy  $path_results/OTU_distribution_tax_normalized_NN.txt >  $path_results/samples/$fichier 2>/dev/null`;
	`sed "s/;/\t/g" $path_results/samples/$fichier | sed "1d" > $path_results/samples/$fichier_csv`;


	$fichier="$ligne_tag[$ct]_normalized_LCA.txt";
	$fichier_csv="$ligne_tag[$ct]_normalized_LCA.csv";
	`cut  "\t" -f$colonne,$taxonomy  $path_results/OTU_distribution_tax_normalized_LCA.txt > $path_results/samples/$fichier 2>/dev/null`;
	`sed "s/;/\t/g" $path_results/samples/$fichier  | sed "1d" > $path_results/samples/$fichier_csv`;


	$ct++;

	}

if (`ls $path_results/samples/*_normalized_*.csv | wc -l 2>/dev/null`> 0)
{
`$ktImportText $path_results/samples/*_normalized_NN.csv -o $path_results/Figures/Fig_taxo_distribution_normalized_NN.html `;
`$ktImportText $path_results/samples/*_normalized_LCA.csv -o $path_results/Figures/Fig_taxo_distribution_normalized_LCA.html`;
}

`rm -R $path_results/samples 2>/dev/null`;

print "Taxonomic distribution can be viewed by any browser in the directory: $path_results/Figures\n";



### End Krona process #####



### Begin  R process #####

#print "Richness, diversity, rarefaction curves and heatmap processed by R with vegan and phyloseq packages...\n";

if(`which R | wc -l` > 0) 
{

`mkdir $path_results/R_output 2>/dev/null`;
open (OTU,"$path_results/OTU_distribution_tax.txt");


# Le 5/7/2017
open (R, ">".$path_results."/R_output/OTU_distribution_phyloseq.txt");

my $li;
my @table=<OTU>;
my $taxo;
my $i; my $j;
my @count;
my $taxo_max=0;
my @lign;
my $identite;
my $best;
my $lca;
my $nn;
my @ligne;
my $taxolapluslongue;




print R $table[0];

for ($j=1; $j<=$#table;$j++) # Recherche de la profondeur taxo max
{
chomp $table[$j];
 @ligne=split("\t", $table[$j]);
 $identite=pop(@ligne);
 $best=pop(@ligne);
 $lca=pop(@ligne);
 $nn=pop(@ligne);

foreach $taxo ($nn, $lca, $best)
    {
    
        if ($taxo =~ /(.*);\s;$/)
        {
        $taxo=$1.";";
        }
        if ($ taxo !~ /.*;$/)
        {
        $taxo=$taxo.";";
        }
    @count = ($taxo =~ /;/g);
    
    if ( @count > $taxo_max){ $taxo_max= @count ;$taxolapluslongue=$taxo;}
   
    }



}

# print "Profondeur taxonomique max : $taxo_max ->  $taxolapluslongue\n";

for ($j=1; $j<=$#table;$j++)
{
chomp $table[$j];
 @ligne=split("\t", $table[$j]);
 $identite=pop(@ligne);
 $best=pop(@ligne);
 $lca=pop(@ligne);
 $nn=pop(@ligne);

print R join("\t",@ligne);
print R "\t";

    foreach $taxo ($nn, $lca, $best)
    { 
        if ($ taxo =~ /(.*);\s;$/)
        {
        $taxo=$1.";";
        }
        if ($ taxo !~ /.*;$/)
        {
        $taxo=$taxo.";";
        }
    @count = ($taxo =~ /;/g);
        for ($i=@count; $i<$taxo_max;$i++)
        {
        $taxo=$taxo."Unclassified;";
        }
    print R $taxo."\t";
    }
    
    print R $identite."\n";
    
}





# R -> files : Richness_diversity.txt
my $command="R --vanilla --args ".$path_results."/R_output/ < ".$path_intallation."/post_process_v4.R 2>> ".$errorlog;# 
print LOG $command."\n";
qx($command);



	if (`ls $path_results/R_output/*.jpg | wc -l 2>/dev/null`> 0)
	{
	# Formating R results
	`sed "s/Observed/\tObserved/g" $path_results/.Richness_diversity.tmp > $path_results/Richness_diversity.txt`;
	
	
		if(-e "$path_results/.Richness_diversity_norm.tmp")
		{
# 		`printf "\n\nNormalized data\n" >> $path_results/Richness_diversity.txt`;
		`sed "s/Observed/\tObserved/g" $path_results/.Richness_diversity_norm.tmp > $path_results/Richness_diversity_normalized.txt`;
		}
	#print "Richness and diversity indices computed after phylogenetic affiliation are in the file: Richness_diversity.txt\n";

	`cp   $path_results/R_output/*.jpg   $path_results/Figures 2>/dev/null `;
	#print "Rarefaction curves and heatmap generated\n";
	print "The R results and the phyloseq objects were saved in the directory: $path_results/R_output\nOpen a terminal and type R for processing other analyses with phyloseq and vegan packages\n\n";
	}
	else
	{
	print "Process R with phyloseq and/or vegan aborted, check the R packages installed\n\n";
	
	}

 `rm $path_results/.*.tmp 2>/dev/null`;
 `rm $path_results/R_output/*.tmp 2>/dev/null`;

}
else
{
print "R software is not installed, process aborted ... \n"; 
}
close TAX; 
close CSV;


### End  R process #####



########### Partie Emilie - 07/06/2017 - Partie affichage de résultats sur page HTML (panam.html) ############
print "HTML page in process..\n";


if (-d "$path_results/.CSS"){
	`rm -r $path_results/.CSS`;
}

if (-d "$path_results/.res"){
	`rm -r $path_results/.res`;
}

if (-d "$path_results/.CSS/.taxonomy"){
	`rm -r $path_results/.CSS/.taxonomy`;
}

unless (-d "$path_results/PhyloDiv_output"){
	print "For displaying missing results (phylogenetic indices) in this HTML report: run phylodiv_panam.pl and postprocess_panam.pl\n\n";
}


system("cp -r $path_panam/.CSS $path_results/.CSS");

`cp $path_panam/.CSS/.panam.html $path_results/panam.html`;

my $command="bash $path_results/.CSS/split-table.sh $path_results";
system($command);


`mkdir $path_results/.CSS/.res`;



#~ scalar
my $forward;
my $reverse;
my $clustering;
my $min_seq;
my $max_seq;
my $overlap;
my $mismatch;
my $chimera;
my $abundance=0;
my $clustering_pourcent;
my $primer=0;
my $n=0;
my $length=0;
my $assemble1=0;
my $assemble2=0;
my $assemble;
my $good=0;
my $nb_bad;
my $nb_good=0;
my $nb;
my $bad;
my $bad_otu;
my $good_otu;
my $i;
my $j;
my $head;
my $ligne;
my $headn;
my $lignen;
my $sample_d;
my $head_clade;
my $ligne_c;
my $p="";
my $seuil;


# tab
my @domain; 
my @tab_bad;
my @tab_good;
my @head;
my @value;
my @sample_d;
my @headn;
my @valuen;
my @sample_dn;
my @head_clade;
my @value_clade;
my @t_c;
my @taxo_c;
my @id_c;
my @chemin_c;
my @clade_c;
my @lien_c;
my @tab_c;
my @c_c;
my @boot_c;
my @otu_c;
my @seq_c;
my @sing_c;
my @sample_c;
my @MMNDa_c;
my @MMNDp_c;
my @depth_c;


# hash
my %obs;
my %schao;
my %se_chao;
my %ace;
my %se_ace;
my %shannon;
my %simpson;
my %inv_simpson;
my %fisher;
my %obsn;
my %schaon;
my %se_chaon;
my %acen;
my %se_acen;
my %shannonn;
my %simpsonn;
my %inv_simpsonn;
my %fishern;
my %clade_c;
my %boot_c;
my %otu_c;
my %seq_c;
my %sing_c;
my %sample_c;
my %MMNDa_c;
my %MMNDp_c;
my %depth_c;
my %NN_c;
my %taxo_c;
my %way_c;


##### Recuperation parametre #####

open(PARA, ">$path_results/.CSS/.res/parametre.html");

my @domain;

foreach my $d (keys %{$parse[32]}){
	$dom{$d}=1;
	push @domain, $d ;
}
$forward=$parse[15];
$reverse=$parse[16];
$clustering=$parse[20];
$min_seq=$parse[4];
$max_seq=$parse[5];
$overlap=$parse[6];
$mismatch=$parse[7];
$chimera=$parse[27];
$abundance=$parse[29];

print PARA "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />
	<center><table width=\"90%\">
	<tr>
		<td>Domain</td>
		<td>";

	for(my $i=0; $i<scalar(@domain); $i++){
		print PARA " $domain[$i] /";
	}
print PARA "</td>
	</tr>
	<tr bgcolor=\"#a9ccfb\">
		<td width=\"50%\">Reverse primer</td>
		<td> $reverse</td>
	</tr>
	<tr bgcolor=\"#a9ccfb\">
		<td>Forward primer</td>
		<td> $forward</td>
	</tr>
	<tr>
		<td>Clustering cutoff</td>
		<td> $clustering</td>
	</tr>
	<tr bgcolor=\"#a9ccfb\">
		<td width=\"auto\">Lower sequence length cutoff </td>
		<td width=\"250px\"> $min_seq</td>
	</tr>
	<tr bgcolor=\"#a9ccfb\">
		<td>Max sequence length cutoff</td>
		<td> $max_seq</td>
	</tr>
	<tr>
		<td>Min overlap length</td>
		<td> $overlap</td>
	</tr>
	<tr>
		<td>Mismatch overlap</td>
		<td> $mismatch</td>
	</tr>
	<tr bgcolor=\"#a9ccfb\">
		<td>Check chimeras</td>
		<td> $chimera</td>
	</tr>
	<tr>
		<td> Abundance cut-off</td>
		<td> $abundance</td>
	</tr>";



$clustering_pourcent=$clustering*100;



##### Graphic quality #####

open(CAMEMBERT, ">$path_results/.CSS/.res/camembert.json");

if($parse[34] eq ""){

	open(PRIMER, $path_results."/quality_output/BAD/seqAll_bad_primers.fasta");
	while (my $lic=<PRIMER>){
		$primer++;
	}
	close PRIMER;

	open(N, $path_results."/quality_output/BAD/seqAll_bad_ATGC.fasta");
	while (my $lic=<N>){
		$n++;
	}
	close N;

	open(LENGTH, $path_results."/quality_output/BAD/seqAll_bad_length.fasta");
	while (my $lic=<LENGTH>){
		$length++;
	}
	close LENGTH;

	open (AS, $path_results."/quality_output/BAD/not_merged_fwd.fasta");
	while (my $lic=<AS>){
		if($lic=~/^>/){
			$assemble1++;
		}
	}
	close AS;

	open (AS, $path_results."/quality_output/BAD/not_merged_rev.fasta");
	while (my $lic=<AS>){
		if($lic=~/^>/){
			$assemble2++;
		}
	}
	close AS;

	if($assemble1<$assemble2){
		$assemble=$assemble2;
	}
	else{
		$assemble=$assemble1;
	}



	open(GOOD, $path_results."/quality_output/seqAll.fasta");
	while (my $lic=<GOOD>){
		if($lic=~/^>/){
			$nb++;
		}
		$good=$nb;
	}
	close GOOD;
}
else{
	
	my @f=`ls $path_results/quality_output/*/quality_output/BAD/*`;

	foreach my $f (@f){

		if($f =~ m/seqAll_bad_primers.fasta/){
			open (FILE, $f);
			while (my $lic=<FILE>){
				$primer++;
			}
		
		}
	
	
		if($f =~ m/seqAll_bad_ATGC.fasta/){
			open (FILE, $f);
			while (my $lic=<FILE>){
				$n++;
			}
		}
	
		if($f =~ m/seqAll_bad_length.fasta/){
			open (FILE, $f);
			while (my $lic=<FILE>){
				$length++;
			}
		}
	
		if ($f =~ m/not_merged_fwd.fasta/){
			open (FILE, $f);
			while (my $lic=<FILE>){
				if($lic=~/^>/){
					$assemble1++;
				}
			}
		}


		if ($f =~ m/not_merged_rev.fasta/){
			open (FILE, $f);
			while (my $lic=<FILE>){
				if($lic=~/^>/){
					$assemble2++;
				}
			}
		}

		$assemble=$assemble2+$assemble1;

	
	}

	open (FILE, "$path_results/quality_output/seqAll.fasta");
	while (my $lic=<FILE>){
		if($lic=~/^>/){
			$nb++;
		}
		$good=$nb;
	}
}

print CAMEMBERT "[{\"k\": \"without primers\" , \"v\":$primer},{\"k\": \"with N\" , \"v\":$n},{\"k\": \"bad length\" , \"v\":$length},{\"k\": \"non-merged\" , \"v\":$assemble},{\"k\": \"sequences for clustering\" , \"v\":$good}]";
close CAMEMBERT;





 #### Tableau quality ####
open(QUAL, ">$path_results/.CSS/.res/seq_nb_quality.html");
print QUAL "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />";
if($parse[34] eq ""){
	my %corres;

	my $sample = `grep -c ">" $path_results/quality_output/seqAll_* | cut -d. -f 1 | cut -d_ -f 3`;

	my $nb = `grep -c ">" $path_results/quality_output/seqAll_* | cut -d: -f 2`;

	my @nb=split("\n", $nb);
	my @sample=split("\n", $sample);



	for (my $i=0; $i<scalar(@nb); $i++) {
		$corres{$sample[$i]}=$nb[$i];
	}

	print QUAL "<table> <caption>Sequences number by sample </caption>
				<tr> 
					<th>Sample</th>
					<th>Number</th>
				</tr>";

	foreach my $key (keys %corres){
		print QUAL "
				<tr>
					<td width=\"auto\" height=\"auto\">$key</td>
					<td width=\"75\"><center>$corres{$key}</center></td>
				</tr>";
	}

	print QUAL "</table>";
	
}
else{

	my @f=`ls $path_results/quality_output/*/quality_output/seqAll_sample.fasta`;
	
	print QUAL "<table> <caption>Sequences number by sample </caption>
				<tr> 
					<th>Sample</th>
					<th>Number</th>
				</tr>
			";
	
	foreach my $f (@f){

	my $sample="";
	my $nb=0;

		open(FILE, $f) or die ("Erreur d'ouverture");


		while (my $lic=<FILE>){
			if ($lic=~/^>(\S+)_/){
				$sample=$1;
				$nb++;
			}

		}
		print QUAL "
				<tr><td>$sample</td>
				<td>$nb</td></tr>";
		close FILE;
	}
	print QUAL "
			</table>";
	
}

close QUAL;


##### Calcul number OTU & sequence + pourcent of clustering & abundance#####

open(CLUSTERING, ">$path_results/.CSS/.res/clustering.json");
open(NB_SEQ, ">$path_results/.CSS/.res/nb_seq.html");

open(C_OTU, ">$path_results/.CSS/.res/otu.json");
open(NB_OTU, ">$path_results/.CSS/.res/nb_otu.html");

open (BAD, "$path_results/preprocess_output/pooled_sample/BAD_OTU");
while(my $lic = <BAD>){
	if($lic!~/^OTU/){
		if($lic=~/^\S+\s+(\d+)/){
			$bad=$1;
		}
		$i++;
	}
	push @tab_bad, $bad;
	
	$bad_otu=$i;
	
}
close BAD;

for(my $i=1; $i<scalar(@tab_bad); $i++){
	$nb_bad+=$tab_bad[$i];
}
	

open (GOOD, "$path_results/preprocess_output/pooled_sample/pooled_sample_OTU");
while(my $lic = <GOOD>){
	if($lic=~/OTU_\d+\t\S+\s(\d+)/){
		$nb=$1;
		$j++;
	}
	push @tab_good, $nb;
		
	$good_otu=$j;
	
}
close GOOD;
	
	for(my $i=1; $i<scalar(@tab_good); $i++){
		$nb_good+=$tab_good[$i];
	}
	

print CLUSTERING "[{\"k\": \"sequences removed \" , \"v\" : $nb_bad},{\"k\": \"Sequences conserved\" , \"v\" : $nb_good}]";
close CLUSTERING;

my $nb_seq=$nb_good+$nb_bad;
print NB_SEQ "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />Sequences number : $nb_seq\n";
close NB_SEQ;


print C_OTU "[{\"k\": \"OTUs removed \" , \"v\" : $bad_otu},{\"k\": \"OTUs conserved\" , \"v\" : $good_otu}]";
close C_OTU;


my $nb_otu=$bad_otu+$good_otu;
print NB_OTU "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />Clusters number : $nb_otu\n";
close NB_OTU;

open (AB, ">$path_results/.CSS/.res/clust.html");
print AB "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" /><table><td class=\"ab\" width=\"40%\"><center><br> Abundance threshold : <b> <$parse[29] </b><br>";
print AB " clustering cutoff : <b> $clustering_pourcent% </b><br><br></center></td><table>";
close AB;





#### Import du tableau Richness_diversity ###

open (DIV, $path_results."/Richness_diversity.txt");

while (my $lic = <DIV>){
		if($lic=~/^\t\w+/){
		$head=$lic;
		@head=split("\t", $head);
	}

	if($lic=~/^\w+/){
	
		$ligne=$lic;
	
		@value=split("\t", $ligne);
	
		push(@sample_d, $value[0]);
		$sample_d=$value[0];
		
		$obs{$sample_d}=$value[1];
		$chao{$sample_d}=$value[2];
		$se_chao{$sample_d}=$value[3];
		$ace{$sample_d}=$value[4];
		$se_ace{$sample_d}=$value[5];
		$shannon{$sample_d}=$value[6];
		$simpson{$sample_d}=$value[7];
		$inv_simpson{$sample_d}=$value[8];
		$fisher{$sample_d}=$value[9];
	}

}
close DIV;

open (DIVERSITY, ">$path_results/.CSS/.res/richness_div.html");

print DIVERSITY "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />

	<table >";
	
foreach my $ent (@head){
	print DIVERSITY "<th><center>$ent</center></th>";
}

print DIVERSITY "</tr>";
	
foreach my $s (@sample_d){
	print DIVERSITY "<tr>
		<td width=\"200\">$s</td>
		<td width=\"100\"><center>$obs{$s}</center></td>
		<td width=\"100\"><center>$chao{$s}</center></td>
		<td width=\"100\"><center>$se_chao{$s}</center></td>			
		<td width=\"100\"><center>$ace{$s}</center></td>
		<td width=\"100\"><center>$se_ace{$s}</center></td>
		<td width=\"100\"><center>$shannon{$s}</center></td>
		<td width=\"100\"><center>$simpson{$s}</center></td>
		<td width=\"100\"><center>$inv_simpson{$s}</center></td>
		<td width=\"100\"><center>$fisher{$s}</center></td>
	</tr>";
}

print DIVERSITY "</table>";

close DIVERSITY;



### Import Richness diversity normalized ###
open (DIVN, $path_results."/Richness_diversity_normalized.txt");

while (my $lic = <DIVN>){
		if($lic=~/^\t\w+/){
		$headn=$lic;
		@headn=split("\t", $head);
	}

	if($lic=~/^\w+/){
	
		$lignen=$lic;
	
		@valuen=split("\t", $lignen);
	
		push(@sample_dn, $valuen[0]);
		$sample_dn=$valuen[0];
		
		$obsn{$sample_dn}=$valuen[1];
		$chaon{$sample_dn}=$valuen[2];
		$se_chaon{$sample_dn}=$valuen[3];
		$acen{$sample_dn}=$valuen[4];
		$se_acen{$sample_dn}=$valuen[5];
		$shannonn{$sample_dn}=$valuen[6];
		$simpsonn{$sample_dn}=$valuen[7];
		$inv_simpsonn{$sample_dn}=$valuen[8];
		$fishern{$sample_dn}=$valuen[9];
	}



}
close DIVN;

open (DIVERSITYN, ">$path_results/.CSS/.res/richness_div_normalized.html");

print DIVERSITYN "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />

	<table >";
	
foreach my $ent (@headn){
	print DIVERSITYN "<th><center>$ent</center></th>";
}

print DIVERSITYN "</tr>";
	
foreach my $s (@sample_dn){
	print DIVERSITYN "<tr>
		<td width=\"200\">$s</td>
		<td width=\"100\"><center>$obsn{$s}</center></td>
		<td width=\"100\"><center>$chaon{$s}</center></td>
		<td width=\"100\"><center>$se_chaon{$s}</center></td>			
		<td width=\"100\"><center>$acen{$s}</center></td>
		<td width=\"100\"><center>$se_acen{$s}</center></td>
		<td width=\"100\"><center>$shannonn{$s}</center></td>
		<td width=\"100\"><center>$simpsonn{$s}</center></td>
		<td width=\"100\"><center>$inv_simpsonn{$s}</center></td>
		<td width=\"100\"><center>$fishern{$s}</center></td>
	</tr>";
}

print DIVERSITYN "</table>";

close DIVERSITYN;



#### Arbre par groupe phylo ####
open (GP, ">$path_results/.CSS/.res/gp_arbre.html");
open (GP1, ">$path_results/.CSS/.res/gp_arbre1.html");

if (-d "$path_results/PhyloDiv_output"){
	@f=`ls $path_results/PhyloDiv_output`;

	my $gp;
	my @gp;
	my $img;
	my @img;
	my @image;



	foreach my $f (@f){

		if (($f =~/^(\w+|\w+\d+)_cp\S+.newick/g) || ($f =~/^(\w+|\w+\d+)_\S+.newick/g) || ($f =~/^(\w+)_\S+.\d+_\S+.newick/g)){
			$gp=$1;
			push @gp, $gp;
		}

		@img=`ls $path_results/PhyloDiv_output/$f`;

		foreach my $i (@img){
			if($i=~/(\S+_unifrac_cluster.svg)/){
				push @image, $i;
			}
		}
	}


	

	print GP "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />
			<ul>
				<div id=\"div3\"><li style=\"list-style-type:none\"><a href=\"PhyloDiv_output/$f[0]/$image[0]\" download><img src=\"http://www.kegg.jp/Fig/get_htext/close.png\"> $gp[0]</a></li>
			
			";

	for(my $j=1; $j<scalar(@gp); $j++){
		print GP "<li style=\"list-style-type:none\"><a href=\"PhyloDiv_output/$f[$j]/$image[$j]\" download><img src=\"http://www.kegg.jp/Fig/get_htext/close.png\"> $gp[$j]</a></li>";
	}

	print GP "</div><div id=\"div4\"><li style=\"list-style-type:none\"><center><img src=\"PhyloDiv_output/$f[0]/$image[0]\" height=\"500\"></center></li></div></ul><div class=\"clear\"></div><br><br><br><br><br>";
}
else{
	print GP1 "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />
			<br><center><b>No result</b></center>";
}
close GP;
close GP1;





#### Phylogenetic indices ######



if(-d "$path_results/PhyloDiv_output"){
	open (CLADE, "$path_results/Clade_Results.txt");

	while (my $lic=<CLADE>){
		if ($lic=~/^clade_/){
			$ligne_c=$lic;
			@value_clade=split("\t", $ligne_c);
		
			my $clade=$value_clade[0];
		
			$clade_c{$clade}=$value_clade[0];
			$boot_c{$clade}=$value_clade[1];
			$otu_c{$clade}=$value_clade[2];
			$seq_c{$clade}=$value_clade[3];
			$sing_c{$clade}=$value_clade[4];
			$sample_c{$clade}=$value_clade[5];
			$MNNDa_c{$clade}=$value_clade[6];
			$MNNDp_c{$clade}=$value_clade[7];
			$depth_c{$clade}=$value_clade[8];
			$tree_c{$clade}=$value_clade[9];
			$NN_c{$clade}=$value_clade[10];
			$taxo_c{$clade}=$value_clade[11];
		
			@taxo_c=$taxo_c{$clade};
			push @NN_c, $NN_c{$clade};
			push @chemin_c, $tree_c{$clade};
			push @clade_c, $clade_c{$clade};
			push @boot_c, $boot_c{$clade};
			push @otu_c, $otu_c{$clade};
			push @seq_c, $seq_c{$clade};
			push @sing_c, $sing_c{$clade};
			push @sample_c, $sample_c{$clade};
			push @MNNDa_c, $MNNDa_c{$clade};
			push @MNNDp_c, $MNNDp_c{$clade};
			push @depth_c, $depth_c{$clade};		

			for (my $i=0; $i<scalar(@taxo_c); $i++){
				$taxo_c[$i]=~ s/;/;\t/g;
				push @t_c, $taxo_c[$i];
			}		
		}
		else{
			$head_clade=$lic;
			@head_clade=split("\t", $head_clade);
		}


		for(my $i=0; $i<scalar(@chemin_c); $i++){
			$way_c{$NN_c[$i]}=$chemin_c[$i];
		}

		foreach my $k (keys %way_c){
			@lien_c = join("/", $way_c{$k}, $k);
		}

		for(my $i=0; $i<scalar(@lien_c); $i++){
			push @tab_c, $lien_c[$i];
		}
	
	}
	close CLADE;

	open (TC, ">$path_results/.CSS/.res/clade.html");

	print TC "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />
			<center><table id=\"\" width=\"1700\"><tr>";


	print TC "<th><center>$head_clade[0]</center></th>
			<th><center>$head_clade[10]</center></th>
			<th><center>$head_clade[1]</center></th>
			<th><center>$head_clade[2]</center></th>
			<th><center>$head_clade[3]</center></th>
			<th><center>$head_clade[4]</center></th>
			<th><center>$head_clade[5]</center></th>
			<th><center>$head_clade[6]</center></th>
			<th><center>$head_clade[7]</center></th>
			<th><center>$head_clade[8]</center></th>
			<th><center>$head_clade[11]</center></th>";
	
	print TC "</tr>";
	
	print TC "<tr>";

	for(my $i=0; $i<scalar(@chemin_c); $i++){
	print TC "		
			<td width=\"70\"><a href=\"../../PhyloDiv_output/$tab_c[$i]_clade.svg\" download style=\"color:blue\"><u> $clade_c[$i] </u></a></td>
			<td width=\"auto\"><center> $NN_c[$i] </center></td>
			<td width=\"auto\"><center> $boot_c[$i] </center></td>
			<td width=\"auto\"><center> $otu_c[$i] </center></td>			
			<td width=\"auto\"><center> $seq_c[$i] </center></td>
			<td width=\"auto\"><center> $sing_c[$i] </center></td>
			<td width=\"auto\"><center> $sample_c[$i] </center></td>
			<td width=\"auto\"><center> $MNNDa_c[$i] </center></td>
			<td width=\"auto\"><center> $MNNDp_c[$i] </center></td>
			<td width=\"auto\"><center> $depth_c[$i] </center></td>
			<td width=\"auto\"> $t_c[$i] </td>
		
		
		</tr>";
	}


	print TC "</table></center><br>";
}
else{
	open (TC1, ">$path_results/.CSS/.res/clade1.html");
	print TC1 "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />
			<br><center><b>No result</b></center>";
}
close TC;
close TC1;




##### Normalized LCA ##### 

my @lca=`ls $path_results/.CSS/.taxonomy/taxonomic_distribution_LCA.*.csv`;

my @sample_lca;

foreach my $file (@lca){

#~ scalar
my $sample;
my $org;
my $ligne;
my $entete;

#~ tab
my @entete;
my @org;
my @value;

#~ hash
my %coverage;
my %seq;
my %otu;
my %schao;
my %shannon;

open(FILE, $file);

	while(my $lic=<FILE>){
	
		if($lic=~/^(\S+)\n/){
			$sample=$1;
		}
	
		if($lic=~/(#\w+\t#\w+\t\w+\t\w+\t\w+)/){
			$entete=$1;
			@entete=split("\t", $entete);
		}
	
		if(($lic=~/^\S+\t/) || ($lic=~/^\w+\s\w+\s\w+\t/)){
		
			$ligne=$lic;
		
			@value=split("\t", $ligne);
		
		
			if($value[1]!=0){
				push(@org, $value[0]);
				$org=$value[0];
				
				$seq{$org}=$value[1];
				$otu{$org}=$value[2];
				$schao{$org}=$value[3];
				$shannon{$org}=$value[4];
				$coverage{$org}=$value[5];
				

			}

		}
	
	}
	close FILE;



	open (OUT, ">>$path_results/.CSS/.res/taxo_LCA.html");

	print OUT "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />
		<hr class=\"ligne1\"><br>	
		<table id=\"$sample._lca\"><caption>$sample</caption>
			<tr><td> </td>";
		
	foreach my $ent (@entete){
		print OUT "<th><center>$ent</center></th>";
	}
	print OUT "</tr>";

	foreach my $org (@org){
		print OUT "<tr>
			<td>$org</td>
			<td width=\"100\"><center>$seq{$org}</center></td>
			<td width=\"100\"><center>$otu{$org}</center></td>
			<td width=\"100\"><center>$schao{$org}</center></td>
			<td width=\"100\"><center>$shannon{$org}</center></td>
			<td width=\"100\"><center>$coverage{$org}</center></td>
		</tr>";
	}

	print OUT "</table>";

	push @sample_lca, $sample;
}
close OUT;






##### No normalized LCA #####

my @lcan=`ls $path_results/.CSS/.taxonomy/taxonomic_distribution_LCA_norm.*.csv`;

my @sample_lcan;

foreach my $file (@lcan){

#~ scalar
my $sample;
my $org;
my $ligne;
my $entete;


#~ tab
my @entete;
my @org;
my @value;


#~ hash
my %coverage;
my %seq;
my %otu;
my %schao;
my %shannon;

open(FILE, $file);
$seuil=0;
	while(my $lic=<FILE>){
	
		if($lic=~/^(\S+)\n/){
			$sample=$1;
		}
	
		if($lic=~/(#\w+\t#\w+\t\w+\t\w+\t\w+)/){
			$entete=$1;
			@entete=split("\t", $entete);
		}
	
		if(($lic=~/^\S+\t/) || ($lic=~/^\w+\s\w+\s\w+\t/)){
		
			$ligne=$lic;
		
			@value=split("\t", $ligne);
		
		
			if(($value[1]ne"NA")&&($value[2]ne"NA")&&($value[3]ne"NA")&&($value[4]ne"NA")&&($value[5]ne"NA")
			&&($value[1]!=0)){
				push(@org, $value[0]);
				$seuil+=$value[1];
				my $org=$value[0];
				$seq{$org}=$value[1];
				$otu{$org}=$value[2];
				$schao{$org}=$value[3];
				$shannon{$org}=$value[4];
				$coverage{$org}=$value[5];
				

			}
		}
	}
	close FILE;


	if(($value[1]ne"NA")){
		open (LCAN, ">>$path_results/.CSS/.res/taxo_LCA_norm.html");

		print LCAN "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />
			<hr class=\"ligne1\"><br>	
			<table id=\"$sample._lcan\"><caption>$sample</caption>
				<tr><td> </td>";
		
		foreach my $ent (@entete){
			print LCAN "<th><center>$ent</center></th>";
		}	
		print LCAN "</tr>";

		foreach my $org (@org){
			push (@seq, $seq{$org});
			print LCAN "<tr>
				<td>$org</td>
				<td width=\"100\"><center>$seq{$org}</center></td>
				<td width=\"100\"><center>$otu{$org}</center></td>
				<td width=\"100\"><center>$schao{$org}</center></td>
				<td width=\"100\"><center>$shannon{$org}</center></td>
				<td width=\"100\"><center>$coverage{$org}</center></td>
			</tr>";
			
		}

		print LCAN "</table>";

		push @sample_lcan, $sample;
	}

}
close LCAN;





##### Normalized NN #####

my @nn=`ls $path_results/.CSS/.taxonomy/taxonomic_distribution_NN.*.csv`;

my @sample_nn;

foreach my $file (@nn){

#~ scalar
my $sample;
my $org;
my $ligne;
my $entete;

#~ tab
my @entete;
my @org;
my @value;

#~ hash
my %coverage;
my %seq;
my %otu;
my %schao;
my %shannon;

open(FILE, $file);

	while(my $lic=<FILE>){
	
		if($lic=~/^(\S+)\n/){
			$sample=$1;
		}
	
		if($lic=~/(#\w+\t#\w+\t\w+\t\w+\t\w+)/){
			$entete=$1;
			@entete=split("\t", $entete);
		}
	
		if(($lic=~/^\S+\t/) || ($lic=~/^\w+\s\w+\s\w+\t/)){
		
			$ligne=$lic;
		
			@value=split("\t", $ligne);
		
			if($value[1]!=0){
				push(@org, $value[0]);
				my $org=$value[0];
				$seq{$org}=$value[1];
				$otu{$org}=$value[2];
				$schao{$org}=$value[3];
				$shannon{$org}=$value[4];
				$coverage{$org}=$value[5];
				
			}
		}
	}
	close FILE;



	open (NN, ">>$path_results/.CSS/.res/taxo_NN.html");

	print NN "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />
		<hr class=\"ligne1\"><br>	
		<table id=\"$sample._nn\"><caption>$sample</caption>
			<tr><td> </td>";
		
	foreach my $ent (@entete){
		print NN "<th><center>$ent</center></th>";
	}	
	print NN "</tr>";

	foreach my $org (@org){
		print NN "<tr>
			<td>$org</td>
			<td width=\"100\"><center>$seq{$org}</center></td>
			<td width=\"100\"><center>$otu{$org}</center></td>
			<td width=\"100\"><center>$schao{$org}</center></td>
			<td width=\"100\"><center>$shannon{$org}</center></td>
			<td width=\"100\"><center>$coverage{$org}</center></td>
		</tr>";
	}

	print NN "</table>";

	push @sample_nn, $sample;
}
close NN;





##### No normalized NN #####

my @nnn=`ls $path_results/.CSS/.taxonomy/taxonomic_distribution_NN_norm.*.csv`;

my @sample_nnn;

foreach my $file (@nnn){

#~ scalar
my $sample;
my $org;
my $ligne;
my $entete;

#~ tab
my @entete;
my @org;
my @value;

#~ hash
my %coverage;
my %seq;
my %otu;
my %schao;
my %shannon;

open(FILE, $file);

	while(my $lic=<FILE>){
	
		if($lic=~/^(\S+)\n/){
			$sample=$1;
		}
	
		if($lic=~/(#\w+\t#\w+\t\w+\t\w+\t\w+)/){
			$entete=$1;
			@entete=split("\t", $entete);
		}
	
		if(($lic=~/^\S+\t/) || ($lic=~/^\w+\s\w+\s\w+\t/)){
		
			$ligne=$lic;
		
			@value=split("\t", $ligne);
		
			if(($value[1]ne"NA")&&($value[2]ne"NA")&&($value[3]ne"NA")&&($value[4]ne"NA")&&($value[5]ne"NA")
			&&($value[1]!=0)){		
				push(@org, $value[0]);
				my $org=$value[0];
				$seq{$org}=$value[1];
				$otu{$org}=$value[2];
				$schao{$org}=$value[3];
				$shannon{$org}=$value[4];
				$coverage{$org}=$value[5];
				
			}
		}
	}
	close FILE;


	if(($value[1]ne"NA")){
		open (NNN, ">>$path_results/.CSS/.res/taxo_NN_norm.html");

		print NNN "<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />
			<hr class=\"ligne1\"><br>	
			<table id=\"$sample._nnn\"><caption>$sample</caption>
				<tr><td> </td>";
		
		foreach my $ent (@entete){
			print NNN "<th><center>$ent</center></th>";
		}
		print NNN "</tr>";

		foreach my $org (@org){
			print NNN "<tr>
				<td>$org</td>
				<td width=\"100\"><center>$seq{$org}</center></td>
				<td width=\"100\"><center>$otu{$org}</center></td>
				<td width=\"100\"><center>$schao{$org}</center></td>
				<td width=\"100\"><center>$shannon{$org}</center></td>
				<td width=\"100\"><center>$coverage{$org}</center></td>
			</tr>";
		}

		print NNN "</table>";

		push @sample_nnn, $sample;

	}
}
close NNN;




##### Construction du menu ####

open (SAMP, ">$path_results/.CSS/.res/sample.html");

print SAMP "
<script type=\"text/javascript\" src=\"../functions.js\"></script>
<link rel=\"stylesheet\" href=\"../panam.css\" type=\"text/css\" />

<div id= \"menuE\">				
					
	<div class=\"main\">
	<br/>

		<li><a name='lvl1' href=\"#parametre\" lvl='1'> Parameters</a></li><br>
			
		<li><a name='lvl1' href=\"#quality\" lvl='1'> Quality</a></li><br>
		<li><a name='lvl1' href=\"#clust\" lvl='1'> Clustering</a></li><br>
						
		<li><a name='lvl_1' href=\"javascript:visibilite('lvl_1')\" lvl='1'> &#8250; &alpha; and &beta; diversities
			
			<ul id='lvl_1'>
				<li lvl='2' style=\"display:none;\"><a href=\"javascript:visibilite('lvl_2');\">&#8250; &alpha; diversity </a>
					
				<ul id='lvl_2'>
					<li lvl='3' style=\"display:none;\"><a href=\"#rare\"> Rarefaction curves</a></li>
					<li lvl='3' style=\"display:none;\"><a href=\"#rich\"> Richness and diversity indices</a></li>
					<li lvl='3' style=\"display:none;\"><a href=\"#rich_n\"> Richness and diversity indices normalized</a></li>";

					if (-d "$path_results/PhyloDiv_output"){
						print SAMP "<li lvl='3' style=\"display:none;\"><a href=\"#indices\"> Phylogenetic indices</a></li>";
					}
				print SAMP "</ul>

				<li lvl='2' style=\"display:none;\"><a href=\"javascript:visibilite('lvl2');\">&#8250; &beta; diversity </a>
					
				<ul id='lvl2'>";
					if (-d "$path_results/PhyloDiv_output"){
						print SAMP "<li lvl='3' style=\"display:none;\"><a href=\"#phyletic\"> By phyletic level</a></li>";
					}
					print SAMP "<li lvl='3' style=\"display:none;\"><a href=\"#otu\"> By OTU level</a></li>
				</ul>
		</ul>

		<a name='lvl_1' href=\"javascript:visibilite('lvl1')\" lvl='1'>&#8250; Taxonomy
				
				<ul id='lvl1'>
					<li lvl='2' style=\"display:none;\"><a href=\"#sim\"> Similarity</a></li>
					<li lvl='2' style=\"display:none;\"><a href=\"javascript:visibilite('lvl.2');\">&#8250; LCA</a> 
					
					<ul id='lvl.2'>
						<li lvl='3' style=\"display:none;\"><a href=\"javascript:visibilite('lvl3');\">&#8250; Non-normalized</a>

						<ul id='lvl3'>
							<li lvl='4' style=\"display:none;\"><a href=\"#krona_lca\">Krona</a></li>";		
				
							for (my $i=0; $i<scalar(@sample_lca); $i++){
								print SAMP "<li lvl='4' style=\"display:none;\"><a href=\"#$sample_lca[$i]._lca\">$sample_lca[$i]</a></li>";				
							}

							print SAMP "
						</ul>

						</li>
						<li lvl='3'  style=\"display:none;\"><a href=\"javascript:visibilite('lvl_3');\">&#8250; Normalized</a>
						
						<ul id='lvl_3'>
							<li lvl='4' style=\"display:none;\"><a href=\"#krona_lcan\">Krona</a></li>";
										
							for (my $i=0; $i<scalar(@sample_lcan); $i++){
								print SAMP "<li lvl='4' style=\"display:none;\"><a href=\"#$sample_lcan[$i]._lcan\">$sample_lcan[$i]</a></li>";				
							}

							print SAMP "
						</ul>
						</li>
					</ul>
					</li>
					
					<li lvl='2' style=\"display:none;\"><a href=\"javascript:visibilite('lvl,2');\">&#8250; NN</a> 
					
					<ul id='lvl,2'>
						<li lvl='3' style=\"display:none;\"><a href=\"javascript:visibilite('lvl,3');\">&#8250; Non-normalized</a>
						<ul id='lvl,3'>
							<li lvl='4' style=\"display:none;\"><a href=\"#krona_nn\">Krona</a></li>";			
				
							for (my $i=0; $i<scalar(@sample_nn); $i++){
								print SAMP "<li lvl='4' style=\"display:none;\"><a href=\"#$sample_nn[$i]._nn\">$sample_nn[$i]</a></li>";				
							}

							print SAMP "
						</ul>
						</li>
						<li lvl='3' style=\"display:none;\"><a href=\"javascript:visibilite('lvl34');\">&#8250; Normalized</a>

						<ul id='lvl34'>
							<li lvl='4' style=\"display:none;\"><a href=\"#krona_nnn\">Krona</a></li>";
							
							for (my $i=0; $i<scalar(@sample_nnn); $i++){
								print SAMP "<li lvl='4' style=\"display:none;\"><a href=\"#$sample_nnn[$i]._nnn\">$sample_nnn[$i]</a></li>";				
							}

							print SAMP "
						</ul>
						</li>
					</ul>
					</li>
				</ul>
				</li>
			</ul>
			</a></li>
		</div>
				
</div>";


if($seuil>1000){
	print PARA "<tr bgcolor=\"#a9ccfb\">
		<td>Threshold for normalizing sequences</td>
		<td> $seuil</td>
	</tr>
	</table></center></body></html>";
}
else {
	print PARA "<tr bgcolor=\"#a9ccfb\">
		<td>Seuil</td>
		<td> <1000 &#8594; Threshold too low</td>
	</tr>
	</table></center></body></html>";
}

close PARA;
