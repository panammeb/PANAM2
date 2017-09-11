#!/usr/bin/perl -w

use strict;
use warnings;

use FindBin;
use Getopt::Long;
use Cwd 'abs_path';

# $path=abs_path($0);
my $path_panam="$FindBin::RealBin/";

my $directory;
my $file1;my $file2;
my $forward; my $reverse;
my $tags;
my $ident;
my $length;
my @domain;
my $path_results;
my $ini;
my $help;
my $panam_ini;
my $div;
my $version;


## Give name to option
GetOptions ('d=s' => \$directory, 'r1=s' => \$file1,  'r2=s' => \$file2, 'f=s' => \$forward,  'r=s' => \$reverse, 'tags=s' => \$tags, 'id=f' => \$ident, 'length=f' => \$length, 'dom=s' => \@domain, 'o=s' => \$path_results,'div=s' => \$div, 'ini=s' => \$ini, 'h!' => \$help, 'v!' => \$version);



if (defined($help)){
# 	if (!defined($ini) || defined($ini)){
	
	print "Usage: perl panam.pl [options]
	perl panam2.pl -h
	perl panam2.pl -v PANAM2 version 
	
	example: perl  panam2.pl -o example -r1 test/R1.fastq.gz -r2 test/R2.fastq.gz -f GTGYCAGCMGCCGCGGTA -r CCCCGYCAATTCMTTTRAGT -tags test/tags.txt -id 0.97 -dom bacteria -dom archaea

Mandatory: 
	-r1\t\t fastq file which contains amorces
	-r2\t\t fastq file which contains amorces
        -tags\t\t file with sample names and tag sequences
        OR
	-d\t\t path of directory with demultiplexed files


	-f\t\t forward primer
	-r\t\t reverse primer	
	-id\t\t clustering threshold for building OTUs
	-dom\t\t microbial domains analyzed

\n
        Optionnal: 
	-div\t\t Y or N (default) phylogenetic indices (R and packages have to be installedin your linux system)
        -length\t\t maximum length of amplicons (default: from 200 to 550bp)
        -o\t\t path for results	(default: Results)
    
OR
	-ini\t\t file .ini with all parameters described above and others explained in the .ini file located in the PANAM2 directory\n";
# 	}
}

elsif(defined($version)){
print "PANAM2 version: ";
system("cat $path_panam/.version");


}

	
elsif (!defined($ini)) { ## 


				if(!defined($file1) && !defined($file2) && !defined($directory)) {
					die ("\nDefine fastq files (-r1 -r2) OR a directory with fastq files (-d)\n");
				  }
				  

				  

				 if(!defined($forward) || !defined($reverse)) {
					die ("\nPrimers not defined\n");
				  }

				  if(!defined($directory) && !defined($tags)) {
				   die ("\nTag file not defined\n");
				  }
			
				if(!defined($domain[0])) {
					die ("\nmicrobial domain not defined\n");	
				}
				  

				  if(!defined($length)) {
					$length=550;
				  }

				  if(!defined($path_results)) {
					$path_results="Results";
				  }

				  if(!defined($div)) {
					$div="N";
				  }

				  if (!defined($ident) || $ident >= 1) {
					die ("\n-id must be defined with a value smaller than 1.\n");
				  }

				if (!(-d $path_results)){
					`mkdir $path_results`;
				}

				
## Create and write in panam.ini 
open (INI, ">".$path_results."/panam.ini");
				
				
                                        print INI "PATH_RESULTS\t$path_results\n";
                                        
                                        if(!defined($directory)){
                                        print INI "INPUT_FILE_FORWARD\t$file1\n";
                                        print INI "INPUT_FILE_REVERSE\t$file2\n";
                                        print INI "BARCODE_FILE\t$tags\n";
                                        }
                                        else
                                        {
                                        print INI "DEMUL_FOLDER\t$directory\n";
                                        print INI "INPUT_FILE_FORWARD\n";
                                        print INI "INPUT_FILE_REVERSE\n";
                                        print INI "BARCODE_FILE\n"
                                        }
				  
					print INI "MAX_SEQUENCE_LENGTH\t$length\n";
					print INI "FORWARD_PRIMER_SEQUENCE\t$forward\n";
					print INI "REVERSE_PRIMER_SEQUENCE\t$reverse\n";

					print INI "CLUSTERING_CUTOFF\t$ident\n";
					
					print INI "REFERENCE_BASE\n";
				
					for (my $i=0; $i<3; $i++){
                                            if(defined($domain[$i])){
						print INI "DOMAIN\t$domain[$i]\n"; # 
						}
					}
					
			
				close INI ;
				$panam_ini=$path_results."/panam.ini";
# 				## Lancement des programmes perl
				system("perl ".$path_panam."quality_panam.pl ".$panam_ini);
				system("perl ".$path_panam."preprocess_panam.pl ".$panam_ini);
				system("perl ".$path_panam."taxonomy_panam.pl ".$panam_ini);
				
					if($div eq "Y"){
                                        system("perl ".$path_panam."phylodiv_panam.pl ".$panam_ini);
					}
					
                                system("perl ".$path_panam."postprocess_panam.pl ".$panam_ini);
                                
}

elsif (defined($ini)) {
			
					$panam_ini=$path_results."/panam.ini";

					## Lancement des programmes perl
					system("perl ".$path_panam."quality_panam.pl ".$panam_ini);
					system("perl ".$path_panam."preprocess_panam.pl ".$panam_ini);
					system("perl ".$path_panam."taxonomy_panam.pl ".$panam_ini);
					system("perl ".$path_panam."phylodiv_panam.pl ".$panam_ini);
					system("perl ".$path_panam."postprocess_panam.pl ".$panam_ini);

			
		
	
}













  
