#!/usr/bin/perl -w
use strict;
use warnings;


my $current_dir=`pwd`;
chomp $current_dir;


chdir "bin";

### Installation des logiciels utils pour PANAM2 ###
install_krona();
install_fasttree();
install_hmmer();
install_vsearch();
install_text() ;
# install_bioperl(); # Must be installed 




### Verification de l'existance de R et de ses packages ####
if (-e "/usr/bin/R"){


my $package=`Rscript -e "list.of.packages <- c('vegan', 'picante', 'cluster','MASS', 'phyloseq' ); new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]; print(length(new.packages))"`;
chomp $package;


    if ($package !~ /\[1\] 0/ ) 
    { 
    print "Some Packages among (vegan, picante, cluster,MASS, phyloseq) were not found and some results vill not be achevied\n";
    }
    else
    {
    print "\n\tR and required packages are installed\n";
    }

}
else{
	print "\n\nR not found in your system\n";
}



chdir ("../");
system("tar -xjf bd_ssrna.tar.bz2");
system ("wget https://www.dropbox.com/s/2y9p1jlfy00bsjq/test.tar.gz");
system ("tar xzf test.tar.gz");



    if  (-z "bin/FastTree") 
    {
    die "\tFastTree has not been installed $!\n" ; 
    }
    elsif ((-z "bin/hmmer-2.3.2/src/hmmalign") and (-z "bin/hmmer-2.3.2/src/hmmbuild")) 
    {
    die "\tHMMER has not been installed $!\n" ; 
    }
    else{
    print "\tPANAM has been successfully installed\n\n";
    system("touch $current_dir/.install");
    }



###############################

sub install_krona{
	print "\nInstall Krona\n";	
	system("wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar");
	system("tar -xf KronaTools-2.7.tar");
 	system("rm KronaTools-2.7.tar");

}

sub install_fasttree{
	print "\nDownload FastTree ...\n";
	system ("wget http://www.microbesonline.org/fasttree/FastTree.c");
	print "\nCompiling FastTree ...\n";
	system ("gcc FastTree.c -DNO_SSE -lm -O3 -finline-functions -funroll-loops -Wall -o FastTree");


}

sub install_hmmer{
	print "\nInstall HMMER v2.3.2 ...\n";	
	system("tar -xzf hmmer-2.3.2.tar.gz");
	chdir ("hmmer-2.3.2");
	system ("./configure");
	system ("make");
	chdir ("../");
}


sub install_vsearch{
	print "\nInstall VSEARCH\n";	
	system("tar xzf v1.9.9.tar.gz");
	chdir("vsearch-1.9.9");
	system("./configure");
        system("make");
        system("mv bin/vsearch ../");
        chdir ("../");
}


sub install_text {
	print "\nInstalling text-csv module ...\n";
	system ("tar -xvzf Text-CSV_XS-1.04.tgz");
	chdir ("Text-CSV_XS-1.04");
	system ("perl Makefile.PL PREFIX=$current_dir/bin/Text-CSV_XS-1.04");
	system ("make test");
	system ("make install");
	chdir ("../");
}

sub install_bioperl {
	system ("unzip bioperl-1.5.2_102.zip");
	chdir "bioperl-1.5.2_102" ;
        system ("perl Build.PL --install_base $current_dir/bin/bioperl-1.5.2_102"); 
	system ("./Build test");
	system ("./Build install");
        chdir ("../");
}


