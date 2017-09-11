path=$1
cd $path

for name in `ls taxo*` 
do 

	
	awk -F "\t\t" '{print $1}' $name | grep -v "^\s\w" | sed '1d'  > .phyla.csv


	if [ ! -d .CSS/.taxonomy ]; then
		`mkdir .CSS/.taxonomy`
	fi 


	name2=`echo $name | cut -d. -f 1`

	number=2
	

	for sample in `head -1 $name | tr "\t" "\n" | sed '/^$/d' | cut -d/ -f 2`
	do 
		echo "grep -v \"^\s\w\" $name | awk -F \"\t\t+\" '{print \$$number}' | sed '1d' > $sample.tmp.csv" > .tmp.sh
		bash .tmp.sh
		echo $sample > .CSS/.taxonomy/$name2.$sample.csv
		paste .phyla.csv $sample.tmp.csv >> .CSS/.taxonomy/$name2.$sample.csv
		rm -f $sample.tmp.csv
		number=$(($number+1))
	done

rm -f .phyla.csv
rm -f .tmp.sh

done
