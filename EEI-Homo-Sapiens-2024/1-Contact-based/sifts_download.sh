# WARNING:
# USE THIS SCRIPT ONLY IF SCRIPT 2_download_sifts_mapping.r DOESN'T PROPERLY DOWNLOAD/UNZIP FILES

list=$(find ../data/SIFTS -name "*.xml.gz")
#echo $list
for l in $list
do
#echo $l
file=${l##*/}
name=${file%.gz}
#echo $file
#echo $name
S=${file:1:2}
#echo $S
if [ ! -f ../../public_data/SIFTS/$name ]; then
    wget -O ../../public_data/SIFTS/$file "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/$S/$file"
    gunzip ../../public_data/SIFTS/$file
fi

done
