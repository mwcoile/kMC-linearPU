# archives current version of code as a previous version that can be easily restored or viewed

CurrentDate=$(date +"%y%m%d")
Comment="backup"
Filename1="molecular_weight.h" 
Filename4="molecular_weight.cpp"
Filename2="chain_update.h"
Filename5="chain_update.cpp"
Filename3="main_KMC_PU.cpp"
mkdir "./OldVersions/${CurrentDate}$Comment"
path="./OldVersions/${CurrentDate}$Comment"
cp "$Filename1" "${path}/${CurrentDate}$Filename1"
cp "$Filename2" "${path}/${CurrentDate}$Filename2"
cp "$Filename3" "${path}/${CurrentDate}$Filename3"
cp "$Filename4" "${path}/${CurrentDate}$Filename4"
cp "$Filename5" "${path}/${CurrentDate}$Filename5"

