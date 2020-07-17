# archives current version of code as a previous version that can be easily restored or viewed

CurrentDate=$(date +"%y%m%d")
Filename1="molecular_weight.h" 
Filename2="chain_update.h"
Filename3="main_KMC_PU.cpp"
mkdir "./OldVersions/${CurrentDate}"
path="./OldVersions/${CurrentDate}"
cp "$Filename1" "${path}/${CurrentDate}$Filename1"
cp "$Filename2" "${path}/${CurrentDate}$Filename2"
cp "$Filename3" "${path}/${CurrentDate}$Filename3"

