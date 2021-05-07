#include"molecular_weight.h"
#include"chain_update.h"

void molecular_weight(double& Mn, double& Mw, chain_pool& all_chains, chain_pool& loops, std::vector<int>& monomerA, std::vector<int>& monomerB, std::vector<double>& monomermassA, std::vector<double>& monomermassB, bool& isloop, bool& isnewchain, bool& ismonomerA,bool& ismonomerB,double& sumNi, double& sumMiNi, double& sumMi2Ni, double& Mi_A, double& Mi_B){
   if (isloop == false) {
        if (isnewchain == true) {
            sumNi++;
            sumMiNi+=Mi_A+Mi_B;
            sumMi2Ni+=(Mi_A+Mi_B)*(Mi_A+Mi_B);
        }
        else if (ismonomerA==true){
            // sumNi does not change
            sumMiNi+=Mi_A; // add weight of additional monomer
            sumMi2Ni+=(Mi_A+Mi_B)*(Mi_A+Mi_B)-Mi_B*Mi_B; // newchainweight^2-oldchainweight^2
        }
        else if (ismonomerB==true){
            // sumNi does not change
            sumMiNi+=Mi_B; // add weight of additional monomer
            sumMi2Ni+=(Mi_A+Mi_B)*(Mi_A+Mi_B)-Mi_A*Mi_A; // newchainweight^2-oldchainweight^2
        } 
        else {
            sumNi--;
            // sumMiNi does not change
            sumMi2Ni+=(Mi_A+Mi_B)*(Mi_A+Mi_B)-Mi_A*Mi_A-Mi_B*Mi_B;
        }
    }
    isloop = false;
    isnewchain=false;
    ismonomerA=false;
    ismonomerB=false;
    Mn=sumMiNi/sumNi;
    Mw=sumMi2Ni/sumMiNi;
}