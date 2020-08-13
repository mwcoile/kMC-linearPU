#pragma once
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include"chain_update.h" // I don't understand why I have to include this line to get chain_pool recognized, given that the main cpp file contains an include"chain update" statement before inclding this file... I thought include was a copy/paste directive, so I don't see why it doens't work here...

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


/*
This code can be used if it is desired to include unreacted monomers in the molecular weight calculation. 
void initialize_molecular_weight(double& Mn, double& Mw,std::vector<int>& monomerA, std::vector<int>& monomerB, std::vector<double>& monomermassA,std::vector<double>& monomermassB, double& sumNi, double& sumMiNi, double& sumMi2Ni) {
    // this function is for the special case at t=0 where only monomers (no chains) are present in the system
    // it calculates initial Mn, Mw, and initializes sumNi, sumMiNi, and sumMi2Ni 
    double Mi_A=0; double Ni_A=0; // Mi_A where A stands for monomer A; Ni_A where A stands for monomer A
    for (int i=0;i<monomerA.size();i++){
        Mi_A=monomermassA[i];
        Ni_A=monomerA[i];
        sumMiNi+=Mi_A*Ni_A;
        sumMi2Ni+=Mi_A*Mi_A*Ni_A;
        sumNi+=Ni_A;
    }
    double Mi_B=0; double Ni_B=0; // Mi_B where B stands for monomer B; Ni_B where B stands for monomer B
    for (int i=0;i<monomerB.size();i++){
        Mi_B=monomermassB[i];
        Ni_B=monomerB[i];
        sumMiNi+=Mi_B*Ni_B;
        sumMi2Ni+=Mi_B*Mi_B*Ni_B;
        sumNi+=Ni_B;
    }
    Mn = sumMiNi/sumNi; // Mn at time t=0, before any chains or loops have formed
    Mw = sumMi2Ni/sumMiNi; // Mw at time t=0, before any chains or loops have formed
}
*/