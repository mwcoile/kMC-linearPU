#pragma once
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include"chain_update.h" // I don't understand why I have to include this line to get chain_pool recognized, given that the main cpp file contains an include"chain update" statement before inclding this file... I thought include was a copy/paste directive, so I don't see why it doens't work here...

void molecular_weight(double& Mn, double& Mw, chain_pool& all_chains,chain_pool& loops, std::array<int,2>& moleculesA, std::array<int,1>& moleculesB, int& A_index, int& B_index, std::array<double,2>& monomermassA,std::array<double,1>& monomermassB){
    double sumMiNi=0;
    double sumNi=0;
    double sumMi2Ni=0;
    // Calculate Mi for each chain, then sum them all up to calculate Mn and Mw
    for (int i=0;i<all_chains.size();i++){
        double Mi=0;
        double Ni=1; // Only 1 chain has this mass so Ni=0

        for (int j=0; j<all_chains[i].v.size();j++) {
            // if monomer is type AA, find the corresponding monomer weight and add it to Mi
            if (all_chains[i].v[j][0]==0) {
                Mi+=monomermassA[all_chains[i].v[j][1]];
            }
            // if monomer is type BB
            else if (all_chains[i].v[j][0]==1) {
                Mi+=monomermassB[all_chains[i].v[j][1]];
            }
        }
        sumMiNi+=Mi*Ni;
        sumMi2Ni+=Mi*Mi*Ni;
        sumNi+=Ni;
    }
    // Do the same for loops
    for (int i=0;i<loops.size();i++){
        double Mi=0;
        double Ni=1; // loops[i].v.size();

        for (int j=0; j<loops[i].v.size();j++) {
            // if monomer is type AA, find the corresponding monomer weight and add it to Mi
            if (loops[i].v[j][0]==0) {
                Mi+=monomermassA[loops[i].v[j][1]];
            }
            // if monomer is type BB
            else if (loops[i].v[j][0]==1) {
                Mi+=monomermassB[loops[i].v[j][1]];
            }
        }
        sumMiNi+=Mi*Ni;
        sumMi2Ni+=Mi*Mi*Ni;
        sumNi+=Ni;
    }
    // chains of length 1 are monomers
    double Mi_A=0; double Ni_A=0; // Mi_A where A stands for monomer A; Ni_A where A stands for monomer A
    for (int i=0;i<moleculesA.size();i++){
        Mi_A=monomermassA[i];
        Ni_A=moleculesA[i];
        sumMiNi+=Mi_A*Ni_A;
        sumMi2Ni+=Mi_A*Mi_A*Ni_A;
        sumNi+=Ni_A;
    }
    double Mi_B=0; double Ni_B=0; // Mi_B where B stands for monomer B; Ni_B where B stands for monomer B
    for (int i=0;i<moleculesB.size();i++){
        Mi_B=monomermassB[i];
        Ni_B=moleculesB[i];
        sumMiNi+=Mi_B*Ni_B;
        sumMi2Ni+=Mi_B*Mi_B*Ni_B;
        sumNi+=Ni_B;
    }
    Mn=sumMiNi/sumNi;
    Mw=sumMi2Ni/sumMiNi;
}