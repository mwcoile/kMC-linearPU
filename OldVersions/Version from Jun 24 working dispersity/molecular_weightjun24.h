#pragma once
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include"chain_updatejun24.h" // I don't understand why I have to include this line to get chain_pool recognized, given that the main cpp file contains an include"chain update" statement before inclding this file... I thought include was a copy/paste directive, so I don't see why it doens't work here...

void molecular_weight(double& Mn, double& Mw, chain_pool& all_chains,chain_pool& loops, int moleculesA, int moleculesB, double* monomermassA,double* monomermassB){
    double sumMiNi=0;
    double sumNi=0;
    double sumMi2Ni=0;
    // Calculate Mi and Ni for each chain, then sum them all up to calculate Mn and Mw
    for (int i=0;i<all_chains.size();i++){
        double Mi=0;
        double Ni=1; // all_chains[i].v.size();

        for (int j=0; j<all_chains[i].v.size();j++) {
            // if monomer is type AA, find the corresponding monomer weight and add it to Mi
            if (all_chains[i].v[j]==0) {
                Mi+=monomermassA[0];
            }
            // if monomer is type BB
            else if (all_chains[i].v[j]==1) {
                Mi+=monomermassB[0];
            }
        }
        sumMiNi+=Mi*Ni;
        sumMi2Ni+=Mi*Mi*Ni;
        sumNi+=Ni;
    }
    for (int i=0;i<loops.size();i++){
        double Mi=0;
        double Ni=1; // loops[i].v.size();

        for (int j=0; j<loops[i].v.size();j++) {
            // if monomer is type AA, find the corresponding monomer weight and add it to Mi
            if (loops[i].v[j]==0) {
                Mi+=monomermassA[0];
            }
            // if monomer is type BB
            else if (isloops[i].v[j]==1) {
                Mi+=monomermassB[0];
            }
        }
        sumMiNi+=Mi*Ni;
        sumMi2Ni+=Mi*Mi*Ni;
        sumNi+=Ni;
    }
    // chains of length 1 are monomers still
    double Mi_A=monomermassA[0]; // Mi_A where A stands for monomer A
    double Ni_A=moleculesA; // Ni_A where A stands for monomer A
    sumMiNi+=Mi_A*Ni_A;
    sumMi2Ni+=Mi_A*Mi_A*Ni_A;
    sumNi+=Ni_A;
    
    double Mi_B=monomermassB[0]; // Mi_A where A stands for monomer A
    double Ni_B=moleculesB; // Ni_A where A stands for monomer A
    sumMiNi+=Mi_B*Ni_B;
    sumMi2Ni+=Mi_B*Mi_B*Ni_B;
    sumNi+=Ni_B;
    Mn=sumMiNi/sumNi;
    Mw=sumMi2Ni/sumMiNi;

    // take 2 do this tomorrow
    
    // sum up all masses, then divide by size
    // for (int i=0;i<all_chains.size();i++) {
    //     double 
    //     double Ni=all_chains[i].v.size();

    //     for (int j=0; j<Ni;j++) {
    //         // if monomer is type AA, find the corresponding monomer weight and add it to Mi
    //         if (all_chains[i].v[j]==0) {
    //             Mi+=monomermassA[0];
    //         }
    //         // if monomer is type BB
    //         else if (all_chains[i].v[j]==1) {
    //             Mi+=monomermassB[0];
    //         }
    //     }
    //     sumMiNi+=Mi*Ni;
    //     sumMi2Ni+=Mi*Mi*Ni;
    //     sumNi+=Ni;
    // }
}