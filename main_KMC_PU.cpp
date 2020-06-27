#include<iomanip>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include<array>
#include"chain_update.h"
#include"molecular_weight.h"

// author: Matthew
void specificfunctionalgroup(int &A_index, int &B_index, std::array<int,2> moleculesA,std::array<int,1> moleculesB, int mu) {
    int count=0;
    for (A_index=0;A_index<moleculesA.size();A_index++) {
        for (B_index=0;B_index<moleculesB.size();B_index++) {
            if (count>=mu) {
                return;
            }
            count++;
        }
    }
}

// This code is designed to generate polyurethane sequences
int main() {
    // inputs //
    double A_tomita = 371.5214/(60*60); double Ea_tomita = 24.6; // L/mol h
    const int M = 2; // number of reaction channels considered
    double A_kf[M] = {A_tomita,A_tomita}; double Ea_kf[M] = {Ea_tomita,Ea_tomita}; // A_kf units: L/mol h. Ea units: kJ/mol (kf = rate constant for forward reaction 1)
    const double Na = 6.02E23;
    // for replicating Krol
    int num_of_molecules_in_simulation=1E4;
    const double moleculesA0=num_of_molecules_in_simulation;
    double V = num_of_molecules_in_simulation/Na; // for 1 mol/L the volume is L (which is also dm^3)
    //double A_kba = 100; double Ea_kr = 10; // reverse reaction 1
    // I think I will have to get the number of monomers from an input file, as well as the number of different reactions
    // for 3 monomers there are 3+2+1 different reactions: 6
    // however for AA BB polymerization, 3 monomers only permits 2 reactions (well 4 if you count forward and reverse separately)
    // In general, the number of possible reactions is given by # of cyclic carbonate monomers * # of amine monomers (ignoring reverse reactions)
    // Can we assume depolymerization is negligible.. only during polymerization step... will need explicit depolymerization for recycling however :( 
    // But that would be in presence of a solvent... hmm so the solvent could attack any chain, anywhere (if you ignore the fact that realistically it could only attack at surface...)
    // For solvolysis, could we ignore the polymerization step? If so, I think this should be fairly straightforward. If not, I will need to include a balance between the
    // polymerization and depolymerization steps.
    const double R = 0.008314; // ideal gas constant, kJ/(K mol)
    double T = 273+80; // temperature (K)
    std::array<int, 2> moleculesA = {num_of_molecules_in_simulation/2,num_of_molecules_in_simulation/2}; // bifunctional monomer containing functional group A
    std::array<int, 1> moleculesB = {num_of_molecules_in_simulation}; // bifunctional monomer B
    // moleculeA = type 0, moleculeB = type 1
    std::array<double,2> monomermassA={1.0,1.0}; // g/mol
    std::array<double,1> monomermassB={1.0}; // g/mol
    double dispersity;
    double Mn;
    double Mw; 
    double over_x; // conversion with respect to A functional groups!
    // also need to keep track of the ends of the chain
    std::array<int, 2> chainsA = {}; // number of chain ends which are an "A", initialize to 0
    std::array<int, 1> chainsB = {}; // number of chain ends which are a "B", initialize to 0
    chain_pool all_chains;  
    chain_pool loops;
    // For bifunctional monomers, the first dyad in a chain
    // will always have frontend = 0 and backend = 1.
    std::ofstream conc;
    conc.open("concentrations.txt");
    conc << "time           AA         BB         polymer \n";
    std::ofstream F;
    F.open("F_tomita2.txt");
    F << "time           F1         F2         F3         F4         F5         F6         F7         F8         F9n\n";
    std::ofstream molwt;
    molwt.open("molecular_weight.txt"); // store time, Mn, Mw, Dispersity
    molwt << "time           Conversion     Mn         Mw         D\n";
    std::ofstream debug;
    debug.open("debugging.txt");
    debug << "time           total rate     Unreacted functional groups\n";
    // calculate all rate constants from kinetic parameters and temperature
    double kf[M];
    for (int i=0;i<M;i++) {
        kf[i]=A_kf[i]*exp(-Ea_kf[i]/(R*T)); // adding a to b
    }
    // ideally rewrite this such that difference in reactivity between A adding to B and B adding to A is taken into account
    
    std::srand(std::time(0));
    double time = 0; // seconds
    double simulation_time = 60*60*36*10; // [=] seconds 

    // KMC loop
    while (time<simulation_time) {

        // calculate propensity functions
        double c[M];
        for (int i=0;i<M;i++) {
            c[i] = kf[i]/V; // 1 / (mol s)
        }
        // update total rate
        double total_rate=0;
        // note: first entries in total rate are moleculesA.size()*moleculesB.size() looping numbering reactions as R1 = A1B1, R2 = A2B1, R3 = A3B1... Rn=A2B1,Rn+1=A2B2,Rn+2=A2B3...
        int v=0;
        int h[M]={};
        double Rv[M]={};
        for (int j=0;j<moleculesA.size();j++) {
            for (int k=0;k<moleculesB.size();k++) {
                h[v]=(2*moleculesA[j]+chainsA[j])*(2*moleculesB[k]+chainsB[k]); // hv is the product of the numbers of molecular reactants involved in the vth reaction channel present at time t
                Rv[v]=c[v]*h[v]/Na; // molecules/s 2*molecules because each monomer is bifunctional
                total_rate+=Rv[v]; // molecules/s 2*molecules because each monomer is bifunctional
                v+=1;
            }
        }
        //double total_rate =c*(2*moleculesA+chainsA)*(2*moleculesB+chainsB)/Na; 
        // choose timestep tau
        double r1 = 1.0*std::rand()/RAND_MAX; 
        double tau = (1/total_rate)*log(1/r1); // calculate timestep tau
        // choose reaction to take place using Gillespie algorithm

        double r2 = 1.0*std::rand()/RAND_MAX;
        //double reaction = r2*total_rate;
        int mu=0;
        double sumRv=0; // Lin Wang eqn 1 multiplied by total rate
        while (sumRv<r2*total_rate) { // this is probably better implemented in a do loop
            sumRv+=Rv[mu];
            if (sumRv<r2*total_rate) mu++;
        }
        // Translation from reaction channel mu to 
        // specific AA monomer type and BB monomer type 
        int A_index=0;
        int B_index=0;
        int count=0;
        specificfunctionalgroup(A_index, B_index, moleculesA, moleculesB, mu);
        
        // pick which A, B functional group
        double r3 = 1.0*std::rand()/RAND_MAX;
        double whichA = r3*(chainsA[A_index]+2*moleculesA[A_index]);
        double r4 = 1.0*std::rand()/RAND_MAX;
        double whichB = r4*(chainsB[B_index]+2*moleculesB[B_index]); 
        // update chains
        explicit_sequence_record(whichA,whichB,moleculesA,moleculesB,A_index,B_index,chainsA,chainsB,all_chains,loops);
        time += tau;
        
        // END KMC CALCULATIONS. 

        // RECORD SEQUENCE STARTS
        // record fractions of length 2, 3, 4 ... eventually add through 8, and dispersity
        int F_krol [10] = {0}; // vector to track fractions to reproduce Krol and Gawdzik

        for (int i=0;i<all_chains.size();i++){
            if (all_chains[i].v.size()<10){
                F_krol[all_chains[i].v.size()-1]++;
            }
            else if (all_chains[i].v.size()>=10) {
                F_krol[9]++;
            }
        }
        for (int i=0;i<loops.size();i++){
            if (loops[i].v.size()<10){
                F_krol[loops[i].v.size()]++;
            }
            else if (loops[i].v.size()>=10){
                F_krol[9]++;
            }
        }
        F << std::left << std::setw(10) << time << "     ";
        // unfortunately the below hard codes in the length of F_krol "i<10"
        for (int i=0;i<10;i++) {
            if (i!=0){
                F << std::setw(6) << F_krol[i]/(Na*V) << "     ";
                if (i==9) F << "\n";
            }
            else {
                int total_monomer = 0; // calculate total amount of monomer
                // add up total number of A molecules
                for (int j=0;j<moleculesA.size();j++) {
                    total_monomer+=moleculesA[j];
                }
                for (int j=0;j<moleculesB.size();j++) {
                    total_monomer+=moleculesB[j];
                }
                F << std::setw(6) << (total_monomer)/(Na*V) << "     ";
            }
            
        }
        conc << std::left << std::setw(10) << time << "     " << std::setw(6) << moleculesA.size()/(Na*V) << "     " << std::setw(6) << moleculesB.size()/(Na*V) << "     " << std::setw(6) << (all_chains.size()+loops.size())/(Na*V) << "\n";
        
        // Record dispersity here at some point
        molecular_weight(Mn, Mw, all_chains, loops, moleculesA, moleculesB, A_index, B_index, monomermassA, monomermassB);
        dispersity=Mw/Mn; // calculate polydispersity index PDI
        // the below overall conversion line could be improved
        over_x=1-(chainsA[0]+2*moleculesA[0]+chainsA[1]+2*moleculesA[1])/(2*moleculesA0); // calculate overall conversion of A functional group
        int unreactedfunctionalgroups=0;
        for (int i=0; i<moleculesA.size();i++) {
            unreactedfunctionalgroups+=2*moleculesA[i]+chainsA[i];
        }
        for (int i=0; i<moleculesB.size();i++) {
            unreactedfunctionalgroups+=2*moleculesB[i]+chainsB[i];
        }
        molwt << std::left << std::setw(10) << time << "     " << std::setw(10) << over_x << "     " << std::setw(6) << Mn << "     " << std::setw(6) << Mw << "     " << std::setw(6) << dispersity << "     " << std::setw(10) << total_rate << "     " << std::setw(6) << unreactedfunctionalgroups << "\n";
        //debug << std::left << std::setw(10) << time << 
        }

    conc.close();
    F.close();
    molwt.close();
    debug.close();
    // could easily print all sequences here too if desired

    // total number of events that can occur is 1 for each reaction 
    return 0;
}