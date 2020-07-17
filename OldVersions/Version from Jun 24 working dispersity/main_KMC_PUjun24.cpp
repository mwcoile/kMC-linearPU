#include<iomanip>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include"chain_updatejun24.h"
#include"molecular_weightjun24.h"

// author: Matthew

// This code is designed to generate polyurethane sequences
int main() {
    // inputs //
    double A_kf = 371.5214/(60*60); double Ea_kf = 24.6; // A_kf units: L/mol h. Ea units: kJ/mol (kf = rate constant for forward reaction 1)
    const double Na = 6.02E23;
    // for replicating Krol
    int num_of_molecules_in_simulation=1E5;
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
    // IGNORE THIS COMMENT need to randomly select whether A-B addition or B-A addition occurs.... then need to randomly select the monomer and the side of the monomer that is added to..
    int moleculesA = num_of_molecules_in_simulation; // bifunctional monomer A
    int moleculesB = num_of_molecules_in_simulation; // bifunctional monomer B
    // moleculeA = type 0, moleculeB = type 1
    double monomermassA[1]={1}; // g/mol
    double monomermassB[1]={1}; // g/mol
    double dispersity;
    double Mn;
    double Mw; 
    double over_x; // conversion with respect to A functional groups!
    // also need to keep track of the ends of the chain
    int chainsA = 0; // number of chain ends which are an "A"
    int chainsB = 0; // number of chain ends which are a "B"
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
    molwt.open("molecular_weight2.txt"); // store time, Mn, Mw, Dispersity
    molwt << "time           Conversion     Mn         Mw         D          total rate     unreactedfunctional\n";
    double kf= A_kf*exp(-Ea_kf/(R*T)); // adding a to b
    // ideally rewrite this such that difference in reactivity between A adding to B and B adding to A is taken into account
    
    std::srand(std::time(0));
    double time = 0; // seconds
    double simulation_time = 60*60*36*10; // [=] seconds (30 hours in seconds)

    // KMC loop
    while (time<simulation_time) {

        // calculate propensity functions
        double c = kf/V; // 1 / (mol s)
        // update total rate
        double total_rate =c*(2*moleculesA+chainsA)*(2*moleculesB+chainsB)/Na; // molecules/s 2*molecules because each monomer is bifunctional
        // choose timestep tau
        double r1 = 1.0*std::rand()/RAND_MAX; 
        double tau = (1/total_rate)*log(1/r1); // calculate timestep tau
        // choose reaction to take place
        double r2 = 1.0*std::rand()/RAND_MAX;
        double reaction = r2*total_rate;
        // pick which A, B
        double r3 = 1.0*std::rand()/RAND_MAX;
        double whichA = r3*(chainsA+2*moleculesA);
        double r4 = 1.0*std::rand()/RAND_MAX;
        double whichB = r4*(chainsB+2*moleculesB); 
        // update chains
        explicit_sequence_record(whichA,whichB,moleculesA,moleculesB,chainsA,chainsB,all_chains,loops);
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
                F << std::setw(6) << (moleculesA+moleculesB)/(Na*V) << "     ";
            }
            
        }
        conc << std::left << std::setw(10) << time << "     " << std::setw(6) << moleculesA/(Na*V) << "     " << std::setw(6) << moleculesB/(Na*V) << "     " << std::setw(6) << (all_chains.size()+loops.size())/(Na*V) << "\n";
        
        // Record dispersity here at some point
        molecular_weight(Mn, Mw, all_chains, loops, moleculesA, moleculesB, monomermassA, monomermassB);
        dispersity=Mw/Mn; // calculate polydispersity index PDI
        over_x=1-(chainsA+2*moleculesA)/(2*moleculesA0); // calculate overall conversion of A functional group
        int unreactedfunctionalgroups;
        unreactedfunctionalgroups=(moleculesA*2+moleculesB*2+chainsA+chainsB);
        molwt << std::left << std::setw(10) << time << "     " << std::setw(10) << over_x << "     " << std::setw(6) << Mn << "     " << std::setw(6) << Mw << "     " << std::setw(6) << dispersity << "     " << std::setw(10) << total_rate << "     " << std::setw(6) << unreactedfunctionalgroups << "\n";
    }

    conc.close();
    F.close();
    molwt.close();
    // total number of events that can occur is 1 for each reaction 
    return 0;
}