#include<iomanip>
#include<stdlib.h>
#include<iostream>
#include<string>
#include<vector>
#include<math.h>

// author: Matthew

// This code is designed to generate polyurethane sequences
int main() {
    // inputs //
    double A_kf = 100; double Ea_kf = 10; // forward reaction 1
    double V = 1;
    //double A_kba = 100; double Ea_kr = 10; // reverse reaction 1
    // I think I will have to get the number of monomers from an input file, as well as the number of different reactions
    // for 3 monomers there are 3+2+1 different reactions: 6
    // however for AA BB polymerization, 3 monomers only permits 2 reactions (well 4 if you count forward and reverse separately)
    // In general, the number of possible reactions is given by # of cyclic carbonate monomers * # of amine monomers (ignoring reverse reactions)
    // Can we assume depolymerization is negligible.. only during polymerization step... will need explicit depolymerization for recycling however :( 
    // But that would be in presence of a solvent... hmm so the solvent could attack any chain, anywhere (if you ignore the fact that realistically it could only attack at surface...)
    // For solvolysis, could we ignore the polymerization step? If so, I think this should be fairly straightforward. If not, I will need to include a balance between the
    // polymerization and depolymerization steps.
    double R = 8.314; // ideal gas constant
    double T = 298; // temperature (K)
    //double simulation_time = 10; // total simulation time
    // For the special case with 2 monomers
    // need to randomly select whether A-B addition or B-A addition occurs.... then need to randomly select the monomer and the side of the monomer that is added to..
    int moleculesA = 1000; // bifunctional monomer A
    int moleculesB = 1000; // bifunctional monomer B
    // moleculeA = type 1, moleculeB = type 2

    struct chain
    {
        std::vector<int> v;
        int end1;
        int end2;
    };
    // Need to keep track of all the polymer chains that have been created
    typedef std::vector<chain> chain_pool;
    typedef std::vector<int> test123;
    // also need to keep track of the ends of the chain
    int chainsA = 0; // number of chain ends which are an "A"
    int chainsB = 0; // number of chain ends which are a "B"
    
    //.. could instead write a loop which goes through chainpool and interrogates 
    // each chain for its ends. However, the above step saves that whole loop
    // as long as properly handled. For difunctional monomers, the first dyad in a chain
    // will always have end1 = 1 and end2 = 0.

    // KMC step
    double kf=A_kf*exp(Ea_kf/(R*T)); // adding a to b
    //double kba=A_kab*exp() 
    double c = kf/V;
    // ideally rewrite this such that difference in reactivity between A adding to B and B adding to A is taken into account
    double total_rate =c*(2*moleculesA+chainsA)*(2*moleculesB+chainsB); // 2*molecules because each monomer is bifunctional
    
    // I should rewrite this in terms of reaction parameter c rather than k
    std::srand(std::time(0));
    double r1 = 1.0*std::rand()/RAND_MAX; 
    double tau = (1/total_rate)*log(1/r1); // calculate tau

    // pick which reaction to have occur
    double r2 = std::rand()/RAND_MAX;
    double reaction = r2*total_rate;

    // if (1==1) //(reaction < moleculesB)
    //     // react B
    //     moleculesB-=1;

        
        // update chain end
    // pick which chain to react...
    // double r3 = std::rand()/RAND_MAX;

    double time = 0;
    double simulation_time = tau*1000;
    chain_pool all_chains;

    // KMC loop
    while (time<simulation_time) {
        // first, case of constant concentrations

        // calculate propensity functions

        // update total rate

        // choose timestep tau
        double r1 = 1.0*std::rand()/RAND_MAX; 
        double tau = (1/total_rate)*log(1/r1); // calculate timestep tau

        // choose reaction to take place
        // initial case: all reactions are equally likely
        //if (1==1) {
        // case for which two monomers react to form new chain
        chain newchain;
        newchain.end1 = 1;
        newchain.end2 = 2;
        newchain.v.push_back(newchain.end1);
        newchain.v.push_back(newchain.end2);
        all_chains.push_back(newchain);
        moleculesA-=1;
        moleculesB-=1;
        //}
        time += tau;
    }
    // total number of events that can occur are 1 for each reaction 
    return 0;
}