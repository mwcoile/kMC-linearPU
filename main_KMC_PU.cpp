#include<iomanip>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include<array>
#include<limits>
#include<numeric>
#include"chain_update.h"
#include"molecular_weight.h"

// author: Matthew

using namespace std::chrono;

void specificmonomertype(int &A_monomer_type, int &B_monomer_type, std::vector<int> monomerA,std::vector<int> monomerB, int mu) {
    int count=0;
    for (A_monomer_type=0;A_monomer_type<monomerA.size();A_monomer_type++) {
        for (B_monomer_type=0;B_monomer_type<monomerB.size();B_monomer_type++) {
            if (count>=mu) {
                return;
            }
            count++;
        }
    }
}

// This code is designed to generate polyurethane sequences
int main() {
    auto start = high_resolution_clock::now(); // time the entire program

    // USER INPUTS //
    
    // 1. Reaction conditions & monomers used
    double T = 273+80; // temperature (K)
    std::vector<double> monomermassA={1700,154.25}; // g/mol
    std::vector<double> concentrationsA={0.7*0.1,0.7*0.9}; // M
    std::vector<double> monomermassB={250}; // g/mol
    std::vector<double> concentrationsB={0.7}; // M
    // should be no need to change this line
    const int M = monomermassA.size()*monomermassB.size(); // number of reaction channels considered
    
    // for volume calculation
    const double total_A_concentration = std::accumulate(concentrationsA.begin(), concentrationsA.end(), decltype(concentrationsA)::value_type(0));
    const double total_B_concentration = std::accumulate(concentrationsB.begin(), concentrationsB.end(), decltype(concentrationsB)::value_type(0));
    double total_monomer_concentration=total_A_concentration+total_B_concentration;

    // 2. Kinetic parameters
    // The parameters below are for primary amine+cyclic carbonate in a DMAc solvent (N,N-dimethylacetamide)
    // Source: Tomita 2001. DOI 10.1002/1099-0518(20010101)39:1<162::AID-POLA180>3.0.CO;2-O
    double A_tomita = 371.5214/(60*60); // L/mol s; 
    double Ea_tomita = 24.6; // kJ/mol. 
    // The below vectors of kinetic parameters should be in the order A1+B1, A1+B2, A1+B3... A1+BN, A2+B1,A2+B2, A2+B3... etc.
    // After parameters for all possible combinations of A+B monomer types are inserted, other reactions follow
    std::vector<double> A_kf(M,A_tomita); // A_kf units: L/mol h.
    std::vector<double> Ea_kf(M,Ea_tomita); // Ea units: kJ/mol (kf = rate constant for forward reaction 1)
    std::string filename = "DevelopingMasterBranch_tirrellparameters_200k.txt"; // description to be appended to filenames

    // 3. Simulation details
    double simulation_time = 60*60*(36+24+100); // [=] seconds 36 h + 24 h 
    int num_of_molecules_in_simulation=200000; // number of monomers simulated
    std::vector<int> monomerA; // starting number of bifunctional monomers containing functional group A
    std::vector<int> monomerB; // starting number of bifunctional monomer containing functional group B
    // moleculeA = type 0, moleculeB = type 1
    for (int i=0; i<concentrationsA.size();i++){
        monomerA.push_back(std::round(1.0*num_of_molecules_in_simulation*concentrationsA[i]/total_monomer_concentration));
    }
    for (int i=0; i<concentrationsB.size();i++){
        monomerB.push_back(std::round(1.0*num_of_molecules_in_simulation*concentrationsB[i]/total_monomer_concentration));
    }

    // for Tirrell calculation
    int startingA0 = monomerA[0];
    int startingA1 = monomerA[1];
    // end Tirrell stuff

    // END USER SUPPLIED INPUTS // 

    // declaring variables
    const double Na = 6.02E23; // avogadro's number: items per mole
    const double R = 0.008314; // ideal gas constant, kJ/(K mol)
    double V = num_of_molecules_in_simulation/(Na*total_monomer_concentration); // Volume in L. For 1 mol/L the volume is L (which is also dm^3)
    double over_x; // overall conversion with respect to A functional groups!
    chain_pool all_chains; // vector of chain objects to store sequences and chain information in
    chain_pool loops; // vector of chain objects which have formed loops
    // Need to keep track of all unreacted functional groups, which are found in either monomers or chain ends
    std::vector<int> chainsA(monomermassA.size(),0); // number of chain ends of each A monomer type, initialize to 0. 
    // chainsA[0] stores # of A1 chain ends, [1] stores # of A2 chain ends, etc. 
    std::vector<int> chainsB(monomermassB.size(),0); // number of chain ends of each B monomer type, initialize to 0. 
    // chainsB[0] stores # of B1 chain ends, [1] stores # of A2 chain ends, etc.

    // For bifunctional monomers, the first dyad in a chain will always have frontend = 0 and backend = 1.
    std::ofstream molwt;
    time_t t = time(0);   // get current time
    struct tm * now = localtime( & t );
    char date[80];
    strftime (date,80,"%Y%m%d_%I%M%S%p_%Z",now); 
    std::string str(date);
    std::string mweight = "molecular_weight";
    std::string path = "/Users/mimivirus/Documents/Research/Polymer_Recycling_Project/MiscellaneousResources/KMC_PU/Output/";
    std::string molwtfilename = path+date+mweight+filename;
    molwt.open(molwtfilename); // store time, Mn, Mw, Dispersity
    molwt << "time           Conversion     Mn         Mw         D\n";

    // Reproduce Tirrell equation 67
    std::ofstream Tirrell;
    std::string Tirl = "Tirrell";
    std::string Tirrellfilename = path+date+Tirl+filename;
    Tirrell.open(Tirrellfilename); // store q1 and average sequence length
    Tirrell << "time        q1         Nnbb       q2         Nncc\n";
    double q1=0;
    double q2=0;

    // declare variables needed to track molecular weight
    double sumMiNi=0;
    double sumNi=0;
    double sumMi2Ni=0;
    double dispersity;
    double Mn=0;
    double Mw=0; 
    bool isloop = false;
    bool isnewchain=false;
    bool ismonomerA=false;
    bool ismonomerB=false;
    double Mi_A; double Mi_B;
    //initialize_molecular_weight(Mn, Mw, monomerA, monomerB, monomermassA, monomermassB, sumNi, sumMiNi, sumMi2Ni);
    
    // calculate all rate constants from kinetic parameters and temperature
    std::vector<double> kf(M);
    for (int i=0;i<M;i++) {
        kf[i]=A_kf[i]*exp(-Ea_kf[i]/(R*T)); // adding a to b
    }

    std::srand(std::time(0));
    double time = 0; // seconds

    // declare timing variables
    auto molwt_time=0;
    auto seq_update_time=0;

    // KMC loop
    while (time<simulation_time) {

        // calculate propensity functions
        double c[M];
        for (int i=0;i<M;i++) {
            c[i] = kf[i]/(V*Na); // 1 / (molecules s)
        }
        // update total rate
        double total_rate=0;
        // note: first entries in total rate are monomerA.size()*monomerB.size() islooping numbering reactions as R1 = A1B1, R2 = A2B1, R3 = A3B1... Rn=A2B1,Rn+1=A2B2,Rn+2=A2B3...
        int v=0;
        signed long long h[M];
        double Rv[M];
        for (int j=0;j<monomerA.size();j++) {
            for (int k=0;k<monomerB.size();k++) {
                //h[v]=(2*monomerA[j]+chainsA[j])*(2*monomerB[k]+chainsB[k]); // hv is the product of the numbers of molecular reactants involved in the vth reaction channel present at time t
                Rv[v]=c[v]*(2*monomerA[j]+chainsA[j])*(2*monomerB[k]+chainsB[k]); // molecules/s 2*molecules because each monomer is bifunctional
                total_rate+=Rv[v]; // molecules/s 2*molecules because each monomer is bifunctional
                v+=1;
            }
        }
        // choose timestep tau
        double r1 = 1.0*std::rand()/RAND_MAX; 
        while (r1==0){
            // this whole loop is to ensure that r1 is not 0 (on my mac, less than 1/2,000,000,000 chance per random number generation event)
            r1 = 1.0*std::rand()/RAND_MAX;
        }
        double tau = (1/total_rate)*log(1/r1); // calculate timestep tau
        
        // choose reaction to take place using Gillespie algorithm
        double r2 = 1.0*std::rand()/RAND_MAX;
        int mu=0;
        double sumRv=0; // Lin Wang eqn 1 multiplied by total rate
        while (sumRv<r2*total_rate) { // this is probably better implemented in a do isloop
            sumRv+=Rv[mu];
            if (sumRv<r2*total_rate) mu++;
        }
        // Translation from reaction channel mu to 
        // specific AA monomer type and BB monomer type 
        int A_monomer_type=0; int B_monomer_type=0; int count=0;
        specificmonomertype(A_monomer_type, B_monomer_type, monomerA, monomerB, mu);
        
        // pick which A, B functional group
        double r3 = 1.0*std::rand()/RAND_MAX;
        while (r3==1) {
            // this whole loop is to ensure that r3 excludes 1, want whichA on range [0, # of functional groups) (on my mac, less than 1/2,000,000,000 chance per random number generation event)
            r3 = 1.0*std::rand()/RAND_MAX;
        }
        double whichA = r3*(chainsA[A_monomer_type]+2*monomerA[A_monomer_type]);
        double r4 = 1.0*std::rand()/RAND_MAX;
        while (r4==1){
            // this whole loop is to ensure that whichB excludes the upper bound, want number on range [0, # of functional groups) (on my mac, less than 1/2,000,000,000 chance per random number generation event)
            r4 = 1.0*std::rand()/RAND_MAX;
        }
        double whichB = r4*(chainsB[B_monomer_type]+2*monomerB[B_monomer_type]); 
        // update chains
        auto start_seq_update = high_resolution_clock::now(); // measure molecular weight calculation time
        explicit_sequence_record(whichA,whichB,monomerA,monomerB,A_monomer_type,B_monomer_type,chainsA,chainsB,all_chains,loops,isloop,isnewchain,ismonomerA,ismonomerB,Mi_A,Mi_B,monomermassA,monomermassB);
        auto stop_seq_update = high_resolution_clock::now();
        auto duration_seq_update = duration_cast<microseconds>(stop_seq_update-start_seq_update);
        seq_update_time += duration_seq_update.count();
        time += tau;
        
        /* END KMC CALCULATIONS. 
           START RECORD SEQUENCE
        */ 
        // Record Mn, Mw, and dispersity
        auto start_molwt = high_resolution_clock::now(); // time molecular weight calculation
        molecular_weight(Mn, Mw, all_chains, loops, monomerA, monomerB, monomermassA, monomermassB,isloop,isnewchain,ismonomerA,ismonomerB,sumNi,sumMiNi,sumMi2Ni,Mi_A,Mi_B);
        dispersity=Mw/Mn; // calculate polydispersity index PDI
        // the below overall conversion line could be improved
        over_x=1-(chainsA[0]+2*monomerA[0]+chainsA[1]+2*monomerA[1])/(2*total_A_concentration); // calculate overall conversion of A functional group

        molwt << std::left << std::setw(10) << time << "     " << std::setw(10) << over_x << "     " << std::setw(6) << Mn << "     " << std::setw(6) << Mw << "     " << std::setw(6) << dispersity << "     " << "\n";
        auto stop_molwt = high_resolution_clock::now(); // stop time molecular weight calculation 
        auto duration_molwt = duration_cast<microseconds>(stop_molwt-start_molwt);
        molwt_time += duration_molwt.count();

        // do Tirrell stuff here
        // what is q1 at any given time? Extent of reaction of 

        // calculate sequence lengths inside all chains in all_chains. These are the sequence lengths defined by Tirrell rather than
        // total chain lengths that correspond to number average degree of polymerization

        int chain_number = 0;
        int monomer_number=0;
        int last_amine=-1;
        int sequence_length=0;
        int num_seqs_A0=0;
        int num_seqs_A1=0;
        double NnA0=0;
        double NnA1=0;
        int total_sequence_length_A0=0;
        int total_sequence_length_A1=0;
        // iterate through each chain in the chain pool
        while (chain_number < all_chains.size()) {
            // iterate through each element in the selected chain
            while (monomer_number<all_chains[chain_number].v.size()) {
                // if sequence_length is 0, initialize new sequence
                if (sequence_length==0) {
                    if (all_chains[chain_number].v[monomer_number][0]==0) {
                        // what is the identity of the first amine in the sequence? i.e. for AA+BB+CC reaction where BB and CC do not react with each other, is the first element in a chain BB or CC?
                        // either A0 or A1
                        last_amine=all_chains[chain_number].v[monomer_number][1]; 
                    }
                    else {
                        // if the first element is not an amine, then the 2nd element must be amine
                        // i.e. if the first element is B, then the second one is A
                        if ((monomer_number+1)<all_chains[chain_number].v.size()) {
                            last_amine=all_chains[chain_number].v[monomer_number+1][1];
                            monomer_number+=1;
                        }
                    }
                }
                if (all_chains[chain_number].v[monomer_number][0]!=0 || last_amine==-1) {
                    // something is wrong
                    printf("Something went wrong");
                }
                // if the next amine is the same as the last amine, increment the sequence length by 1
                if (all_chains[chain_number].v[monomer_number][1]==last_amine){
                    sequence_length+=1;
                }
                // if the next amine is not the same as the last amine, or if the end of the chain is reached, END THE SEQUENCE
                // i.e. update the total sequence length and the num of sequences
                if (all_chains[chain_number].v[monomer_number][1]!=last_amine || (monomer_number+2)>=all_chains[chain_number].v.size()) {
                    if (last_amine==0) {
                        total_sequence_length_A0+=sequence_length;
                        num_seqs_A0+=1;
                    }
                    else if (last_amine==1) {
                        total_sequence_length_A1+=sequence_length;
                        num_seqs_A1+=1;
                    }
                    else {
                        printf("something went wrong");
                    }
                    // start a new sequence
                    sequence_length=0;
                }
                monomer_number+=2;
            }
            monomer_number=0;
            chain_number+=1;
        }
        if (num_seqs_A0 != 0) {
            NnA0=(1.0)*(total_sequence_length_A0+monomerA[0])/(num_seqs_A0+monomerA[0]);
        }
        if (num_seqs_A1 != 0) {
            NnA1=((1.0)*total_sequence_length_A1+monomerA[1])/(num_seqs_A1+monomerA[1]);
        }

        // calculate q1 and q2
        q1 = (2.0*startingA0-(2.0*monomerA[0]+chainsA[0]))/(2.0*startingA0);
        q2 = (2.0*startingA1-(2.0*monomerA[1]+chainsA[1]))/(2.0*startingA1);
        //    Tirrell << "time       q1         Nnbb       q2         Nncc\n";
        Tirrell << std::left << std::setw(7) << time << "     " << std::setw(6) << q1 << "     " << std::setw(6) << NnA0 << "     " << std::setw(6) << q2 << "     " << std::setw(6) << NnA1 << "     \n";
        }

    molwt.close();
    Tirrell.close();

    // could easily print all sequences here too if desired
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop-start);
    auto totalprogramtime = duration.count();
    // total number of events that can occur is 1 for each reaction 
    return 0;
}