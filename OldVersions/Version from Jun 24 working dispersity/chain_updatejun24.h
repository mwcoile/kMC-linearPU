#pragma once
#include<iomanip>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
struct chain
{
    std::vector<int> v;
    int frontend; // identity of species at start of vector
    int backend; // identity of species at end of vector
};
// Need to keep track of all the polymer chains that have been created
typedef std::vector<chain> chain_pool;

void explicit_sequence_record(int whichA, int whichB, int& moleculesA, int& moleculesB, int& chainsA, int& chainsB, chain_pool& all_chains, chain_pool& loops) {
    bool front = false;
    bool back = false;
    bool something_went_wrong = false;
    // case for which 1 A monomer reacts with 1 B monomer
    if (whichA<2*moleculesA && whichB<2*moleculesB) {
        chain newchain;
        newchain.frontend = 0; // i.e. front end is A
        newchain.backend = 1; // backend is defined arbitrarily as species B on new polymer chain
        newchain.v.push_back(newchain.frontend);
        newchain.v.push_back(newchain.backend);
        all_chains.push_back(newchain);
        moleculesA--;
        moleculesB--;
        chainsA++;
        chainsB++;  
    }
    // case for which 1 A monomer reacts with 1 B chain
    else if (whichA<2*moleculesA && whichB>2*moleculesB) {
        // STEP 1: Select B chain
        int Bselect = (int) (whichB-2*moleculesB);
        int selected_chain = 0;
        int countB = 0;

        while (selected_chain<all_chains.size()){
            // select which B chain and which B end
            if (all_chains[selected_chain].backend==1){
                countB++;
            }
            if (Bselect<countB) {
                back=true;
                break;
            }
            if (all_chains[selected_chain].frontend==1){
                countB++;
            }
            if (Bselect<countB) {
                front=true;
                break;
            }
            selected_chain++;
        }

        // STEP 2: which end of the chain has B?
        // CASE 1: if the B is at the end of the chain
        if (1 == all_chains[selected_chain].v.back() && back) {
            // then add A to the end
            all_chains[selected_chain].v.push_back(0); // can I store a pointer to all_chains[selected_chain] as something shorter?
            // and update the trackers accordingly
            all_chains[selected_chain].backend=0;
            chainsB--;
            chainsA++;
            moleculesA--;
        }
        // CASE 2: if the B is at the front of the chain
        else if (1 == all_chains[selected_chain].v.front() && front) {
            // then add A to the front
            all_chains[selected_chain].v.insert(all_chains[selected_chain].v.begin(),0);
            // and update the trackers accordingly
            all_chains[selected_chain].frontend=0;
            chainsB--;
            chainsA++;
            moleculesA--;
        }
        // could add a line to check that this is occuring properly
    }
    // case for which 1 B monomer reacts with 1 A chain
    else if (whichA>2*moleculesA && whichB<2*moleculesB) {
        // STEP 1: Select A chain
        int Aselect = (int) (whichA-2*moleculesA);
        int selected_chain = 0;
        int countA = 0;
        while (selected_chain<all_chains.size()) {
            // select which A chain and which A end
            if (all_chains[selected_chain].backend==0){
                countA++;
            }
            if (Aselect<countA) {
                back=true;
                break;
            }
            if (all_chains[selected_chain].frontend==0){
                countA++;
            }
            if (Aselect<countA) {
                front=true;
                break;
            }
            selected_chain++;
        }

        // STEP 2: which end of the chain has A?
        // CASE 1: if the A is at the end of the chain .. but what to do if both ends have A?
        if (0 == all_chains[selected_chain].v.back() && back) {
            // then add B to the end
            all_chains[selected_chain].v.push_back(1); // can I store a pointer to all_chains[selected_chain] as something shorter?
            // and update the trackers accordingly
            all_chains[selected_chain].backend=1;
            chainsB++;
            chainsA--;
            moleculesB--;
        }
        // CASE 2: if the A is at the front of the chain
        else if (0 == all_chains[selected_chain].v.front() && front) {
            // then add B to the front
            all_chains[selected_chain].v.insert(all_chains[selected_chain].v.begin(),1);
            // and update the trackers accordingly
            all_chains[selected_chain].frontend=1;
            chainsB++;
            chainsA--;
            moleculesB--;
        }
        // could add a line to check that this is occuring properly
    }
    // case for which 1 A chain reacts with 1 B chain
    else if (whichA>2*moleculesA && whichB>2*moleculesB) {
        // STEP 1: select A chain
        int Aselect = (int) (whichA-2*moleculesA);
        int selected_A_chain = 0;
        int countA = 0;
        bool frontA=false;
        bool backA=false;

        while (selected_A_chain<all_chains.size()) {
            // select which A chain and which A end
            if (all_chains[selected_A_chain].backend==0){
                countA++;
            }
            if (Aselect<countA) {
                backA=true;
                break;
            }
            if (all_chains[selected_A_chain].frontend==0){
                countA++;
            }
            if (Aselect<countA) {
                frontA=true;
                break;
            }
            selected_A_chain++;
        }
        // STEP 2: select B chain
        int Bselect = (int) (whichB-2*moleculesB);
        int selected_B_chain = 0;
        int countB = 0;
        bool backB=false;
        bool frontB=false;
        while (selected_B_chain<all_chains.size()){
            // select which B chain and which B end
            if (all_chains[selected_B_chain].backend==1){
                countB++;
            }
            if (Bselect<countB) {
                backB=true;
                break;
            }
            if (all_chains[selected_B_chain].frontend==1){
                countB++;
            }
            if (Bselect<countB) {
                frontB=true;
                break;
            }
            selected_B_chain++;
        }
        // which chain is getting added to and which chain is getting deleted?
        bool add_to_chain_A = false;
        if (selected_A_chain<selected_B_chain) {
            add_to_chain_A = true;
        }
        // case 0: loop formation
        if (selected_A_chain==selected_B_chain){
            // delete from vector of chains
            // add to loops vector. This does not currently consider whether sterically it is possible for this to occur (i.e. loops consisting of 2 monomers are permitted)
            loops.push_back(all_chains[selected_A_chain]);
            all_chains.erase(all_chains.begin()+selected_A_chain); // what the heck is an iterator, and what's the difference between it and a const_iterator
        }
        // case 1: front of A chain to front of B chain
        else if (frontA && frontB) {
            if (add_to_chain_A) {
                // reverse A chain, append B chain to end of A chain, update the ends, delete B chain
                std::reverse(all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_A_chain].v.insert(all_chains[selected_A_chain].v.end(),all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end());
                all_chains[selected_A_chain].frontend=all_chains[selected_A_chain].backend;
                all_chains[selected_A_chain].backend=all_chains[selected_B_chain].backend;
                all_chains.erase(all_chains.begin()+selected_B_chain);
            }
            else {
                // reverse B chain, append A to end of B chain, update the ends, delete A chain
                std::reverse(all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end()); 
                all_chains[selected_B_chain].v.insert(all_chains[selected_B_chain].v.end(),all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_B_chain].frontend=all_chains[selected_B_chain].backend;
                all_chains[selected_B_chain].backend=all_chains[selected_A_chain].backend;
                all_chains.erase(all_chains.begin()+selected_A_chain);
            }
        }
        // case 2: front of A chain to back of B chain
        else if (frontA && backB) {
            if (add_to_chain_A) {
                // add A chain to B chain, then put it in the A chain location, then delete chain B
                all_chains[selected_B_chain].v.insert(all_chains[selected_B_chain].v.end(),all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_B_chain].backend=all_chains[selected_A_chain].backend;
                std::iter_swap(all_chains.begin()+selected_A_chain,all_chains.begin()+selected_B_chain); 
                // I'm not sure that the above and below lines are doing what I am intending
                all_chains.erase(all_chains.begin()+selected_B_chain);
            }
            else {
                // add A chain to back of B chain, delete A chain
                all_chains[selected_B_chain].v.insert(all_chains[selected_B_chain].v.end(),all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_B_chain].backend=all_chains[selected_A_chain].backend;
                all_chains.erase(all_chains.begin()+selected_A_chain); // modified this line 2:44pm on Jun 15
            }
        }
        // case 3: back of A chain to front of B chain
        else if (backA && frontB) {
            if (add_to_chain_A) {
                // add B chain to back of A chain
                all_chains[selected_A_chain].v.insert(all_chains[selected_A_chain].v.end(),all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end());
                all_chains[selected_A_chain].backend=all_chains[selected_B_chain].backend;
                all_chains.erase(all_chains.begin()+selected_B_chain);
            }
            else {
                // add B chain to back of A chain, then put it in the B chain location, then
                all_chains[selected_A_chain].v.insert(all_chains[selected_A_chain].v.end(),all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end());
                all_chains[selected_A_chain].backend=all_chains[selected_B_chain].backend;
                //all_chains.at(selected_B_chain)=all_chains[selected_A_chain]; 
                std::iter_swap(all_chains.begin()+selected_A_chain,all_chains.begin()+selected_B_chain); 
                // I'm not sure that the above and below lines are doing what I am intending
                all_chains.erase(all_chains.begin()+selected_A_chain);
            }
        }
        // case 4: back of A chain to back of B chain
        else if (backA && backB) {
            if (add_to_chain_A) {
                // reverse chain B, add to end of chain A, delete chain B
                std::reverse(all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end());
                all_chains[selected_A_chain].v.insert(all_chains[selected_A_chain].v.end(),all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end());
                all_chains[selected_A_chain].backend=all_chains[selected_B_chain].frontend;
                all_chains.erase(all_chains.begin()+selected_B_chain);
            }
            else {
                // reverse chain A, add to end of chain B, delete chain A
                std::reverse(all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_B_chain].v.insert(all_chains[selected_B_chain].v.end(),all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_B_chain].backend=all_chains[selected_A_chain].frontend;
                all_chains.erase(all_chains.begin()+selected_A_chain);
            }
        }
        else
        {
            something_went_wrong = true;
        }
        chainsA--;
        chainsB--;
    }
    
}
