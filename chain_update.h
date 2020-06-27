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
    std::vector<std::vector<int>> v;
    int frontend; // identity of species at start of vector
    int backend; // identity of species at end of vector
};
// Need to keep track of all the polymer chains that have been created
typedef std::vector<chain> chain_pool;

void explicit_sequence_record(int whichA, int whichB, std::array<int,2>& moleculesA, std::array<int,1>& moleculesB, int& A_index, int& B_index, std::array<int, 2>& chainsA, std::array<int, 1>& chainsB, chain_pool& all_chains, chain_pool& loops) {
    bool front = false;
    bool back = false;
    bool something_went_wrong = false;
    // case for which 1 A monomer reacts with 1 B monomer
    if (whichA<2*moleculesA[A_index] && whichB<2*moleculesB[B_index]) {
        chain newchain; // make  a new chain
        newchain.frontend = 0; // i.e. front end is A
        newchain.backend = 1; // backend is defined arbitrarily as species B on new polymer chain
        newchain.v.push_back({newchain.frontend,A_index});
        newchain.v.push_back({newchain.backend,B_index});
        all_chains.push_back(newchain);
        moleculesA[A_index]--; // update the trackers
        moleculesB[B_index]--; // update the trackers
        chainsA[A_index]++; // update the trackers
        chainsB[B_index]++;  // update the trackers
    }
    // case for which 1 A monomer reacts with 1 B chain
    else if (whichA<2*moleculesA[A_index] && whichB>2*moleculesB[B_index]) {
        // STEP 1: Select B chain
        int Bselect = (int) (whichB-2*moleculesB[B_index]);
        int selected_chain = 0;
        int countB = 0;

        while (selected_chain<all_chains.size()){
            // select which B chain and which B end
            if (all_chains[selected_chain].backend==1 && all_chains[selected_chain].v.back()[1]==B_index) { // if B functional group on back AND correct type of B monomer
                countB++;
            }
            if (Bselect<countB) {
                back=true;
                break;
            }
            if (all_chains[selected_chain].frontend==1 && all_chains[selected_chain].v.front()[1]==B_index){ // if B functional group on front AND correct type of B monomer
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
        if (back) { // these two conditions are redundant! Just use if (back) {} 1 == all_chains[selected_chain].v.back()[0] && 
            // then add A to the end
            all_chains[selected_chain].v.push_back({0,A_index}); // can I store a pointer to all_chains[selected_chain] as something shorter?
            // and update the trackers accordingly
            all_chains[selected_chain].backend=0;
            chainsB[B_index]--;
            chainsA[A_index]++;
            moleculesA[A_index]--;
        }
        // CASE 2: if the B is at the front of the chain
        else if (front) { //1 == all_chains[selected_chain].v.front()[0] && 
            // then add A to the front
            all_chains[selected_chain].v.insert(all_chains[selected_chain].v.begin(),{0,A_index});
            // and update the trackers accordingly
            all_chains[selected_chain].frontend=0;
            chainsB[B_index]--;
            chainsA[A_index]++;
            moleculesA[A_index]--;
        }
        // could add a line to check that this is occurring properly
    }
    // case for which 1 B monomer reacts with 1 A chain
    else if (whichA>2*moleculesA[A_index] && whichB<2*moleculesB[B_index]) {
        // STEP 1: Select A chain
        int Aselect = (int) (whichA-2*moleculesA[A_index]);
        int selected_chain = 0;
        int countA = 0;
        while (selected_chain<all_chains.size()) {
            // select which A chain and which A end
            if (all_chains[selected_chain].backend==0 && all_chains[selected_chain].v.back()[1]==A_index) { // if A functional group on back AND correct type of A monomer
                countA++;
            }
            if (Aselect<countA) {
                back=true;
                break;
            }
            if (all_chains[selected_chain].frontend==0 && all_chains[selected_chain].v.front()[1]==A_index) { // if A functional group on front AND correct type of A monomer
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
        if (back) { //0 == all_chains[selected_chain].v.back()[0] && 
            // then add B to the end
            all_chains[selected_chain].v.push_back({1,B_index}); // can I store a pointer to all_chains[selected_chain] as something shorter?
            // and update the trackers accordingly
            all_chains[selected_chain].backend=1;
            chainsB[B_index]++;
            chainsA[A_index]--;
            moleculesB[B_index]--;
        }
        // CASE 2: if the A is at the front of the chain
        else if (front) { //0 == all_chains[selected_chain].v.front()[0] && 
            // then add B to the front
            all_chains[selected_chain].v.insert(all_chains[selected_chain].v.begin(),{1,B_index});
            // and update the trackers accordingly
            all_chains[selected_chain].frontend=1;
            chainsB[B_index]++;
            chainsA[A_index]--;
            moleculesB[B_index]--;
        }
        // could add a line to check that this is occuring properly
    }
    // case for which 1 A chain reacts with 1 B chain
    else if (whichA>2*moleculesA[A_index] && whichB>2*moleculesB[B_index]) {
        // STEP 1: select A chain
        int Aselect = (int) (whichA-2*moleculesA[A_index]);
        int selected_A_chain = 0;
        int countA = 0;
        bool frontA=false;
        bool backA=false;

        while (selected_A_chain<all_chains.size()) {
            // select which A chain and which A end
            if (all_chains[selected_A_chain].backend==0 && all_chains[selected_A_chain].v.back()[1]==A_index){ // if A functional group on back AND correct type of A monomer
                countA++;
            }
            if (Aselect<countA) {
                backA=true;
                break;
            }
            if (all_chains[selected_A_chain].frontend==0 && all_chains[selected_A_chain].v.front()[1]==A_index){ // if A functional group on back AND correct type of A monomer
                countA++;
            }
            if (Aselect<countA) {
                frontA=true;
                break;
            }
            selected_A_chain++;
        }
        // STEP 2: select B chain
        int Bselect = (int) (whichB-2*moleculesB[B_index]);
        int selected_B_chain = 0;
        int countB = 0;
        bool backB=false;
        bool frontB=false;
        while (selected_B_chain<all_chains.size()){
            // select which B chain and which B end
            if (all_chains[selected_B_chain].backend==1 && all_chains[selected_B_chain].v.back()[1]==B_index){ // if B functional group on back AND correct type of B monomer
                countB++;
            }
            if (Bselect<countB) {
                backB=true;
                break;
            }
            if (all_chains[selected_B_chain].frontend==1 && all_chains[selected_B_chain].v.front()[1]==B_index){ // if B functional group on front AND correct type of B monomer
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
                all_chains.erase(all_chains.begin()+selected_B_chain);
            }
            else {
                // add A chain to back of B chain, delete A chain
                all_chains[selected_B_chain].v.insert(all_chains[selected_B_chain].v.end(),all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_B_chain].backend=all_chains[selected_A_chain].backend;
                all_chains.erase(all_chains.begin()+selected_A_chain); 
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
                std::iter_swap(all_chains.begin()+selected_A_chain,all_chains.begin()+selected_B_chain); 
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
        chainsA[A_index]--;
        chainsB[B_index]--;
    }
    
}
