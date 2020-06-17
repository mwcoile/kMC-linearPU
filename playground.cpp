#include<iomanip>
#include<stdlib.h>
#include<iostream>
#include<string>
#include <vector>
#include<math.h>

// author: Matthew
// playground file for understanding C++
//template class std::vector<MyClass>;

int main() {
    std::srand(std::time(NULL));
    double r1 = 1.0*std::rand()/RAND_MAX; 
    double r2 = 1.0*std::rand()/RAND_MAX;
    int foo = RAND_MAX;
    std::vector<int> test123 {1,2,3,4,5,6};

    // At some point, make sure I know where to find the output from printf... debug folder?
    printf("the 2nd value in the array test123 is %i", test123[1]);
    //double x=test123[1];
    std::reverse(test123.begin(),test123.end());
    std::iter_swap(test123.begin()+1,test123.begin()+4);
    test123.erase(test123.begin()+4);
    struct chain
    {
        std::vector<int> v;
        int frontend; // identity of species at start of vector
        int backend; // identity of species at end of vector
    };
    typedef std::vector<chain> chain_pool;
    chain_pool all_chains;
    chain chain1;
    chain1.frontend=1;
    chain1.backend=0;
    chain1.v = {0,0,1,6};
    all_chains.push_back(chain1);
    all_chains.push_back(chain1);
    all_chains.push_back(chain1);
    all_chains.push_back(chain1);
    all_chains.erase(all_chains.begin()+1);

    bool bar = (all_chains[0].v[0]==1);
    bool barend = (all_chains[0].v[all_chains[0].v.back()]==0);
    bool hmmm = (all_chains[0].v.back()==3);
    int sizer = all_chains[0].v.size(); // this is 4
    int lol = all_chains[0].v.front(); // I understand why this is 0
    int wtf = all_chains[0].v.back(); // but why is this also 0?
    //int testagain=all_chains[0].v.begin();

    double Na=6.02;
    int num_of_molecules_in_simulation=1E4;
    double V = num_of_molecules_in_simulation/Na;
    return 0;

}