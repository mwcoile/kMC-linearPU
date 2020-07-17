#include<iomanip>
#include<stdlib.h>
#include<iostream>
#include<string>
#include <vector>
#include<math.h>
#include<array>
#include<chrono>
// author: Matthew
// playground file for understanding C++
//template class std::vector<MyClass>;

using namespace std::chrono;

int main() {
    auto start = high_resolution_clock::now();
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
    int tester[5] = {0};
    const int bignumber=28;
    double testmore[bignumber]={};

    std::vector<std::vector<int>> foo123={};
    foo123.push_back({1,80});
    int testermatthew=foo123[0][1];
    std::array<int,3> helloarrays = {1,2,8};
    int number = helloarrays.size();

    unsigned long long int a[65],i;
    a[0]=1;
    for(i=1;i<65;i++) {
        a[i]=2<<(i-1);
        printf("i=%d a[i]=%lld\n",i, a[i]);
    }
    int q = 2<<12;
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );

    char buffer [80];
    strftime (buffer,80,"%Y%m%d_%I%M%S%p_%Z",now);
    std::string phooey = "lolz";
    phooey.append(buffer);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop-start);
    
    auto sixtytwo = duration.count();
    int amistupid=RAND_MAX;
    double tauy = log(1/0);
    return 0;

}