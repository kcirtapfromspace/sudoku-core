// Randomized exact-cover completion and independent givens uniqueness discovery.
// Input: mode seed limit node_limit; then81 masks(mode1) or puzzle(mode0).
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <numeric>
struct Cover {
    std::array<int, 4000> l{},r{},u{},d{},c{},row{};
    std::array<int, 325> size{};
    int next=325, limit=8;
    uint64_t nodes=0, maxnodes=5000000;
    bool cutoff=false;
    std::vector<int> selected;
    std::vector<std::string> solutions;
    Cover(const std::array<int,81>& domains, uint64_t seed) {
        for(int i=0;i<=324;i++) {l[i]=i-1;r[i]=i+1;u[i]=d[i]=i;}
        l[0]=324;r[324]=0;
        std::vector<int> choices(729);std::iota(choices.begin(),choices.end(),0);
        if(seed) {std::mt19937_64 rng(seed);std::shuffle(choices.begin(),choices.end(),rng);}
        for(int choice:choices) {
            int cell=choice/9,digit=choice%9;
            int rr=cell/9,cc=cell%9,b=rr/3*3+cc/3;
            {
                if(!(domains[cell] & (1<<digit))) continue;
                int cols[]={cell+1,82+rr*9+digit,163+cc*9+digit,244+b*9+digit};
                int first=next;
                for(int col:cols) {
                    int n=next++; c[n]=col; row[n]=cell*9+digit;
                    u[n]=u[col];d[n]=col;d[u[col]]=n;u[col]=n;size[col]++;
                    l[n]=n-1;r[n]=n+1;
                }
                l[first]=next-1;r[next-1]=first;
            }
        }
    }
    void cover(int col) {
        r[l[col]]=r[col];l[r[col]]=l[col];
        for(int i=d[col];i!=col;i=d[i])
            for(int j=r[i];j!=i;j=r[j]) {d[u[j]]=d[j];u[d[j]]=u[j];size[c[j]]--;}
    }
    void uncover(int col) {
        for(int i=u[col];i!=col;i=u[i])
            for(int j=l[i];j!=i;j=l[j]) {size[c[j]]++;d[u[j]]=j;u[d[j]]=j;}
        r[l[col]]=col;l[r[col]]=col;
    }
    void search() {
        if(solutions.size()>=static_cast<unsigned>(limit) || cutoff) return;
        if(++nodes>maxnodes) {cutoff=true;return;}
        if(r[0]==0) {
            std::string solution(81,'.');
            for(int id:selected) solution[id/9]='1'+id%9;
            solutions.push_back(solution);return;
        }
        int col=r[0];
        for(int j=r[col];j!=0;j=r[j]) if(size[j]<size[col]) col=j;
        if(size[col]==0) return;
        cover(col);
        for(int i=d[col];i!=col;i=d[i]) {
            selected.push_back(row[i]);
            for(int j=r[i];j!=i;j=r[j]) cover(c[j]);
            search();
            for(int j=l[i];j!=i;j=l[j]) uncover(c[j]);
            selected.pop_back();
            if(solutions.size()>=static_cast<unsigned>(limit) || cutoff) break;
        }
        uncover(col);
    }
};
int main() {
    int mode,limit;uint64_t seed,maxnodes;
    while(std::cin>>mode>>seed>>limit>>maxnodes) {
        std::array<int,81> domains{};
        if(mode==0) {
            std::string puzzle;std::cin>>puzzle;if(puzzle.size()!=81)return 2;
            for(int c=0;c<81;c++)domains[c]=(puzzle[c]>='1'&&puzzle[c]<='9')?(1<<(puzzle[c]-'1')):511;
        } else {
            for(int& mask:domains)if(!(std::cin>>mask)||mask<1||mask>511)return 2;
        }
        if(limit<1)return 2;
        Cover solver(domains,seed);solver.limit=limit;solver.maxnodes=maxnodes;solver.search();
        std::cout<<"{\"nodes\":"<<solver.nodes<<",\"cutoff\":"<<(solver.cutoff?"true":"false")<<",\"solutions\":[";
        for(unsigned i=0;i<solver.solutions.size();i++) {
            if(i)std::cout<<',';std::cout<<'"'<<solver.solutions[i]<<'"';
        }
        std::cout<<"]}"<<std::endl;
    }
}
