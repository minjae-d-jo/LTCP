#include <fstream>
#include <sstream>
#include <iostream>
#include <RandomNumber.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace std;
using namespace Snu::Cnrc;
using Size = unsigned int;
using Time = unsigned int;

void branching(vector<Size>& upSpinSite, vector<vector<bool>>& state, const Size L, Size X, Size Y) {
    if(state[X][Y]==0) {
        state[X][Y] = 1;
        upSpinSite.push_back(L*X+Y);
    }
}

void LongRangeBranching(vector<Size>& upSpinSite, vector<vector<bool>>& state, const Size L, Size X, Size Y, const double sigma) {
    double x, y, M;
    long long dX, dY;
    long long r;
    Size XLR, YLR;
    RandomRealGenerator rnd(0.0, 1.0);
    RandomRealGenerator rnd2(-1.0, 1.0);
    r = pow(rnd(),-1./sigma);
    if(r<sqrt(2)) { // Nearest neghbor interaction
        double toss_1=rnd();
        Size XLeft = X == 0 ? L-1 : X-1;
        Size XRight = X == L-1 ? 0 : X+1;
        Size YBottom = Y == 0 ? L-1 : Y-1;
        Size YTop = Y == L-1 ? 0 : Y+1;
        if( toss_1<0.25 ) branching(upSpinSite, state, L, XLeft, Y);
        else if( toss_1<0.5 ) branching(upSpinSite, state, L, XRight, Y);
        else if( toss_1<0.75 ) branching(upSpinSite, state, L, X, YTop);
        else branching(upSpinSite, state, L, X, YBottom);
    }
    else {
        do {
            x = rnd2(), y = rnd2(); // (1)
        } while(x*x+y*y>1);
        M = sqrt(x*x+y*y);
        dX = round(x*r/M), dY = round(y*r/M); // (2)
        if(dX>0) XLR = X+dX > L-1 ? (X+dX)%L : X+dX;
        else XLR = X+dX < 0 ? L-((-(dX+X))%L) : X+dX;
        if(XLR==L) XLR=0;
        if(dY>0) YLR = Y+dY > L-1 ? (Y+dY)%L : Y+dY;
        else YLR = Y+dY < 0 ? L-((-(dY+Y))%L) : Y+dY;
        if(YLR==L) YLR=0;
        if(state[XLR][YLR]==0) {
            state[XLR][YLR] = 1;
            upSpinSite.push_back(L*XLR+YLR);
        }
    }
}

/*
void LongRangePairBranching(vector<Size>& upSpinSite, vector<vector<bool>>& state, const Size L, 
    Size XNeigh, Size YNeigh, Size X, Size Y) {
    RandomRealGenerator rnd(0.0, 1.0);
    double toss_1=rnd();
    if(XNeigh==X) {
        Size XLeft = X == 0 ? L-1 : X-1;
        Size XRight = X == L-1 ? 0 : X+1;
        if( toss_1<1./3 ) branching(upSpinSite, state, L, XLeft, Y);
        else if( toss_1<2./3 ) branching(upSpinSite, state, L, XRight, Y);
        else { 
            Size YBottom = Y == 0 ? L-1 : Y-1;
            Size YTop = Y == L-1 ? 0 : Y+1;
            if(YBottom!=YNeigh) branching(upSpinSite, state, L, X, YBottom);
            else branching(upSpinSite, state, L, X, YTop);
        }
    }
    else if(YNeigh==Y) {
        Size YBottom = Y == 0 ? L-1 : Y-1;
        Size YTop = Y == L-1 ? 0 : Y+1;
        if( toss_1<1./3 ) branching(upSpinSite, state, L, X, YTop);
        else if( toss_1<2./3 ) branching(upSpinSite, state, L, X, YBottom);
        else { 
            Size XLeft = X == 0 ? L-1 : X-1;
            Size XRight = X == L-1 ? 0 : X+1;
            if(XLeft!=XNeigh) branching(upSpinSite, state, L, XLeft, Y);
            else branching(upSpinSite, state, L, XRight, Y);
        }
    }
}
*/


void LongRangePairBranching(vector<Size>& upSpinSite, vector<vector<bool>>& state, const Size L, 
    Size XNeigh, Size YNeigh, Size X, Size Y, const double sigma) {
    double x, y, M;
    long long dX, dY;
    long long r;
    Size XLR, YLR;
    RandomRealGenerator rnd(0.0, 1.0);
    double toss_1=rnd();
    RandomRealGenerator rnd2(-1.0, 1.0);
    r = pow(rnd(),-1./sigma);
    if(r<sqrt(2)) { // Nearest neghbor interaction
        if(XNeigh-X ==0) {
            Size XLeft = X == 0 ? L-1 : X-1;
            Size XRight = X == L-1 ? 0 : X+1;
            if( toss_1<1./3 ) branching(upSpinSite, state, L, XLeft, Y);
            else if( toss_1<2./3 ) branching(upSpinSite, state, L, XRight, Y);
            else { 
                Size YBottom = Y == 0 ? L-1 : Y-1;
                Size YTop = Y == L-1 ? 0 : Y+1;
                if(YBottom!=YNeigh) branching(upSpinSite, state, L, X, YBottom);
                else branching(upSpinSite, state, L, X, YTop);
            }
        }
        else if(YNeigh-Y ==0) {
            Size YBottom = Y == 0 ? L-1 : Y-1;
            Size YTop = Y == L-1 ? 0 : Y+1;
            if( toss_1<1./3 ) branching(upSpinSite, state, L, X, YTop);
            else if( toss_1<2./3 ) branching(upSpinSite, state, L, X, YBottom);
            else { 
                Size XLeft = X == 0 ? L-1 : X-1;
                Size XRight = X == L-1 ? 0 : X+1;
                if(XLeft!=XNeigh) branching(upSpinSite, state, L, XLeft, Y);
                else branching(upSpinSite, state, L, XRight, Y);
            }
        }
    }
    else { // Long-range interaction
        do {
            x = rnd2(), y = rnd2(); // (1)
        } while(x*x+y*y>1);
        M = sqrt(x*x+y*y);
        dX = round(x*r/M), dY = round(y*r/M); // (2)
        if(dX>0) XLR = X+dX > L-1 ? (X+dX)%L : X+dX;
        else XLR = X+dX < 0 ? L-((-dX-X)%L) : X+dX;
        if(XLR==L) XLR=0;
        if(dY>0) YLR = Y+dY > L-1 ? (Y+dY)%L : Y+dY;
        else YLR = Y+dY < 0 ? L-((-dY-Y)%L) : Y+dY;
        if(YLR==L) YLR=0;
        
        if(state[XLR][YLR]==0) {
            state[XLR][YLR] = 1;
            upSpinSite.push_back(L*XLR+YLR);
        }
    }
}
void LongRangeCPUpdate(vector<Size>& upSpinSite, vector<vector<bool>>& state, const Size L, const double p, 
    const Size X, const Size Y, const Size order, const double sigma) {
    RandomRealGenerator rnd(0.0, 1.0);
    double toss_1=rnd(); //, toss_2=rnd(), zL, zR;
    
    if( toss_1<p ) {
        LongRangeBranching(upSpinSite, state, L, X, Y, sigma);
    }
    else {
        state[X][Y] = 0;
        upSpinSite.erase(upSpinSite.begin()+order); 
    }
}

/*
void CPUpdateGivenNeigh(vector<Size>& upSpinSite, vector<vector<bool>>& state, const Size L, const double p, 
    const Size X, const Size Y, const Size XNeigh, const Size YNeigh, const Size order, const double sigma) {
    RandomRealGenerator rnd(0.0, 1.0);
    double toss_1=rnd();

    if( toss_1<p ) {
        double x, y, M;
        long long dX, dY;
        long long r;
        Size XLR, YLR;
        RandomRealGenerator rnd(0.0, 1.0);
        RandomRealGenerator rnd2(-1.0, 1.0);
        r = pow(rnd(),-1./sigma);
        if(r<sqrt(2)) { // Nearest neghbor interaction
            state[XNeigh][YNeigh] = 1;
            upSpinSite.push_back(L*XNeigh+YNeigh);
        }
        else {
            do {
                x = rnd2(), y = rnd2(); // (1)
            } while(x*x+y*y>1);
            M = sqrt(x*x+y*y);
            dX = round(x*r/M), dY = round(y*r/M); // (2)
            if(dX>0) XLR = X+dX > L-1 ? (X+dX)%L : X+dX;
            else XLR = X+dX < 0 ? L-((-(dX+X))%L) : X+dX;
            if(XLR==L) XLR=0;
            if(dY>0) YLR = Y+dY > L-1 ? (Y+dY)%L : Y+dY;
            else YLR = Y+dY < 0 ? L-((-(dY+Y))%L) : Y+dY;
            if(YLR==L) YLR=0;
            if(state[XLR][YLR]==0) {
                state[XLR][YLR] = 1;
                upSpinSite.push_back(L*XLR+YLR);
            }
        }
    }
    else {
        state[X][Y] = 0;
        upSpinSite.erase(upSpinSite.begin()+order); 
    }
}
*/

void LongRangeTCPUpdate(vector<Size>& upSpinSite, vector<vector<bool>>& state, const Size L, const double p, 
    const double q, const double sigma) {
    RandomRealGenerator rnd(0.0, 1.0);
    RandomIntGenerator rndInt(0, upSpinSite.size()-1);
    Size order = rndInt();
    Size X = upSpinSite[order]/L;
    Size Y = upSpinSite[order]%L;
    double toss_1=rnd();
    
        
    if( toss_1<q ) { // pair creation
        double toss_p=rnd();
        Size XLeft = X == 0 ? L-1 : X-1;
        Size XRight = X == L-1 ? 0 : X+1;
        Size YBottom = Y == 0 ? L-1 : Y-1;
        Size YTop = Y == L-1 ? 0 : Y+1;
        if( toss_p<p ) {
            double toss_2=rnd(), toss_3=rnd();
            if( toss_2<0.25 ) { 
                if(state[XLeft][Y]==1) { // if left is up
                    if( toss_3<0.5 ) LongRangePairBranching(upSpinSite, state, L, X, Y, XLeft, Y, sigma);
                    else LongRangePairBranching(upSpinSite, state, L, XLeft, Y, X, Y, sigma);
                }
            }
            else if( toss_2<0.5 ) { 
                if(state[XRight][Y]==1) { // if left is up
                    if( toss_3<0.5 ) LongRangePairBranching(upSpinSite, state, L, X, Y, XRight, Y, sigma);
                    else LongRangePairBranching(upSpinSite, state, L, XRight, Y, X, Y, sigma);
                }
            }
            else if( toss_2<0.75 ) { 
                if(state[X][YBottom]==1) { // if left is up
                    if( toss_3<0.5 ) LongRangePairBranching(upSpinSite, state, L, X, Y, X, YBottom, sigma);
                    else LongRangePairBranching(upSpinSite, state, L, X, YBottom, X, Y, sigma);
                }
            }
            else {  
                if(state[X][YTop]==1) { // if right is up
                    if( toss_3<0.5 ) LongRangePairBranching(upSpinSite, state, L, X, Y, X, YTop, sigma);
                    else LongRangePairBranching(upSpinSite, state, L, X, YTop, X, Y, sigma);
                }
            }
        }
        else {
            double toss_pp=rnd();
            if( toss_pp<0.25 ) { 
                if(state[XLeft][Y]==0) {
                    state[X][Y] = 0;
                    upSpinSite.erase(upSpinSite.begin()+order); 
                }
            }
            else if( toss_pp<0.5 ) { 
                if(state[XRight][Y]==0) {
                    state[X][Y] = 0;
                    upSpinSite.erase(upSpinSite.begin()+order); 
                }
            }
            else if( toss_pp<0.75 ) { 
                if(state[X][YBottom]==0) {
                    state[X][Y] = 0;
                    upSpinSite.erase(upSpinSite.begin()+order); 
                }
            }
            else {
                if(state[X][YTop]==0) {
                    state[X][Y] = 0;
                    upSpinSite.erase(upSpinSite.begin()+order); 
                }
            }
        }
    }
    else {
        LongRangeCPUpdate(upSpinSite, state, L, p, X, Y, order, sigma);
    }
}

int main(int argc, char *argv[]) {
    const double p = stod(argv[1]); // cp=1.64877, p=0.622466
    const double q = stod(argv[2]); 
    const double sigma = stod(argv[3]); 
    const Time maxStep = stoul(argv[4]);
    const Size L = stoul(argv[5]);
    const Size Nensemble = stoul(argv[6]);
    vector<double> density(maxStep, 0);
    Size NensembelOverSuviving = 0;
    
    for(Size ensemble=0; ensemble<Nensemble; ++ensemble) { // loop over ensemble
        vector<vector<bool>> state(L, vector<bool>(L, 1));
        vector<Size> upSpinSite;
        for(Size i=0; i<L*L; ++i) upSpinSite.push_back(i);
        vector<double> densityPerConf(maxStep, 0);
        for(Size t=0; t<maxStep; ++t) { // temporal loop.
            densityPerConf[t] += upSpinSite.size();
            Size iterNum = upSpinSite.size();
            for(Size trial=0; trial<iterNum; ++trial) { // number of trials during unit time = system size.
                LongRangeTCPUpdate(upSpinSite, state, L, p, q, sigma);
            }
            if(t==maxStep-1) {
                ++NensembelOverSuviving;
                for(Size t=0; t<maxStep; ++t) density[t] += densityPerConf[t];
            }
            if(upSpinSite.size() == 0) {
                break;
            }
        }   
    }
    for(double t=0; t<maxStep; ++t) {
        if(density[t]!=0 && NensembelOverSuviving!=0) cout << t << '\t' << density[t]/L/L/NensembelOverSuviving << endl;
        else break;
    }
    return 0;
}
