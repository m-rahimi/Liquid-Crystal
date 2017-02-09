/* 
 * File:   Function.cpp
 * Author: amin
 * 
 * Created on May 4, 2015, 9:05 AM
 */

#include <iostream>
#include <numeric>
#include "Function.h"
#include "Points_bulk.h"

using namespace std;


Function::Function() {
}

Function::Function(const Function& orig) {
}

Function::~Function() {
}

void Function::trace(Points_bulk * const pbulk){
    for (int ii=0; ii<6; ++ii){Q2[ii] = pbulk->QQ[ii] * pbulk->QQ[ii];}
}

void Function::gradient(Points_bulk* const pbulk, float dx2){
    int indx = 0;
    for (int ii=0; ii<6; ii+=2){ // loop for directions
        for (int jj=0; jj<6; ++jj){ // loop for components 
            grad[indx] = (pbulk->plist[ii]->QQ[jj] - pbulk->plist[ii+1]->QQ[jj])*dx2;
            ++indx;
        }
    }
}
void Function::LdG(Points_bulk * const pbulk){
    tr1 = inner_product(Q2,Q2+6,fac,0.0);
    tr2 = tr1 * tr1;
    
    a0 = -pbulk->QQ[0] - pbulk->QQ[3];
    a3 =  pbulk->QQ[0] - pbulk->QQ[3];
    
    tr3 = 0.75*a0*a0*a0 - 3.0*a0*pbulk->QQ[1]*pbulk->QQ[1] + 1.5*a0*pbulk->QQ[2]*pbulk->QQ[2] 
            + 1.5*a3*pbulk->QQ[2]*pbulk->QQ[2] - 0.75*a0*a3*a3 + 6*pbulk->QQ[1]*pbulk->QQ[2]*pbulk->QQ[4]
            + 1.5*a0*pbulk->QQ[4]*pbulk->QQ[4] - 1.5*a3*pbulk->QQ[4]*pbulk->QQ[4];
}

void Function::grad_LdG(Points_bulk* const pbulk, float UU){
    Function::trace(pbulk);
    doubledot = inner_product( Q2, Q2+6, fac, 0.0);
    dotd[0] = pbulk->QQ[0]*pbulk->QQ[0] + pbulk->QQ[1]*pbulk->QQ[1] + pbulk->QQ[2]*pbulk->QQ[2];
    dotd[1] = pbulk->QQ[0]*pbulk->QQ[1] + pbulk->QQ[1]*pbulk->QQ[3] + pbulk->QQ[2]*pbulk->QQ[4];
    dotd[2] = pbulk->QQ[0]*pbulk->QQ[2] + pbulk->QQ[1]*pbulk->QQ[4] + pbulk->QQ[2]*pbulk->QQ[5];
    dotd[3] = pbulk->QQ[1]*pbulk->QQ[1] + pbulk->QQ[3]*pbulk->QQ[3] + pbulk->QQ[4]*pbulk->QQ[4];
    dotd[4] = pbulk->QQ[1]*pbulk->QQ[2] + pbulk->QQ[3]*pbulk->QQ[4] + pbulk->QQ[4]*pbulk->QQ[5];
    dotd[5] = pbulk->QQ[2]*pbulk->QQ[2] + pbulk->QQ[4]*pbulk->QQ[4] + pbulk->QQ[5]*pbulk->QQ[5];
    
    for (int ii=0; ii<6; ++ii){
        dLdG[ii] = (1.0 - UU/3.0) * pbulk->QQ[ii] - UU * dotd[ii] + UU * doubledot * (pbulk->QQ[ii] + delta[ii]);
    }
}

void Function::grad_El1(Points_bulk* const pbulk, float ddx){
    for (int ii=0; ii<6; ++ii){ dEl1[ii] = 0.0;}
    for (int ii=0; ii<6; ii+=2){ // loop for directions
        for (int jj=0; jj<6; ++jj){ // loop for components 
            dEl1[jj] += (pbulk->plist[ii]->QQ[jj] -2.f*pbulk->QQ[jj] + pbulk->plist[ii+1]->QQ[jj])*ddx;
       
        }
    }
}

void Function::grad_chiral(Points_bulk* const pbulk, float qch){
    dchiral[0] = 2.0*qch*(grad[8]- grad[13]); 
    dchiral[1] = qch*(grad[10]-grad[15]-grad[2]+grad[12]);
    dchiral[2] = qch*(grad[11]-grad[16]+grad[1]-grad[6]);
    dchiral[3] = 2.0*qch*(grad[13]-grad[4]);
    dchiral[4] = qch*(grad[14]-grad[5]+grad[3]-grad[7]);
    dchiral[5] = 2.0*qch*(grad[4]-grad[8]);
}



