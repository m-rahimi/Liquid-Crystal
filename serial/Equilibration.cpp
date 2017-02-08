/* 
 * File:   Equilibration.cpp
 * Author: amin
 * 
 * Created on May 5, 2015, 4:41 PM
 */

#include <iostream>
#include "Equilibration.h"
#include "Points_bulk.h"
#include "Function.h"
#include "Parameter.h"

using namespace std;

Equilibration::Equilibration() {
}

Equilibration::Equilibration(const Equilibration& orig) {
}

Equilibration::~Equilibration() {
}

void Equilibration::GL(const Parameter* in, Equilibration * const pnew, Points_bulk* const pbulk){
    Function fun;
    Equilibration * pn = pnew;
    Points_bulk * pb = pbulk;
    for (int ii=0; ii<in->Bpoints; ++ii){
        fun.grad_LdG(pb,in->UldG); 
        fun.grad_El1(pb,in->ddx);
        
        if (in->qch!=0) {
            fun.gradient(pb,in->dx2);
            fun.grad_chiral(pb,in->qch);
        }
         
        for (int jj=0; jj<6; ++jj){
            pn->QQ[jj] = pb->QQ[jj] - in->dt * ( fun.dLdG[jj] - fun.dEl1[jj] + 2.0f*fun.dchiral[jj]);
        }
        ++pn;
        ++pb;
    } 
    
    //assign new to the old
    pn = pnew;
    pb = pbulk;
    for (int ii=0; ii<in->Bpoints; ++ii){
        for (int jj=0; jj<6; ++jj){
            pb->QQ[jj] = pn->QQ[jj] ;
        }
        ++pn;
        ++pb;
    } 
}

