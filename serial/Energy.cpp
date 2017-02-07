/* 
 * File:   Energy.cpp
 * Author: amin
 * 
 * Created on May 4, 2015, 11:19 AM
 */

#include <iostream>
#include <numeric>
#include <math.h>
#include "Energy.h"
#include "Parameter.h"
#include "Function.h"
#include "Points_bulk.h"

using namespace std; 

Energy::Energy() {
}

Energy::Energy(const Energy& orig) {
}

Energy::~Energy() {
}

void Energy::free(const Parameter* in, Points_bulk * const pbulk){
    LdG=0.0;El1=0.0;chiral=0.0;
    int indx;
    float grad2[18];
    Function fun;
    Points_bulk * pp = pbulk;
    for (int ii=0; ii<in->Bpoints; ++ii){
        //Calculate LdG Energy
        fun.trace(pp);fun.LdG(pp);
        LdG += 0.5*(1.0 - in->UldG/3.0)*fun.tr1 - in->UldG*fun.tr3/3.0 + in->UldG*fun.tr2/4.0;     
        //Calculate Elastic Energy
        fun.gradient(pp,in->dx2);
        for ( indx=0; indx<18; ++indx){grad2[indx] = fun.grad[indx] * fun.grad[indx];}
        for ( indx=0; indx<18; indx+=6){
            El1 += 0.5 * inner_product( grad2+indx, grad2+indx+6, fac, 0.0);
        }
        //Calculate Chiral Energy
        if (in->qch!=0){
            chiral +=(pp->QQ[0] * (-fun.grad[13] + fun.grad[8]) +
                      pp->QQ[1] * (fun.grad[12] - fun.grad[15] - fun.grad[2] + fun.grad[10]) +
                      pp->QQ[2] * (-fun.grad[6] + fun.grad[1] - fun.grad[16] + fun.grad[11]) +
                      pp->QQ[3] * (-fun.grad[4] + fun.grad[13]) +
                      pp->QQ[4] * (-fun.grad[7] + fun.grad[14] + fun.grad[3] - fun.grad[5] ) +
                      pp->QQ[5] * (-fun.grad[8] + fun.grad[4]))*2.0*in->qch;
        }
        ++pp;
    }
    
    dF = fabs((LdG + El1 + chiral - TF)*100/TF);
    TF = LdG + El1 + chiral;  
}

