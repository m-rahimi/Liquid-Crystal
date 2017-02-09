/* 
 * File:   Function.h
 * Author: amin
 *
 * Created on May 4, 2015, 9:05 AM
 */

#ifndef FUNCTION_H
#define	FUNCTION_H
#include "Points_bulk.h"
#include "Parameter.h"



class Function {
public:
    float Q2[6];
    float tr1,tr2,tr3;
    float grad[18];
    float doubledot;
    float dotd[6];
    float dLdG[6], dEl1[6], dchiral[6];
    Function();
    Function(const Function& orig);
    virtual ~Function();
    void trace(Points_bulk * const pbulk);
    void gradient(Points_bulk * const pbulk, float dx2);
    void LdG(Points_bulk * const pbulk);
    void grad_LdG(Points_bulk * const pbulk, float UU);
    void grad_El1(Points_bulk * const pbulk, float ddx);
    void grad_chiral(Points_bulk * const pbulk, float qch);
private:
    const float fac[6] = {1.0,2.0,2.0,1.0,2.0,1.0};
    const float delta[6] = {0.333333333, 0.0, 0.0, 0.333333333, 0.0, 0.333333333};
    float a0,a3;
};

#endif	/* FUNCTION_H */

