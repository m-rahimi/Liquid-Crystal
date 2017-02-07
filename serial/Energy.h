/* 
 * File:   Energy.h
 * Author: amin
 *
 * Created on May 4, 2015, 11:19 AM
 */

#ifndef ENERGY_H
#define	ENERGY_H
#include "Parameter.h"
#include "Points_bulk.h"

class Energy {
public:
    float LdG, El1, chiral;
    float TF=0;
    float dF=10;
    Energy();
    Energy(const Energy& orig);
    virtual ~Energy();
    void free(const Parameter *in, Points_bulk * const pbulk);
private:
    const float fac[6] = {1.0,2.0,2.0,1.0,2.0,1.0};

};

#endif	/* ENERGY_H */

