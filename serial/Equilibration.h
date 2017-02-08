/* 
 * File:   Equilibration.h
 * Author: amin
 *
 * Created on May 5, 2015, 4:41 PM
 */

#ifndef EQUILIBRATION_H
#define	EQUILIBRATION_H

#include "Points.h"
#include "Points_bulk.h"
#include "Parameter.h"

class Equilibration : public Points {
public:
    Equilibration();
    Equilibration(const Equilibration& orig);
    virtual ~Equilibration();
    void GL(const Parameter* in, Equilibration * const pnew, Points_bulk * const pbulk);
private:

};

#endif	/* EQUILIBRATION_H */

