/* 
 * File:   Points.h
 * Author: amin
 *
 * Created on May 1, 2015, 9:10 AM
 */

#ifndef POINTS_H
#define	POINTS_H
#include "Parameter.h"


class Points {
public:
    float QQ[6];
    Points();
    Points(const Points& orig);
    virtual ~Points();
    Points& operator= (const Points &pnew);
    
private:

};

#endif	/* POINTS_H */

