/* 
 * File:   Points_bulk.h
 * Author: amin
 *
 * Created on May 1, 2015, 11:04 AM
 */

#ifndef POINTS_BULK_H
#define	POINTS_BULK_H
#include "Points.h"
#include "Parameter.h"


class Points_bulk : public Points {
public:
    Points_bulk *plist[6];
    Points_bulk();
    Points_bulk(const Parameter *in);
    Points_bulk(const Points_bulk& orig);
    virtual ~Points_bulk();
    void initial_Q(const Parameter *in, Points_bulk * const pbulk);
    void list(const Parameter *in, Points_bulk * const bulk);
private:
    void tensor(float Sinit,float vv[3], Points_bulk *pp);

};

#endif	/* POINTS_BULK_H */

