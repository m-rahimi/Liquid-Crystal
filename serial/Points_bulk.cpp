/* 
 * File:   Points_bulk.cpp
 * Author: amin
 * 
 * Created on May 1, 2015, 11:04 AM
 */

#include <iostream>
#include "Points_bulk.h"
#include "Parameter.h"
#include "Points.h"
#include <stdlib.h> 
#include <math.h>
#include <vector>

using namespace std;

Points_bulk::Points_bulk() { 
    //srand (time(NULL));
  //  cout << "  " << endl;
 }

Points_bulk::Points_bulk(const Parameter *in) { 
    cout << " Bulk Points construct " << endl;
 }

Points_bulk::Points_bulk(const Points_bulk& orig) {
}

Points_bulk::~Points_bulk() {
}
void Points_bulk::initial_Q (const Parameter *in, Points_bulk * const pbulk){
    Points_bulk *pp = pbulk;
    float vv[3];
    int indx;
    float BI2 = 1.41421356237 * in->qch * in->red;
    for (int ii=0; ii<in->Bpoints; ++ii){
        // Initial Configuration Uniform
        if (strcmp("Unix",in->Init)==0){
            vv[0]=1;vv[1]=0;vv[2]=0;
            Points_bulk::tensor(in->Sinit, vv, pp);}
        else if (strcmp("Uniy",in->Init)==0){
            vv[0]=0;vv[1]=1;vv[2]=0;
            Points_bulk::tensor(in->Sinit, vv, pp);}
        else if (strcmp("Uniz",in->Init)==0){
            vv[0]=0;vv[1]=0;vv[2]=1;
            Points_bulk::tensor(in->Sinit, vv, pp);}
        // Initial Configuration Random
        if (strcmp("rand",in->Init)==0){
            for (int ii=0; ii<3; ++ii){vv[ii] = float(rand()%1000)/999;}
            float dvv = sqrt(vv[0]*vv[0]+vv[1]*vv[1]+vv[2]*vv[2]);
            for (int ii=0; ii<3; ++ii){vv[ii] /= dvv;
            Points_bulk::tensor(in->Sinit, vv, pp);}
        }
        if (strcmp("BI2",in->Init)==0){
            indx = 3*ii;
            float s2 = sqrt(2.0);
            float xx = (float(in->Position[indx]) -   float(in->box[0]-1)/2.0) * BI2;
            float yy = (float(in->Position[indx+1]) - float(in->box[1]-1)/2.0) * BI2;
            float zz = (float(in->Position[indx+2]) - float(in->box[2]-1)/2.0) * BI2;

            pp->QQ[0] = 0.2*(-sin(yy)*cos(xx) - sin(xx)*cos(zz) + 2.0*sin(zz)*cos(yy));
            pp->QQ[1] = 0.2*(-s2*sin(xx)*sin(zz) - s2*cos(yy)*cos(zz) + sin(xx)*cos(yy));
            pp->QQ[2] = 0.2*(-s2*sin(zz)*sin(yy) - s2*cos(xx)*cos(yy) + sin(zz)*cos(xx));
            pp->QQ[3] = 0.2*(-sin(zz)*cos(yy) - sin(yy)*cos(xx) + 2.0*sin(xx)*cos(zz));
            pp->QQ[4] = 0.2*(-s2*sin(yy)*sin(xx) - s2*cos(zz)*cos(xx) + sin(yy)*cos(zz));
            pp->QQ[5] = -pp->QQ[0] - pp->QQ[3];
                    
        }
        ++pp;   
    }
    
}

void Points_bulk::tensor(float Sinit,float vv[3], Points_bulk *pp){
    pp->QQ[0] = Sinit*(vv[0]*vv[0]-1.0/3.0);
    pp->QQ[1] = Sinit*(vv[0]*vv[1]);
    pp->QQ[2] = Sinit*(vv[0]*vv[2]);
    pp->QQ[3] = Sinit*(vv[1]*vv[1]-1.0/3.0);
    pp->QQ[4] = Sinit*(vv[1]*vv[2]);
    pp->QQ[5] = Sinit*(vv[2]*vv[2]-1.0/3.0);
}

void Points_bulk::list(const Parameter *in, Points_bulk * const pbulk){
    Points_bulk *pp, *pr, *pl;
        int box = in->box[0], boy = in->box[1], boz = in->box[2];
        for (int kk=0; kk<box; ++kk){
            for (int jj=0; jj<boy; ++jj){
                for (int ii=0; ii<boz; ++ii){
                    int indx = ii + jj*box + kk*box*boy;
                    pp = pbulk + indx;
 //                   cout << indx << "   " ; //<< pp.QQ[0] << "   ";
                    for (int xx=0; xx<3; ++xx){
                        int pos_r[3] = {ii,jj,kk};
                        ++pos_r[xx];
                        if (pos_r[xx] >= in->box[xx]) pos_r[xx]=0;
                        int indx_r = pos_r[0] + pos_r[1]*box + pos_r[2]*box*boy;
                        pr = pbulk + indx_r;
//                        cout << indx_r << "   ";
                        int pos_l[3] = {ii,jj,kk};
                        --pos_l[xx];
                        if (pos_l[xx] < 0) pos_l[xx]=in->box[xx]-1;
                        int indx_l = pos_l[0] + pos_l[1]*box + pos_l[2]*box*boy;
                        pl = pbulk + indx_l;
 //                       cout << indx_l << "   ";

                        pp->plist[xx*2] = pr;
                        pp->plist[xx*2+1] = pl;
                    }
//                    cout << endl;                    
                }
            }
        }
}
