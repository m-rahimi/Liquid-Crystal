#include <iostream>
#include "Points.h"
#include "Points_bulk.h"
#include "Parameter.h"
#include "Eigen.h"
#include <fstream>

using namespace std;

void output(const Parameter* in, Points_bulk * const pbulk){
    float a[3][3];
    float v[3][3];
    float d[3];
    int points = in->Bpoints;
    float *eigen = new float[3*points];
    float *dir   = new float[3*points];
    
    Points_bulk *pp = pbulk;
    for (int ii=0; ii<3*points; ii+=3){
        a[0][0]=pp->QQ[0]; a[0][1]=pp->QQ[1]; a[0][2]=pp->QQ[2];
        a[1][0]=pp->QQ[1]; a[1][1]=pp->QQ[3]; a[1][2]=pp->QQ[4];
        a[2][0]=pp->QQ[2]; a[2][1]=pp->QQ[4]; a[2][2]=pp->QQ[5];
        
        eigen_decomposition(a, v, d);
        
        eigen[ii]=d[2];  eigen[ii+1]=d[1];  eigen[ii+2]=d[2];
        dir[ii]=v[0][2]; dir[ii+1]=v[1][2]; dir[ii+2]=v[2][2];
        ++pp;
    }
    
    cout << "Write vtk file to the bulk.vtk" << endl;
    ofstream myfile ("bulk.vtk");
    if (myfile.is_open()){
        myfile << "# vtk DataFile Version 3.0\n";
        myfile << "vtk output\n";
        myfile << "ASCII\n";
        myfile << "DATASET STRUCTURED_GRID\n";
        myfile << "DIMENSIONS " << in->box[0] << " " << in->box[1] << " " << in->box[2] << endl;
        myfile << "POINTS " << points << "  float" << endl;
        
        for (float kk=0; kk<in->box[2]; kk+=in->dis[2]){
            for (float jj=0; jj<in->box[1]; jj+=in->dis[1]){
                for (float ii=0; ii<in->box[0]; ii+=in->dis[0]){
                    myfile << ii << "  " << jj << "  " << kk << endl;
                }
            }
        }
        
        int nn=0;
        while (strcmp("",in->vtk[nn])!=0){
            if (strcmp("ord",in->vtk[nn])==0){
                myfile << endl;
                myfile << "POINT_DATA " << points << endl;
                myfile << "SCALARS order float\n";
                myfile << "LOOKUP_TABLE default\n";
                for (int ii=0; ii<3*points; ii+=3){
                    myfile << eigen[ii]*1.5 << endl;
                }
            }
        
            if (strcmp("bia",in->vtk[nn])==0){
                myfile << endl;
                myfile << "SCALARS biaxial float\n";
                myfile << "LOOKUP_TABLE default\n";
                for (int ii=0; ii<3*points; ii+=3){
                    myfile << 3.0*eigen[ii+1]+eigen[ii]*1.5 << endl;
                }
            }
        
            if (strcmp("dir",in->vtk[nn])==0){
                myfile << endl;
                myfile << "VECTORS dirZ float\n";
                    for (int ii=0; ii<3*points; ii+=3){
                        if (dir[ii+2]<0.0) dir[ii+2]*=-1;
                        myfile << dir[ii] << "  " << dir[ii+1] << "  " << dir[ii+2] << endl;
            
                    }
                }
            ++nn;
        }
        
        myfile.close();
    }
    else cout << "Unable to open file";
    
    
}
