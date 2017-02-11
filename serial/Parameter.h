/* 
 * File:   Parameter.h
 * Author: amin
 *
 * Created on April 30, 2015, 7:08 PM
 */

#ifndef PARAMETER_H
#define	PARAMETER_H
#include <string.h>
#include <vector>

using namespace std;

class Parameter {
public:
    int box[3];
    float dis[3],dx,dx2,ddx;
    float UldG, Sinit, Sbulk;
    char Init[5];
    float dt;
    float dF;
    int Fsam;
    char vtk[4][5];
    int Tpoints, Bpoints, Spoints;
    vector<int> Position;
    float pitch, red, qch;
    Parameter(char *input);
    Parameter();
    Parameter(const Parameter& orig);
    virtual ~Parameter();
private:
    void initial();

};

#endif	/* PARAMETER_H */

