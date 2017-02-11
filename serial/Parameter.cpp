/* 
 * File:   Parameter.cpp
 * Author: amin
 * 
 * Created on April 30, 2015, 7:08 PM
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "Parameter.h"

using namespace std;

Parameter::Parameter(char *input) {
    string str;
    char   seg[256];
    char  *pch;
    for (int ii=0; ii<4; ++ii) strcpy(vtk[ii],"");
    ifstream infile;
    infile.open(input);
    while(!infile.eof()) {
	getline(infile, str);
	if (str!="") {
            strncpy(seg, str.c_str(), sizeof(seg));
            seg[sizeof(seg) - 1] = 0;
            pch = strtok (seg," ,:=");
            string word(pch);
	    if (word == "box") {
		for ( int ii=0; ii<3; ++ii){box[ii] = stoi(strtok (NULL, " ,:="));}}
            else if (word == "dis"){
                for ( int ii=0; ii<3; ++ii){
                    dis[ii] = stof(strtok (NULL, " ,:="));
                    dx=1.0/dis[0];dx2=0.5*dx;ddx=dx*dx;}}
            else if (word == "UldG"){
		UldG = stof(strtok (NULL, " ,:="));
		Sbulk = 0.25 * (1 + 3 * sqrt(1 - 8 / (3 * UldG)));
                Sinit = 0.25 * (1 + 3 * sqrt(1 - 8 / (3 * 2.80)));
        	}
            else if (word == "Init"){
                char *pch = strtok (NULL, " ,:="); 
                strcpy(Init,pch);}
            else if (word == "dt"){
                dt = stof(strtok (NULL, " ,:="));
            }
            else if (word == "dF"){
                dF = stof(strtok (NULL, " ,:="));
            }
            else if (word == "pitch"){
                pitch = stof(strtok (NULL, " ,:="));
                qch = 0.0;
                if (pitch!=0) {qch = 6.28318530718 / pitch;}
            }
            else if (word == "red"){
                red = stof(strtok (NULL, " ,:="));
            }
            else if (word == "Fsam"){
                Fsam = stof(strtok (NULL, " ,:="));
            }
            else if (word == "vtk"){
                char *pch = strtok (NULL, " ,:=");
                int ii=0;
                while (pch != NULL){
                    strcpy(vtk[ii],pch);
                    pch = strtok (NULL, " ,:=");
                    ++ii;
                }       
            }
            
	}
    }
    infile.close();
    initial();
}

void Parameter::initial(){
        int boxx = box[0], boxy = box[1], boxz = box[2];
        Tpoints = boxx*boxy*boxz;
        for (int kk=0; kk<boxx; ++kk){
            for (int jj=0; jj<boxy; ++jj){
                for (int ii=0; ii<boxz; ++ii){
                    Position.push_back(ii);
                    Position.push_back(jj);
                    Position.push_back(kk);
                }
            }
        }
        Bpoints = Position.size()/3;
        Spoints = 0;
}

Parameter::Parameter(const Parameter& orig) {
}

Parameter::~Parameter() {
}

