/* 
 * File:   main.cpp
 * Author: amin
 *
 * Created on April 30, 2015, 6:56 PM
 */

#include <cstdlib>
#include <iostream>
#include <new>
#include <numeric>
#include <vector>
#include <time.h>
#include <stdio.h>
#include "Parameter.h"
#include "Points.h"
#include "Points_bulk.h"
#include "Function.h"
#include "Energy.h"
#include "Equilibration.h"
#include "Output.h"


using namespace std;

int main(int argc, char *argv[]) {
    time_t t0 = time(0);
    char* cht = ctime(&t0);
    cout << "The local date and time is: " << cht ;
    
    char* temp="control";
//    cout << "The Input File :  " << argv[1] << endl;
    cout << "The Input File :  " << temp << endl;
//    Parameter in(argv[1]);
    Parameter in(temp);
        
//    while (argc>2){++argv;--argc;}
    cout << "Total Number of points      :  " << in.Tpoints << endl;
    cout << "Total Number of bulk points :  " << in.Bpoints << endl;
    Points_bulk *bulk = new Points_bulk[in.Bpoints];
    
    bulk->initial_Q(&in, bulk); // initial bulk Q
    in.Position.clear();  // Clear the Position vector 
    bulk->list(&in, bulk); // Make the list of neighbors 
    
    Energy Ef;
    Equilibration *newQ = new Equilibration[in.Bpoints];
    
    int step=0;
    while (Ef.dF > in.dF){
        if (step%in.Fsam==0) {
            Ef.free(&in,bulk);
            printf("%6i %15.5F %15.5f %15.5f %15.5f %15.5f \n", step, Ef.LdG, Ef.El1, Ef.chiral, Ef.TF, Ef.dF);
        }  
        newQ->GL(&in, newQ, bulk);
        ++step;
    }
    
    // Write a vtk file as output
    if (strcmp("",in.vtk[0])!=0) { output(&in, bulk); }
    
    
    // Simulation Time
    time_t t1 = time(0);
    cht = ctime(&t1);
    cout << "The local date and time is: " << cht ; 
    time_t t2 = t1 - t0;
    cout << "Simulation Time: " << t2 << " seconds" << endl;
    struct tm *tt = gmtime ( &t2 );
    printf ("Simulation Time:  %2d:%02d:%02d\n", (tt->tm_hour)%24, tt->tm_min, tt->tm_sec);
    cout << "END PROGRAM" << endl; 
    
    return 0;
}
