/* 
 * File:   Points.cpp
 * Author: amin
 * 
 * Created on May 1, 2015, 9:10 AM
 */

#include <iostream>
#include "Points.h"

using namespace std;

Points::Points() { 
}

Points::Points(const Points& orig) {
}

Points::~Points() {
}

Points& Points::operator =(const Points& pnew){
    this->QQ[0] = pnew.QQ[0];
    return *this;
         
}





