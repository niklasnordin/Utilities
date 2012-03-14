#include <iostream>

using namespace std;

#ifndef VECTOR_H
#define VECTOR_H

#define SMALL 1.0e-15
typedef double scalar;
typedef int label;

class Vector
{
private:
    double x_;
    double y_;
    double z_;

public:

    //constructor
    Vector();

    //destructor
    ~Vector();

    //access functions
};

#endif
