#ifndef TYPES_H
#define TYPES_H

struct details
{
    label featurei;
    label pointi;
    label pointo;
    scalar distance;
    scalar distanceo;
    scalar lambda;
    scalar lambdaOther;
    scalar cos;
    vector p;
    vector po;
};

struct edgeDetails
{
    details d0;
    details d1;
};

#endif
