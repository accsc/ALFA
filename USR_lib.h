//--- DATA AND STRUCTURES ---//
#ifndef _ZZZ_

#define _ZZZ_

#define DIMENSION 4
#define NUM_CENTROIDS 5
#define NUM_MOMENTS 15

#define MAX_BUFFER 1024

#define SCALING_CHARGE 25

extern double centroids[NUM_CENTROIDS][DIMENSION];

int getCentroids( double * xs, double * ys, double * zs, double * qs, int n_atoms );
double scalarProduct(double * v1,double * v2);
double * crossProduct(double * v1,double * v2, double * crossPro);
int getCentroid1(double * xs, double * ys, double * zs, double * qs, int n_atoms);
int getCentroid2(double * xs, double * ys, double * zs, double * qs, int n_atoms);
int getCentroid3(double * xs, double * ys, double * zs, double * qs, int n_atoms);
int getCentroid4_5(double * xs, double * ys, double * zs, double * qs, int n_atoms);
int getCentroids( double * xs, double * ys, double * zs, double * qs, int n_atoms);
double * calculateMoments(double * x, double * y, double * z, double * q, int n_atoms,double * M);
double shapeDistance(double * u, double * v);

#endif
