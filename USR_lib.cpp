#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "USR_lib.h"

double centroids[NUM_CENTROIDS][DIMENSION];

double scalarProduct(double * v1,double * v2){

	double scaProd =0;
	int i;

	for (i=0;i<DIMENSION-1;i++){
		scaProd = scaProd + v1[i]*v2[i];
	}

	scaProd = sqrt(scaProd);
	return scaProd;
}

double * crossProduct(double * v1,double * v2, double * crossPro){

	crossPro[0] = v1[1]*v2[2]-v1[2]*v2[1];
	crossPro[1] = v1[2]*v2[0]-v1[0]*v2[2];
	crossPro[2] = v1[0]*v2[1]-v1[1]*v2[0];

	return crossPro;

}

int getCentroid1(double * xs, double * ys, double * zs, double * qs, int n_atoms){

	int i;
	double sumX = 0 , sumY = 0, sumZ = 0, sumQ = 0;

	for( i = 0; i < n_atoms; i++ ){

		sumX = sumX + xs[i];
		sumY = sumY + ys[i];
		sumZ = sumZ + zs[i];
		sumQ = sumQ + SCALING_CHARGE*qs[i];
	}

	centroids[0][0] = (double) sumX / (double) n_atoms;
	centroids[0][1] = (double) sumY / (double) n_atoms;
	centroids[0][2] = (double) sumZ / (double) n_atoms;
	centroids[0][3] = (double) sumQ / (double) n_atoms;

	return 1;

}

int getCentroid2(double * xs, double * ys, double * zs, double * qs, int n_atoms){

	int i;
	int id_biggest=0;

	double dst_tmp = 0, dst_biggest = 0;
	
	for (i = 0; i< n_atoms; i++){

		dst_tmp = sqrt((xs[i] - centroids[0][0])*(xs[i] - centroids[0][0]) + (ys[i] - centroids[0][1])*(ys[i] - centroids[0][1]) + (zs[i] - centroids[0][2])*(zs[i] - centroids[0][2]) + (SCALING_CHARGE*qs[i] - centroids[0][3])*(SCALING_CHARGE*qs[i] - centroids[0][3]));

		if( dst_tmp > dst_biggest ){

			dst_biggest = dst_tmp;	
			id_biggest = i;

		}
	}

	centroids[1][0] = xs[id_biggest];
	centroids[1][1] = ys[id_biggest];
	centroids[1][2] = zs[id_biggest];
	centroids[1][3] = SCALING_CHARGE*qs[id_biggest];

	return 1;
}

int getCentroid3(double * xs, double * ys, double * zs, double * qs, int n_atoms){

	int i;
	int id_biggest=0;

	double dst_tmp = 0, dst_biggest = 0;
	
	for (i = 0; i< n_atoms; i++){

		dst_tmp = sqrt((xs[i] - centroids[1][0])*(xs[i] - centroids[1][0]) + (ys[i] - centroids[1][1])*(ys[i] - centroids[1][1]) + (zs[i] - centroids[1][2])*(zs[i] - centroids[1][2])
										+ (SCALING_CHARGE*qs[i] - centroids[1][3])*(SCALING_CHARGE*qs[i] - centroids[1][3]));

		if( dst_tmp > dst_biggest ){

			dst_biggest = dst_tmp;	
			id_biggest = i;

		}
	}

	centroids[2][0] = xs[id_biggest];
	centroids[2][1] = ys[id_biggest];
	centroids[2][2] = zs[id_biggest];
	centroids[2][3] = qs[id_biggest]*SCALING_CHARGE;

	return 1;
}

int getCentroid4_5(double * xs, double * ys, double * zs, double * qs, int n_atoms){

	double A[DIMENSION-1],B[DIMENSION-1];
	double moduleA, moduleCrossAB, crossAB[DIMENSION-1],K;
	int i;
	double tmp[3];
	double qMin=0, qMax=0;
	

	for (i =0;i < DIMENSION-1;i++){
		A[i] = centroids[1][i]-centroids[0][i];
		B[i] = centroids[2][i]-centroids[0][i];
	}

	moduleA = sqrt((centroids[1][0]-centroids[0][0])*(centroids[1][0]-centroids[0][0])+(centroids[1][1]-centroids[0][1])*(centroids[1][1]-centroids[0][1])+(centroids[1][2]-centroids[0][2])*(centroids[1][2]-centroids[0][2])+(centroids[1][3]-centroids[0][3])*(centroids[1][3]-centroids[0][3]));
	crossProduct(A,B,crossAB);
	moduleCrossAB = scalarProduct(crossAB,crossAB);

  K	= moduleA/(2*moduleCrossAB);

	tmp[0] = K*crossAB[0];
	tmp[1] = K*crossAB[1];
	tmp[2] = K*crossAB[2];

	for (i=0; i < n_atoms; i++){
		if (SCALING_CHARGE*qs[i] > qMax){
			qMax = SCALING_CHARGE*qs[i];
		}else if(SCALING_CHARGE*qs[i] < qMin){
			qMin = SCALING_CHARGE*qs[i];
		}
	}

	centroids[3][0] = centroids[0][0] + tmp[0];
	centroids[3][1] = centroids[0][1] + tmp[1];
	centroids[3][2] = centroids[0][2] + tmp[2];
	centroids[3][3] = qMax;

	centroids[4][0] = centroids[0][0] + tmp[0];
	centroids[4][1] = centroids[0][1] + tmp[1];
	centroids[4][2] = centroids[0][2] + tmp[2];
	centroids[4][3] = qMin;

	return 1;
	
}

int getCentroids( double * xs, double * ys, double * zs, double * qs, int n_atoms){

	if (getCentroid1( xs, ys, zs, qs, n_atoms) == 0) return 0;
	if (getCentroid2( xs, ys, zs, qs, n_atoms) == 0) return 0;
	if (getCentroid3( xs, ys, zs, qs, n_atoms) == 0) return 0;
	if (getCentroid4_5( xs, ys, zs, qs, n_atoms) == 0) return 0;
	
	return 1;
}


double * calculateMoments(double * x, double * y, double * z, double * q, int n_atoms, double * M){

	int i;
	double d0, d1, d2, d3, d4;

	double sum_d0_1=0,sum_d0_2=0,sum_d0_3=0;
	double sum_d1_1=0,sum_d1_2=0,sum_d1_3=0;
	double sum_d2_1=0,sum_d2_2=0,sum_d2_3=0;
	double sum_d3_1=0,sum_d3_2=0,sum_d3_3=0;
	double sum_d4_1=0,sum_d4_2=0,sum_d4_3=0;
		
	for(i=0; i < n_atoms; i++){

//! Calulating distances to the centroids (if more centroids are needed make this block into a for loop)
		d0 = sqrt( (x[i]-centroids[0][0])*(x[i]-centroids[0][0]) +
						(y[i]-centroids[0][1])*(y[i]-centroids[0][1]) + (z[i]-centroids[0][2])*(z[i]-centroids[0][2]) + (q[i]*SCALING_CHARGE-centroids[0][3])*(q[i]*SCALING_CHARGE-centroids[0][3]));
						
		d1 = sqrt( (x[i]-centroids[1][0])*(x[i]-centroids[1][0]) +
						(y[i]-centroids[1][1])*(y[i]-centroids[1][1]) + (z[i]-centroids[1][2])*(z[i]-centroids[1][2]) + (q[i]*SCALING_CHARGE-centroids[1][3])*(q[i]*SCALING_CHARGE-centroids[1][3]));

		d2 = sqrt( (x[i]-centroids[2][0])*(x[i]-centroids[2][0]) +
						(y[i]-centroids[2][1])*(y[i]-centroids[2][1]) + (z[i]-centroids[2][2])*(z[i]-centroids[2][2]) + (q[i]*SCALING_CHARGE-centroids[2][3])*(q[i]*SCALING_CHARGE-centroids[2][3]));

		d3 = sqrt( (x[i]-centroids[3][0])*(x[i]-centroids[3][0]) +
						(y[i]-centroids[3][1])*(y[i]-centroids[3][1]) + (z[i]-centroids[3][2])*(z[i]-centroids[3][2]) + (q[i]*SCALING_CHARGE-centroids[3][3])*(q[i]*SCALING_CHARGE-centroids[3][3]) );

		d4 = sqrt( (x[i]-centroids[4][0])*(x[i]-centroids[4][0]) +
						(y[i]-centroids[4][1])*(y[i]-centroids[4][1]) + (z[i]-centroids[4][2])*(z[i]-centroids[4][2])  + (q[i]*SCALING_CHARGE-centroids[4][3])*(q[i]*SCALING_CHARGE-centroids[4][3]));
						

//! Acumulators for Caluculating moments
		sum_d0_1 += d0;
		sum_d0_2 += d0*d0; 
		sum_d0_3 += d0*d0*d0;

		sum_d1_1 += d1;
		sum_d1_2 += d1*d1; 
		sum_d1_3 += d1*d1*d1;

		sum_d2_1 += d2;
		sum_d2_2 += d2*d2; 
		sum_d2_3 += d2*d2*d2;

		sum_d3_1 += d3;
		sum_d3_2 += d3*d3; 
		sum_d3_3 += d3*d3*d3;

		sum_d4_1 += d4;
		sum_d4_2 += d4*d4; 
		sum_d4_3 += d4*d4*d4;
	}

//! Moments Calculation
  M[0]  = (double) sum_d0_1/n_atoms;
	M[1]  = (double) (sum_d0_2/n_atoms) - M[0]*M[0];
	M[2]  = (double) (sum_d0_3/n_atoms) - 3*M[0]*M[1] - M[0]*M[0]*M[0];
	
	M[3]  = (double) sum_d1_1/n_atoms;
	M[4]  = (double) (sum_d1_2/n_atoms) - M[3]*M[3];
	M[5]  = (double) (sum_d1_3/n_atoms) - 3*M[3]*M[4] - M[3]*M[3]*M[3];
	
	M[6]  = (double) sum_d2_1/n_atoms;
	M[7]  = (double) (sum_d2_2/n_atoms) - M[6]*M[6];
	M[8]  = (double) (sum_d2_3/n_atoms) - 3*M[6]*M[7] - M[6]*M[6]*M[6];

	M[9]  = (double) sum_d3_1/n_atoms;
	M[10] = (double) (sum_d3_2/n_atoms) - M[9]*M[9];
	M[11] = (double) (sum_d3_3/n_atoms) - 3*M[9]*M[10] - M[9]*M[9]*M[9];

	M[12] = (double) sum_d4_1/n_atoms;
	M[13] = (double) (sum_d4_2/n_atoms) - M[12]*M[12];
	M[14] = (double) (sum_d4_3/n_atoms) - 3*M[12]*M[13] - M[12]*M[12]*M[12];

	if( M[2] < 0)
		M[2] = -1.0*(pow(-1.*M[2],1./3.));
	else
		M[2] = pow(M[2],1./3.);
  if( M[5] < 0)
		M[5] = -1.0*(pow(-1.0*M[5],1.0/3.0));
  else
		M[5] = pow(M[5],1.0/3.0);
  if( M[8] < 0)
		M[8] = -1.0*(pow(-1.0*M[8],1.0/3.0));
  else
		M[8] = pow(M[8],1.0/3.0);
  if( M[11] < 0)
    M[11] = -1.0*(pow(-1.0*M[11],1.0/3.0));
  else
		M[11] = pow(M[11],1.0/3.0);
  if( M[14] < 0)
		M[14] = -1.0*(pow(-1.0*M[14],1.0/3.0));
  else
		M[14] = pow(M[14],1.0/3.0);

	//printf("calculate Moments: %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", M[0], M[1], M[2], M[3],
	//M[4], M[5], M[6], M[7], M[8], M[9], M[10], M[11], M[12], M[13], M[14]);

	return M;
}

double shapeDistance(double * u, double * v){

	int i;
	double tmp, dist;

	tmp=0;

	for (i=0; i < NUM_MOMENTS ;i++){
		tmp += fabs(u[i]-v[i]);
	}

	dist = (double) 1/((double) 1 + (double)((double)1/(double)NUM_MOMENTS)*tmp);

	//printf("Moments u:%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", u[0], u[1], u[2], u[3],
		//u[4], u[5], u[6], u[7], u[8], u[9], u[10], u[11], u[12], u[13], u[14]);
	//printf("Moments v:%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", v[0], v[1], v[2], v[3],
		//v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13], v[14]);
		//printf ("\n\n-------------> %f  \n\n",dist);

	return dist;
}

