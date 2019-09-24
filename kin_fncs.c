/*Chudamani Aryal
	1000692493
	Due: Feb 28, 2012
	Question 1 and 2
*/
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415927
#endif

void matrixMul(double X[4][4], double Y[4][4], double W[4][4]);
void copyOut(double J[4][4], double K[4][4] );

fwd_kin(theta, x)
double theta[6];
double x[3];
{
double Dx0[4][4] = {{1,0,0,0.15},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
double Rx0[4][4] = {{1,0,0,0},{0,cos(theta[4]),-sin(theta[4]),0}, {0,sin(theta[4]),cos(theta[4]),0},{0,0,0,1}};
double Dz0[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,-0.05},{0,0,0,1}};
double Ry0[4][4] = {{cos(theta[3]),0,sin(theta[3]),0},{0,1,0,0},{-sin(theta[3]),0,cos(theta[3]),0},{0,0,0,1}};
double Dx1[4][4] = {{1,0,0,0.2},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
double Ry1[4][4] = {{cos(theta[2]),0,sin(theta[2]),0},{0,1,0,0},{-sin(theta[2]),0,cos(theta[2]),0},{0,0,0,1}};
double Dx2[4][4] = {{1,0,0,0.3},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
double Ry2[4][4] = {{cos(theta[1]),0,sin(theta[1]),0},{0,1,0,0},{-sin(theta[1]),0,cos(theta[1]),0},{0,0,0,1}};
double Dy0[4][4] = {{1,0,0,0},{0,1,0,0.05},{0,0,1,0},{0,0,0,1}};
double Dz1[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0.25},{0,0,0,1}};
double Ry3[4][4] = {{cos(theta[0]),-sin(theta[0]),0,0},{sin(theta[0]),cos(theta[0]),0,0},{0,0,1,0},{0,0,0,1}};
double result[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
double backup[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};

matrixMul(Dx0,Rx0,result);
copyOut(result, backup);

matrixMul(backup,Dz0,result);
copyOut(result, backup);

matrixMul(backup,Ry0,result);
copyOut(result, backup);

matrixMul(backup,Dx1,result);
copyOut(result, backup);

matrixMul(backup,Ry1,result);
copyOut(result, backup);

matrixMul(backup,Dx2,result);
copyOut(result, backup);

matrixMul(backup,Ry2,result);
copyOut(result, backup);

matrixMul(backup,Dy0,result);
copyOut(result, backup);

matrixMul(backup,Dz1,result);
copyOut(result, backup);

matrixMul(backup,Ry3,result);

x[0] = result[0][3];
x[1] = result[1][3];
x[2] = result[2][3];

}


void matrixMul(double A[4][4], double B[4][4], double C[4][4])
{

int i, j, k;
	
	for (i=0; i<4;i++)
	{
		for (j = 0;j<4;j++)
		{
			C[i][j] = 0;
		}
	}

for (i=0 ; i<4; ++i)
	{
		for (j=0 ; j<4; ++j)
		{
			for (k=0; k<4; ++k)
			{
				C[j][i] += A[k][i] * B[j][k];
			}
		}
		
	}
}

void copyOut(double D[4][4], double E[4][4] )
{
	int a,b;
	for (a=0; a<4;a++)
	{
		for (b = 0;b<4;b++)
		{
			E[a][b] = D[a][b];
		}
	}
}



inv_kin(x, theta)
double x[3];
double theta[6];
{
double L0 = 0.25;
double L1 = 0.3;
double L2 = 0.2;
double L3 = 0.15;
double d1 = 0.05;
double d2 = 0.05;

double alpha = atan2(x[1], x[0]);
double beta = atan2(d1, sqrt(pow(x[0], 2) + pow(x[1], 2) - pow(d1, 2)));
theta[0] = alpha - beta;

double x1 = x[0] + d1* sin(theta[0]) + d2 * cos(theta[0]);
double y1 = x[1] - d1* cos(theta[0]) + d2 * sin(theta[0]);
double z1 = x[2] + L3 - L0;

double gamma = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2));

double phi = acos((pow(gamma,2) - pow(L1,2) - pow(L2,2)) / (-2*L1*L2));

theta[2] = M_PI - phi;

 
double z2;
double epsilon, aama;
if((x[2] + L3) > L0)
{
z2 = x[2] + L3 - L0;
 epsilon = atan2(z2, sqrt(pow(gamma,2)-pow(z2,2)));
 aama = acos((pow(L2, 2) - pow(L1, 2) - pow(gamma, 2)) / (-2*L1*gamma));
theta[1] = (epsilon + aama)*-1;
}
else
{
z2 = L0 - x[2] - L3;
epsilon = atan2(z2, sqrt(pow(gamma,2)-pow(z2,2)));
 aama = acos((pow(L2, 2) - pow(L1, 2) - pow(gamma, 2)) / (-2*L1*gamma));
theta[1] = (epsilon - aama)*1;
}

theta[3] = M_PI/2 - theta[2] -theta[1];


}








