#include <iostream>
#include <cmath>
#include <fstream>
#include <math.h>
using namespace std;
const int dim=2;
void stepfunc (double* p, double* q, const double dt, double& H);

int main ()
{
		
	//const double pi=M_PI;
	const double e=0.6;
	const double dt= 0.05;
	const double tend= 20*M_PI;
	double t=0.0;
	double H;
	double p[dim];
	double q[dim];

		q[0]=1.0-e;
		q[1]=0;
		p[0]=0;
		p[1]=sqrt((1.0+e)/(1.0-e));

			H=0.5*(pow(p[0],2)+pow(p[1],2))-(1.0/(sqrt(pow(q[0],2)+(pow(q[1],2)))));

ofstream out("out");
//out<<0<<"\t"<<q[0]<<"\t"<<q[1]<<"\t"<< H<<endl;

		while(t<=tend)
{
	stepfunc(p,q,dt,H);
	out<<t<<"\t"<<q[0]<<"\t"<<q[1]<<"\t"<< H<<endl;
	t+=dt;
}
out.close();
return 0;
}
void stepfunc (double* p, double* q, const double dt, double& H)
{
	double rad=sqrt(pow(q[0],2)+(pow(q[1],2)));
	for (int i=0; i<dim;i++)
	{	
	p[i]-=dt*q[i]/(pow(rad,3));
	q[i]+=dt*p[i];
	}
H=0.5*(pow(p[0],2)+pow(p[1],2))-1.0/(sqrt(pow(q[0],2)+(pow(q[1],2))));
}


