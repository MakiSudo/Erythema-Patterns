#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void main(){
	int i, x, y, N;
	int GridEdge = 200;
	int GridCenter = GridEdge/2;
	double p[200][200];
	double a[200][200];
	double ptemp[200][200];
	double atemp[200][200];

	double t = 0.0;
	double dt;
	double Tmax =1000.0;//Calculation Time
	//Tmax = 1000.0
	N=10000; 
	dt = Tmax / N;//TimeStep
	int div = N / 10;

	double k1p, k1a;
	double k2p, k2a;
	double k3p, k3a;
	double k4p, k4a;
	double fp(double, double, double), p0;
	double fa(double, double, double), a0;
	double Dp = 0.3; //<0.5
	double Da = 0.3;

	double thr = 0.5;
	int Area = 0;
	int range[100];

	FILE *fout, *fout2;
	fout=fopen("AIGreen.csv","w");
	//fout2=fopen("AIGreenAreaFade2.csv","w");

	//InitialCondition
	for(x=0; x<GridEdge; x++){
		for(y=0; y<GridEdge; y++){
			p[x][y] = 0.01;
			a[x][y] = 0.01;
			if(p[x][y] > thr)
				Area++;
		}
	}
	//Output of InitialCondition
	for(x=0;x<GridEdge; x++){
		for(y=0; y<GridEdge; y++){
			fprintf(fout,"%.4f",p[x][y]);
			if(y!=GridEdge-1){
				fprintf(fout," ,");
			}
		}
		fprintf(fout,"\n");
	}

	//Iteration using Runge-Kutta Method
	for (i = 0; i <= N; i++) {

		//Stimulation
		if (i == 1000){
			for(x=0; x<GridEdge; x++){
				for(y=0; y<GridEdge; y++){
					//p[x][y] = 1.0*(exp(-((x-GridCenter)*(x-GridCenter)+(y-GridCenter)*(y-GridCenter))/2500));
					p[x][y] = 1.0*(exp(-((x-70)*(x-70)+(y-60)*(y-60))/750)+exp(-((x-150)*(x-150)+(y-80)*(y-80))/450)+exp(-((x-100)*(x-100)+(y-150)*(y-150))/450));
				
				}
			}
		}

		//BoundaryCondition.TopLeft
		x=0;
		y=0;
			k1p = dt * fp(p[x][y], a[x][y], t);
			k1a = dt * fa(p[x][y], a[x][y], t);
			k2p = dt * fp(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k2a = dt * fa(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k3p = dt * fp(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k3a = dt * fa(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k4p = dt * fp(p[x][y] + k3p, a[x][y] + k3a, t + dt);
			k4a = dt * fa(p[x][y] + k3p, a[x][y] + k3a, t + dt);

			ptemp[x][y] = p[x][y]+((k1p + 2 * k2p + 2 * k3p + k4p) / 6) + dt * Dp * (p[GridEdge-1][0] + p[x+1][0] + p[0][GridEdge-1] + p[0][y+1] - 4*p[0][0]);
			atemp[x][y] = a[x][y]+((k1a + 2 * k2a + 2 * k3a + k4a) / 6) + dt * Da * (a[GridEdge-1][0] + a[x+1][0] + a[0][GridEdge-1] + a[0][y+1] - 4*a[0][0]);


		//BoundaryCondition.BottomLeft
		x=0;
		y=GridEdge-1;
			k1p = dt * fp(p[x][y], a[x][y], t);
			k1a = dt * fa(p[x][y], a[x][y], t);
			k2p = dt * fp(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k2a = dt * fa(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k3p = dt * fp(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k3a = dt * fa(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k4p = dt * fp(p[x][y] + k3p, a[x][y] + k3a, t + dt);
			k4a = dt * fa(p[x][y] + k3p, a[x][y] + k3a, t + dt);

			ptemp[x][y] = p[x][y]+((k1p + 2 * k2p + 2 * k3p + k4p) / 6) + dt * Dp * (p[GridEdge-1][GridEdge-1] +p[x+1][GridEdge-1] + p[0][y-1] + p[0][0] - 4*p[0][GridEdge-1]);
			atemp[x][y] = a[x][y]+((k1a + 2 * k2a + 2 * k3a + k4a) / 6) + dt * Da * (a[GridEdge-1][GridEdge-1] +a[x+1][GridEdge-1] + a[0][y-1] + a[0][0] - 4*a[0][GridEdge-1]);


		//BoundaryCondition.TopRight
		x=GridEdge-1; 
		y=0;
			k1p = dt * fp(p[x][y], a[x][y], t);
			k1a = dt * fa(p[x][y], a[x][y], t);
			k2p = dt * fp(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k2a = dt * fa(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k3p = dt * fp(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k3a = dt * fa(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k4p = dt * fp(p[x][y] + k3p, a[x][y] + k3a, t + dt);
			k4a = dt * fa(p[x][y] + k3p, a[x][y] + k3a, t + dt);

			ptemp[x][y] = p[x][y]+((k1p + 2 * k2p + 2 * k3p + k4p) / 6) + dt * Dp * (p[x-1][0] + p[0][0] + p[GridEdge-1][GridEdge-1] + p[GridEdge-1][y+1] - 4*p[GridEdge-1][0]);
			atemp[x][y] = a[x][y]+((k1a + 2 * k2a + 2 * k3a + k4a) / 6) + dt * Da * (a[x-1][0] + a[0][0] + a[GridEdge-1][GridEdge-1] + a[GridEdge-1][y+1] - 4*a[GridEdge-1][0]);


		//BoundaryCondition.BottomRight
		x=GridEdge-1; 
		y=GridEdge-1;
			k1p = dt * fp(p[x][y], a[x][y], t);
			k1a = dt * fa(p[x][y], a[x][y], t);
			k2p = dt * fp(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k2a = dt * fa(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k3p = dt * fp(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k3a = dt * fa(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k4p = dt * fp(p[x][y] + k3p, a[x][y] + k3a, t + dt);
			k4a = dt * fa(p[x][y] + k3p, a[x][y] + k3a, t + dt);

			ptemp[x][y] = p[x][y]+((k1p + 2 * k2p + 2 * k3p + k4p) / 6) + dt * Dp * (p[x-1][GridEdge-1] + p[0][GridEdge-1] + p[GridEdge-1][y-1] + p[GridEdge-1][0] - 4*p[GridEdge-1][GridEdge-1]);
			atemp[x][y] = a[x][y]+((k1a + 2 * k2a + 2 * k3a + k4a) / 6) + dt * Da * (a[x-1][GridEdge-1] + a[0][GridEdge-1] + a[GridEdge-1][y-1] + a[GridEdge-1][0] - 4*a[GridEdge-1][GridEdge-1]);


		//BoundaryCondition.Left
		x=0;
		for(y=1; y<GridEdge-1; y++){
			k1p = dt * fp(p[x][y], a[x][y], t);
			k1a = dt * fa(p[x][y], a[x][y], t);
			k2p = dt * fp(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k2a = dt * fa(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k3p = dt * fp(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k3a = dt * fa(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k4p = dt * fp(p[x][y] + k3p, a[x][y] + k3a, t + dt);
			k4a = dt * fa(p[x][y] + k3p, a[x][y] + k3a, t + dt);

			ptemp[x][y] = p[x][y]+((k1p + 2 * k2p + 2 * k3p + k4p) / 6) + dt * Dp * (p[GridEdge-1][y] + p[x+1][y] + p[0][y-1] + p[0][y+1] - 4*p[0][y]);
			atemp[x][y] = a[x][y]+((k1a + 2 * k2a + 2 * k3a + k4a) / 6) + dt * Da * (a[GridEdge-1][y] + a[x+1][y] + a[0][y-1] + a[0][y+1] - 4*a[0][y]);
		}


		//BoundaryCondition.Right
		x=GridEdge-1;
		for(y=1; y<GridEdge-1; y++){
			k1p = dt * fp(p[x][y], a[x][y], t);
			k1a = dt * fa(p[x][y], a[x][y], t);
			k2p = dt * fp(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k2a = dt * fa(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k3p = dt * fp(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k3a = dt * fa(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k4p = dt * fp(p[x][y] + k3p, a[x][y] + k3a, t + dt);
			k4a = dt * fa(p[x][y] + k3p, a[x][y] + k3a, t + dt);

			ptemp[x][y] = p[x][y]+((k1p + 2 * k2p + 2 * k3p + k4p) / 6) + dt * Dp * (p[x-1][y] + p[0][y] + p[GridEdge-1][y-1] + p[GridEdge-1][y+1] -4*p[GridEdge-1][y]);
			atemp[x][y] = a[x][y]+((k1a + 2 * k2a + 2 * k3a + k4a) / 6) + dt * Da * (a[x-1][y] + a[0][y] + a[GridEdge-1][y-1] + a[GridEdge-1][y+1] -4*a[GridEdge-1][y]);
		}


		//BoundaryCondition.Top
		y=0;
		for(x=1; x<GridEdge-1; x++){
			k1p = dt * fp(p[x][y], a[x][y], t);
			k1a = dt * fa(p[x][y], a[x][y], t);
			k2p = dt * fp(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k2a = dt * fa(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k3p = dt * fp(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k3a = dt * fa(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k4p = dt * fp(p[x][y] + k3p, a[x][y] + k3a, t + dt);
			k4a = dt * fa(p[x][y] + k3p, a[x][y] + k3a, t + dt);

			ptemp[x][y] = p[x][y]+((k1p + 2 * k2p + 2 * k3p + k4p) / 6) + dt * Dp * (p[x-1][0] + p[x+1][0] + p[x][GridEdge-1] + p[x][y+1] -4*p[x][0]);
			atemp[x][y] = a[x][y]+((k1a + 2 * k2a + 2 * k3a + k4a) / 6) + dt * Da * (a[x-1][0] + a[x+1][0] + a[x][GridEdge-1] + a[x][y+1] -4*a[x][0]);
		}


		//BoundaryCondition.Bottom
		y=GridEdge-1;
		for(x=1; x<GridEdge-1; x++){
			k1p = dt * fp(p[x][y], a[x][y], t);
			k1a = dt * fa(p[x][y], a[x][y], t);
			k2p = dt * fp(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k2a = dt * fa(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
			k3p = dt * fp(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k3a = dt * fa(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
			k4p = dt * fp(p[x][y] + k3p, a[x][y] + k3a, t + dt);
			k4a = dt * fa(p[x][y] + k3p, a[x][y] + k3a, t + dt);
			
			ptemp[x][y] = p[x][y]+((k1p + 2 * k2p + 2 * k3p + k4p) / 6) + dt * Dp * (p[x-1][GridEdge-1] + p[x+1][GridEdge-1] + p[x][y-1] + p[x][0] -4*p[x][GridEdge-1]);
			atemp[x][y] = a[x][y]+((k1a + 2 * k2a + 2 * k3a + k4a) / 6) + dt * Da * (a[x-1][GridEdge-1] + a[x+1][GridEdge-1] + a[x][y-1] + a[x][0] -4*a[x][GridEdge-1]);
		}


		//Others
		for(x=1; x<GridEdge-1; x++){
			for(y=1; y<GridEdge-1; y++){
				k1p = dt * fp(p[x][y], a[x][y], t);
				k1a = dt * fa(p[x][y], a[x][y], t);
				k2p = dt * fp(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
				k2a = dt * fa(p[x][y] + k1p / 2, a[x][y] + k1a / 2, t + dt / 2);
				k3p = dt * fp(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
				k3a = dt * fa(p[x][y] + k2p / 2, a[x][y] + k2a / 2, t + dt / 2);
				k4p = dt * fp(p[x][y] + k3p, a[x][y] + k3a, t + dt);
				k4a = dt * fa(p[x][y] + k3p, a[x][y] + k3a, t + dt);
				
				ptemp[x][y] = p[x][y]+((k1p + 2 * k2p + 2 * k3p + k4p) / 6) + dt * Dp * (p[x-1][y] + p[x+1][y] + p[x][y-1]+ p[x][y+1] - 4*p[x][y]);
				atemp[x][y] = a[x][y]+((k1a + 2 * k2a + 2 * k3a + k4a) / 6) + dt * Da * (a[x-1][y] + a[x+1][y] + a[x][y-1]+ a[x][y+1] - 4*a[x][y]);
			}
		}


		//Update
		for(x=0;x<GridEdge; x++){
			for(y=0; y<GridEdge; y++){
				p[x][y]=ptemp[x][y];
				a[x][y]=atemp[x][y];
			}
		}
		//Output
		if(i%div == 0){
			for(x=0;x<GridEdge; x++){
				for(y=0; y<GridEdge; y++){
					fprintf(fout,"%.4f",p[x][y]);
					if(p[x][y] > thr)
						Area++;
					if(y!=GridEdge-1){
						fprintf(fout," ,");
					}
				}
				fprintf(fout,"\n");
			}
		range[i / div] = Area;
		//fprintf(fout2, "%d\n", range[i / div]);
		printf("%d", i);
		Area = 0;
		}
		t = (i + 1) * dt;
	}
}

double fp(double p, double a, double t){
	double pa = 0.05;
	double qa = 3.4;//
	double ra = 0.8;
	return pa + (qa*p*p/((1+p*p)*(1+a))) - ra*p;
}
double fa(double p, double a, double t){
	double pi = 0.06;//
	double qi = 6.0;
	double ri = 1.0;//Fixed
	
	//return pi + (qi*p*a/((1+p)*(1+a))) - ri*a; //Pattern RED
	//return pi + (qi*p*a*a/((1+p)*(1+a*a))) - ri*a; //Pattern BLUE
	return pi + (qi*p*p*a*a/((1+p*p)*(1+a*a))) - ri*a; //Pattern GREEN
}


