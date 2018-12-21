#include <iostream>
#include <vector>
#include <ctime>
#include "lambda.cpp"
#include <fstream>
//#include "my_random.h"
using namespace std;

int main(){
    cout << "Enter main\n";
    double a(2), b(1), sigma(0.1), alpha(2);
    double beta(1.5), X0(3), c(1.5), h(0.01), path(10000);
    int iter(10), step(1), T(400);
    double Plambda1[iter], Plambda2[iter], Plambda3[iter], Plambda4[iter];
	double v(0.3927), eta(0.2946);
// Q measurement
	double aTran = a-eta*pow(sigma,2);
	double tempTran = 1+a*eta-eta^2*sigma^2/2;
	double bTran = a*b*tempTran / (a-eta*pow(sigma,2));
	double sigmaTran = sigma * sqrt(tempTran);
	double alphaTran = (alpha-eta) / tempTran;
	double betaTran = beta - v;
	double m = eta/tempTran;
// for random distribution
    random_device rd;
    mt19937 gen(rd());
    exponential_distribution<> exp_beta(beta);
	exponential_distribution<> exp_betaTran(betaTran);
	ofstream fout("record.csv");

    for (int k = 1; k <= iter; k++){
        clock_t begin = clock();
// P measurement
		double lambda0 = k*step;
        int Num(0);
        int T = k*step;
        int J = T/h;
        for (int i = 0; i < path; i++){
            vector<int> Nt = Lambda_Generator(h,T,lambda0,a,b,sigma,alpha);
            int pos(0), Ntlen = Nt.size();
            double Xt(X0);
            for (int j=0; j<J; j++){
                Xt += c*h;
                if (pos >= Ntlen){break;}
                if (j+1 == Nt[pos]){
                    pos ++;
                    Xt -= exp_beta(gen);}
                if (Xt <= 0){
                    Num ++;
                    break;}
            }
        }

        Plambda1[k] = Num/path;
		Plambda2[k] = exp(eta*lambda0 - v*X0);
		Plambda3[k] = Plambda2[k]*(beta-v)/beta;

// Q measurement
		double lambda0Tran = lambda0*tempTran;
		double lambdatau[path];
		for (int i = 0; i < path; i++){
			vector<int> Nt = Lambda_Generator(h,T,lambda0Tran,aTran,bTran,sigmaTran,alphaTran);
			int pos(0), j(0), Ntlen(Nt.size());
			double Xt(X0);
			do{
				Xt += c*h;
				if (pos >= Ntlen){break;}
				if (j+1 == Nt[pos]){
					pos ++;
					Xt -= exp_betaTran(gen);}
				j++;
			}
			while (Xt > 0);
			lambdatau[i] = lambda[];
		}
		Plambda4[k] = ;


        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout <<"No." << k <<" iteration, PT: "<< PT[k] << ", Elapsed time: " << elapsed_secs << endl;
    	fout << k+1 << "," << PT[k] << "," << elapsed_secs << endl;
	}
	fout.close();
	return 0;
}
