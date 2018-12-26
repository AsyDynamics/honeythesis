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
    double beta(1.5), X0(3), c(1.5), h(0.01);
    int iter(10), step(1), T(400), path(10000);
    double Plambda[iter], Plambda1[iter], Plambda2[iter], Plambda3[iter], Plambda4[iter];
	double v(0.3927), eta(0.2946);
// Q measurement
	double aTran = a-eta*pow(sigma,2);
	double tempTran = 1 + a*eta - pow(eta,2)*pow(sigma,2) / 2;
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
    normal_distribution<> normal(0,1);
    exponential_distribution<> exp_1(1);
    exponential_distribution<> exp_alpha(alpha);
    exponential_distribution<> exp_alphaTran(alphaTran);
    double eps = exp_1(gen);

	ofstream fout("record.csv");

    for (int k = 1; k <= iter; k++){
        clock_t begin = clock();


// Real P measurement
		double lambda0 = k*step;
        int Num(0);
        int J = T/h;

		for (int i = 0; i < path; i++){
			vector<int> Nt = Lambda_Generator(h,T,lambda0,a,b,sigma,alpha);
			int pos(0), j(0), Ntlen(Nt.size());
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
        Plambda[k]  = lambda0;
        Plambda1[k] = Num/path;
		Plambda2[k] = exp(eta*lambda0 - v*X0);
		Plambda3[k] = Plambda2[k]*(beta-v)/beta;
        clock_t middle = clock();
// Real P measurement finished


// Real Q measurement
        double tempMeanExp(0);
        double lambda0Tran = lambda0*tempTran;
        double lambdatau[path]; // double lambda0Tran = lambda0*tempTran
        for (int i = 0; i < path; i++){
    // this line should return lambda and Nt, instead of calling function
            double lambda[J], temp(0);
            lambda[0] = lambda0Tran; // used in i-path loop
            vector<int> Nt; // used in i-path loop
            for (int j=1; j<J; j++){
                lambda[j] = lambda[j-1] + aTran*(bTran-lambda[j-1])*h + sigmaTran*sqrt(lambda[j-1])*sqrt(h)*normal(gen);
                temp += lambda[j];
                if (temp>eps){
                    double skip = exp_alphaTran(gen);
                    lambda[j] += skip;
                    eps = exp_1(gen);
                    temp = 0;
                    Nt.push_back(j+1);
                }
            }
    // lmabda_generation finished
    // insted of j-J for loop, the do-while should be used
            int pos(0), j(0), Ntlen(Nt.size());
            double Xt(X0);
            do {
                Xt += c*h;
                if (pos >= Ntlen){break;}
                if (j+1 == Nt[pos]){
                    pos ++;
                    Xt -= exp_betaTran(gen);}
                j++;
            }
            while (Xt > 0);
            lambdatau[i] = lambda[j];

            tempMeanExp += -m*exp(lambdatau[i]);
        }

        Plambda4[k] = exp( m*lambda0*tempTran - v*X0) * (beta-v) / beta*(alpha-eta)/alpha*tempMeanExp/path; 
// Real Q measurement finished


        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        double elapsed_middle_secs = double(end - middle) / CLOCKS_PER_SEC;
        cout <<"No." << k <<" iteration, PT: "<< Plambda1[k] << ", full time: " << elapsed_secs << << ", half time: " << elapsed_middle_secs << endl;
// output result
//    	fout << k+1 << "," << PT[k] << "," << elapsed_secs << endl;
        fout << Plambda[k] << "," << Plambda1[k]*100 << "," << Plambda2[k]*100 << "," << Plambda[3]*100 << "," << Plambda4[k]*100 << endl;
	}
	fout.close();
	return 0;
}
