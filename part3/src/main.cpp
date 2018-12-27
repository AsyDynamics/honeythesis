
#include <iostream>
#include <vector>
#include <ctime>
#include "lambda.cpp"
#include <fstream>
//#include "my_random.h"
using namespace std;

int main(){
    cout << "Enter main\n";
    double a(3), b(0.7), sigma(0.1), alpha(1.5);
    double beta(1), X0(10), c(1.5), h(0.01), path(10000);
    int iter(10), step(1), T(100);
    double Plambda[iter], Plambda1[iter], Plambda2[iter], Plambda3[iter], Plambda4[iter];
	double v(0.2771), eta(0.1979);
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

    for (int k = 0; k < iter; k++){
        clock_t begin = clock();


// Real P measurement
		double lambda0 = (k+1)*step;
        int Num(0);
        int J = T/h;

		for (int i = 0; i < path; i++){
			vector<int> Nt = Lambda_Generator(h,T,lambda0,a,b,sigma,alpha);
			int pos(0), Ntlen(Nt.size());
			double Xt(X0);

            for (int j=0; j<J; j++){
                Xt += c*h;
                if (pos >= Ntlen){break;}
                if (j == Nt[pos]){
                    pos ++;
                    Xt -= exp_beta(gen);}
                if (Xt <= 0){
                    Num ++;
                    break;}
            }
        }
//		cout << "Num: " << Num <<"," << Num/path <<endl;
        Plambda[k]  = lambda0;
        Plambda1[k] = Num/path;
		Plambda2[k] = exp(eta*lambda0 - v*X0);
		Plambda3[k] = exp(eta*lambda0 - v*X0)*(beta-v)/beta;
        clock_t middle = clock();
// Real P measurement finished


// Real Q measurement
        double tempMeanExp(0);
        double lambda0Tran = lambda0*tempTran;
        double lambdatau[int(path)]; // double lambda0Tran = lambda0*tempTran

    // lmabda_generation finished
    // insted of j-J for loop, the do-while should be used
    	for (int i=0; i<path; i++){
            double Xt(X0), temp(0), lambda(lambda0Tran);
            do {
		        lambda = lambda + aTran*(bTran-lambda)*h + sigmaTran*sqrt(lambda)*sqrt(h)*normal(gen);
                Xt += c*h;
		        temp += h*lambda;
		        if (temp>eps){
			        double skip = exp_alphaTran(gen);
			        Xt -= exp_betaTran(gen);
			        if (Xt <= 0){
                        tempMeanExp += exp(-m*lambda);
				        lambdatau[i] = lambda;
				        break;}
                    else {
			        lambda += skip;
			        eps = exp_1(gen);
			        temp = 0;
                    }
                }
// whether record lambdatau here
            }while (Xt > 0);
        }
        Plambda4[k] = exp( m*lambda0*tempTran - v*X0) * (beta-v) / beta*(alpha-eta)/alpha*tempMeanExp/path;
// Real Q measurement finished
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        double elapsed_middle_secs = double(end - middle) / CLOCKS_PER_SEC;
//        cout <<"No." << k <<" iter, PT: "<< Plambda1[k] << ", Q time: " << elapsed_secs - elapsed_middle_secs << ", P time: " << elapsed_middle_secs << endl;
        cout <<"No." << k <<" iteration, lambda 1,2,3,4: "<< Plambda1[k]<< "," << Plambda2[k] <<"," << Plambda3[k]<< "," << Plambda4[k] << ", full time: " << elapsed_secs << endl;
// output result
//    	fout << k+1 << "," << PT[k] << "," << elapsed_secs << endl;
        fout << Plambda[k] << "," << Plambda1[k]*100 << "," << Plambda2[k]*100 << "," << Plambda[3]*100 << "," << Plambda4[k]*100 << endl;
//        fout << Plambda[k] << "," << Plambda1[k]*100 << "," << Plambda2[k]*100 << "," << Plambda[3]*100 << endl;

	}
	fout.close();
	return 0;
}
