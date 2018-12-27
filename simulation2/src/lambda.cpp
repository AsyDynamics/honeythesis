#include "lambda.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip> // set precision

using namespace std;

vector<int> Lambda_Generator(double h, int T, double lambda0, double a, double b, double sigma, double alpha){
    // variable for random distribution
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> normal(0,1);
    exponential_distribution<> exp_1(1);
    exponential_distribution<> exp_alpha(alpha);
    // variable for main function
    double eps = exp_1(gen);
    int J = T/h;
    double temp = 0;
//    vector<double> lambda;
    double lambda[J];
    vector<int> Nt;
    lambda[0] = lambda0;
    // main
    for (int j=1; j<J; j++){
        double lambda_temp = lambda[j-1] + a*(b-lambda[j-1])*h + sigma*sqrt(lambda[j-1])*sqrt(h)*normal(gen);
        lambda[j] = lambda_temp;
        temp += h*lambda[j];
        if (temp>eps){
            double skip = exp_alpha(gen);
            lambda[j] += skip;
            eps = exp_1(gen);
            temp = 0;
            Nt.push_back(j);
        }
    }
    return Nt;
}
