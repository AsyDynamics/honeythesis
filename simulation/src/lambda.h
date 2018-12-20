#ifndef LAMBDA_H
#define LAMBDA_H

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip> // set precision

using namespace std;

vector<int> Lambda_Generator(double h, int T, double lambda0, double a, double b, double sigma, double alpha, mt19937 gen);

#endif
