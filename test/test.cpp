#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

int main(){
	double alpha = 1.5;
	double sum = 0;
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> normal(0,1);
    exponential_distribution<> exp_1(1);
    exponential_distribution<> exp_alpha(alpha);
	for (int i=0; i<50000; i++){
//		cout << normal(gen) << ", " << exp_1(gen) << ", " << exp_alpha(gen) << endl;
//		sumq += 
		sum += exp_alpha(gen);
	}
	cout << sum/50000 << endl;
	return 0;
}
