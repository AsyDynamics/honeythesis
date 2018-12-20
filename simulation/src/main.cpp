#include <iostream>
#include <vector>
#include <ctime>
#include "lambda.cpp"
using namespace std;

int main(){
    cout << "Enter main\n";

    double a(3), b(0.7), sigma(0.1), alpha(1.5), lambda0(0.7);
    double beta(1), X0(10), c(1.5), h(0.01), path(10000);
    int iter(10), step(5);
    double PT[iter];

// for random distribution
    random_device rd;
    mt19937 gen(rd());
    exponential_distribution<> exp_beta(beta);

    for (int k=0; k<iter; k++){
//        cout << "Enter " << k << " iteration\n";
        clock_t begin = clock();
        int Num(0);
        int T = k*step;
        int J = T/h;
        for (int i=0; i<path; i++){
            vector<int> Nt = Lambda_Generator(h,T,lambda0,a,b,sigma,alpha,gen);
            int pos(0), Ntlen = Nt.size();
            double Xt(X0);
            for (int j=0; j<J; j++){
                Xt += c*h;
                if (pos>=Ntlen) {break;}
                if (j+1 == Nt[pos]){
                    pos ++;
                    Xt -= exp_beta(gen);}
                if (Xt <= 0){
                    Num ++;
                    break;}
            }
        }
        PT[k] = Num/path;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout <<"Current PT: "<< PT[k] << ", Elapsed time: " << elapsed_secs << endl;
    }
	return 0;
}
