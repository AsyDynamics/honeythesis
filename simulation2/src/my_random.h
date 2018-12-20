#ifndef LAMBDA_H
#define LAMBDA_H

#include <stdlib.h>
#include <math.h>

double gaussrand_NORMAL() {
	static double V1, V2, S;
	static int phase = 0;
	double X;
	if (phase == 0) {
		do {
			double U1 = (double) rand() / RAND_MAX;
			double U2 = (double) rand() / RAND_MAX;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);
		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);
	phase = 1 - phase;
	return X;
}

double gaussrand(double mean, double stdc) {
	return mean + gaussrand_NORMAL() * stdc;
}


double cls_random::randomExponential(double gamma){
    double pV = 0.0;
    while(true)
    {
        pV = (double)rand()/(double)RAND_MAX;
        if (pV != 1){
            break;}
    }
    pV = (-1.0/gamma)*log(1-pV);
    return pV;
}


#endif
