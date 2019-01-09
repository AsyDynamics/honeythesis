#include <vector>

double mystderr(vector<double> array, double meanvalue){
	double err(0);
	for (int i=0; i<array.size(); i++){
		err += pow(array[i] - meanvalue, 2);
	}
	return  sqrt(err/array.size());
}
