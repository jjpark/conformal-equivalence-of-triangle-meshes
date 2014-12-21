#include "main.hpp"

using namespace std;
using namespace Eigen;

void vectorXdToArray (VectorXd v, GLfloat* g) {
  int size = v.size();
  for (int i = 0; i < size; i++) {
    g[i] = v(i);
  }
}

void vectorToArray (vector<float> v, GLfloat*g) {
  int size = v.size();
  for (int i = 0; i < size; i++) {
    g[i] = v.at(i);
  }
}

VectorXd ArrayToVectorXd (GLfloat *a, int n) {
	VectorXd v(n);
	for (int i = 0; i < n; i++) {
		v(i) = a[i];
	}
	return v;
}

VectorXd vectorToVectorXd (vector<float> vec) {
	VectorXd v(3);
	v << vec[0],vec[1],vec[2];
	return v;
}
