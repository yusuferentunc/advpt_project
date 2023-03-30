#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include <omp.h>
#include "solver.hpp"
#include "Matrix.hpp"
#include "FlowField.hpp"
#include "ImgDer.hpp"
#include "mg.hpp"

using namespace std;


const float alpha = 1.0f;

// TODO DELETE LATER - DEBUG PURPOSE
#include <stdlib.h>
double getTimeStamp()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.e-9;
}


int main(int argc, char* argv[])
{
	if (argc != 3 && argc != 5)
		cout << "Wrong arguments!" << endl;

	Matrix<float> a(argv[1]);
	Matrix<float> b(argv[2]);

	omp_set_num_threads(6);

	// TODO DELETE LATER - DEBUG PURPOSE
    double startStamp = getTimeStamp();

	//calculate Ix, Iy and It
	IStorage I(a, b);

	//set up vectors
	UV phi (a.getShape(), 0.0, a.getShape());
	UV f (	((I(0).x * I(0).t) * -1.f),
			((I(0).y * I(0).t) * -1.f)	);

	//start calculation
	UV res (a.getShape(), 0.0);
	double resNorm;
	if(true) {
		for(size_t iteration = 0; iteration < 10000; iteration++)
		{
			fCycle(phi, f, I, alpha, 0);

			//norm testing
			res = calcResidual(phi, f, I(0), alpha);
			resNorm = res.u.l2Norm() + res.v.l2Norm();
			std::cout << "residual norm: " << resNorm << "\n";
			if(resNorm < 0.0005) {
				break;
			}
		}
	}
	else {
		for(size_t i = 0; i < 1000000; i++) {
			for(size_t j = 1; j < (phi.u.cols() - 1); j++) {
				for(size_t i = 1; i < (phi.u.rows() - 1); i++) {
					phi.u(i, j) = (	- I(0).x(i, j) * I(0).t(i, j)
									+ alpha * ( phi.u((i + 1), j) + phi.u((i - 1), j) + phi.u(i, (j + 1)) + phi.u(i, (j - 1)) )
									- (I(0).x(i,j) * I(0).y(i,j) * phi.v(i,j))
								)
								/ ((I(0).x(i,j) * I(0).x(i,j)) + (4.0 * alpha));

					phi.v(i, j) = ( - I(0).y(i, j) * I(0).t(i, j)
									+ alpha * ( phi.v((i + 1), j) + phi.v((i - 1), j) + phi.v(i, (j + 1)) + phi.v(i, (j - 1)))
									- (I(0).x(i, j) * I(0).y(i, j) * phi.u(i, j))
								)
								/ ((I(0).y(i, j) * I(0).y(i, j)) + (4.0 * alpha));
				}
			}

			//norm testing
			if(i%100 == 0) {
				res = calcResidual(phi, f, I(0), alpha);
				resNorm = res.u.l2Norm() + res.v.l2Norm();
				std::cout << "residual norm: " << resNorm << "\n";
				if(resNorm < 0.0005) {
					break;
				}
			}
		}
	}


	// TODO DELETE LATER - DEBUG PURPOSE
    double endStamp = getTimeStamp();
    double timeDifference = endStamp-startStamp;
    std::cout << "Total time is "<< timeDifference << std::endl;

	phi.normalize();

	//print
	string nameU = "resultU.bmp";
	string nameV = "resultV.bmp";
	if (argc > 4) {
		nameU = argv[3];
		nameV = argv[4];
	}
	phi.writeToImage(nameU, nameV);

	//compare
	if(false)
		std::cout << "Difference to reference: "
				  << phi.compare(argv[1], true) << " (L2), "
				  << phi.compare(argv[1], false) << "(Linf)" << std::endl;
}
