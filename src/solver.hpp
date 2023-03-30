#pragma once
#ifndef SOLVER
#define SOLVER

#include <vector>
#include "Matrix.hpp"
#include "FlowField.hpp"
#include "ImgDer.hpp"

#include <omp.h>

using namespace std;


//ITERATIVE SOLVER

inline float iterationFormulaU(const Matrix<float> &u, float v, float Ix, float Iy, float alpha, float f, size_t i, size_t j)
{
    return (    + f
                + alpha * ( u((i + 1), j) + u((i - 1), j) + u(i, (j + 1)) + u(i, (j - 1)) )
                - (Ix * Iy * v)
            )
            / ((Ix * Ix) + (4.0 * alpha));
}

inline float iterationFormulaV(const Matrix<float> &v, float u, float Ix, float Iy, float alpha, float f, size_t i, size_t j)
{
    return  (   + f
                + alpha * ( v((i + 1), j) + v((i - 1), j) + v(i, (j + 1)) + v(i, (j - 1)))
                - (Ix * Iy * u)
            )
            / ((Iy * Iy) + (4.0 * alpha));
}


inline void gaussSeidel(UV &phi, const UV &f, const I &I, float alpha)
{
    for(size_t j = 1; j < (phi.u.cols() - 1); j++)
        for(size_t i = 1; i < (phi.u.rows() - 1); i++) {
            phi.u(i, j) = iterationFormulaU(phi.u, phi.v(i, j), I.x(i, j), I.y(i, j), alpha, f.u(i, j), i, j);
            phi.v(i, j) = iterationFormulaV(phi.v, phi.u(i, j), I.x(i, j), I.y(i, j), alpha, f.v(i, j), i, j);
        }
}           

inline void rbgs(UV &phi, const UV &f, const I &I, float alpha)
{
   //update u
    for(size_t offset = 0; offset < 2; offset++)
    {
        #pragma omp parallel for schedule(static)
        for(size_t j = 1; j < (phi.u.cols() - 1); j++)
            for(size_t i = 1 + ((j + offset) % 2); i < (phi.u.rows() - 1); i += 2) {
                phi.u(i, j) = iterationFormulaU(phi.u, phi.v(i, j), I.x(i, j), I.y(i, j), alpha, f.u(i, j), i, j);
            }
    }

    //update v
    for(size_t offset = 0; offset < 2; offset++)
    {
        #pragma omp parallel for schedule(static)
        for(size_t j = 1; j < (phi.u.cols() - 1); j++)
            for(size_t i = 1 + ((j + offset) % 2); i < (phi.u.rows() - 1); i += 2) {
                phi.v(i, j) = iterationFormulaV(phi.v, phi.u(i, j), I.x(i, j), I.y(i, j), alpha, f.v(i, j), i, j);
            }
    }          
}


//RESIDUAL

inline float residualU(const Matrix<float> &u, float v, float Ix, float Iy, float f, float alpha, size_t i, size_t j)
{
    return  + f
            - ((Ix * Ix) + (4.0 * alpha)) * u(i, j)
            + alpha * ( u((i + 1), j) + u((i - 1), j) + u(i, (j + 1)) + u(i, (j - 1)) )
            - (Ix * Iy * v);
}

inline float residualV(const Matrix<float> &v, float u, float Ix, float Iy, float f, float alpha,  size_t i, size_t j)
{
    return  + f
            - ((Iy * Iy) + (4.0 * alpha)) * v(i, j)
            + alpha * ( v((i + 1), j) + v((i - 1), j) + v(i, (j + 1)) + v(i, (j - 1)) )
            - (Ix * Iy * u);
}


inline UV calcResidual(const UV &phi, const UV &f, const I &I, float alpha)
{
    UV res(phi.u.getShape(), 0.0);

    #pragma omp parallel for schedule(static)
    for(size_t j = 1; j < (phi.u.cols() - 1); j++)
        for(size_t i = 1; i < (phi.u.rows() - 1); i++) {
            res.u(i, j) = residualU(phi.u, phi.v(i, j), I.x(i, j), I.y(i, j), f.u(i, j), alpha, i, j);
            res.v(i, j) = residualV(phi.v, phi.u(i, j), I.x(i, j), I.y(i, j), f.v(i, j), alpha, i, j);
        }

    return res;
}

#endif