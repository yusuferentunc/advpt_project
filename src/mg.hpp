#pragma once

#ifndef MG
#define MG

#include <tuple>
#include <vector>
#include "Matrix.hpp"
#include "FlowField.hpp"
#include "ImgDer.hpp"
#include "solver.hpp"

void vCycle(UV &phi, UV &f, const IStorage &II, float alpha, size_t level);
void fCycle(UV &phi, UV &f, const IStorage &II, float alpha, size_t level);
void wCycle(UV &phi, UV &f, const IStorage &II, float alpha, size_t level);

inline void checkMultithreading(size_t rows, size_t cols)
{
    if(cols < 25 || rows < 25)
        omp_set_num_threads(1);
    else
        omp_set_num_threads(6);
}

#endif