#include "mg.hpp"
#include <iostream>
#include <utility>

using namespace std;

const size_t preSmooting = 5;
const size_t coarsestSmooting = 5;
const size_t postSmooting = 5;

void vCycle(UV &phi, UV &f, const IStorage &II, float alpha, size_t level)
{
    //Pre-Smoothing
    checkMultithreading(phi.u.rows(), phi.u.cols());
    for (size_t i = 0; i < preSmooting; i++)
        rbgs(phi, f, II(level), alpha);

    //Compute Residual Error
    UV residual = calcResidual(phi, f, II(level), alpha);

    //Restrict
    residual.restrict();

    UV eps(residual.u.getShape(), 0.0, phi.u.getShape());

    //recursion
    if((residual.u.rows() < 5) || (residual.u.cols() < 5)) {
        for (size_t i = 0; i < coarsestSmooting; i++)
            rbgs(eps, residual, II(level + 1), alpha);
    }
    else {
        vCycle(eps, residual, II, alpha, (level + 1));
    }
    checkMultithreading(phi.u.rows(), phi.u.cols());

    //Prolongation and Correction
    eps.prolongateInPlace();
    phi += eps;

    //Post-Smoothing
    for (size_t i = 0; i < postSmooting; i++)
        rbgs(phi, f, II(level), alpha);
}

void fCycle(UV &phi, UV &f, const IStorage &II, float alpha, size_t level)
{
    //Pre-Smoothing
    checkMultithreading(phi.u.rows(), phi.u.cols());
    for (size_t i = 0; i < preSmooting; i++)
        rbgs(phi, f, II(level), alpha);

    //Compute Residual Error
    UV residual = calcResidual(phi, f, II(level), alpha);

    //Restrict
    residual.restrict();

    UV eps(residual.u.getShape(), 0.0, phi.u.getShape());

    //F-Cycle Recursion
    if((residual.u.rows() < 5) || (residual.u.cols() < 5)) {
        for (size_t i = 0; i < coarsestSmooting; i++)
            rbgs(eps, residual, II(level + 1), alpha);
    }
    else {
        fCycle(eps, residual, II, alpha, (level + 1));
    }
    checkMultithreading(phi.u.rows(), phi.u.cols());

    //Prolongation and Correction
    phi += eps.prolongate();

    //Re-Smoothing
    for (size_t i = 0; i < postSmooting; i++)
        rbgs(phi, f, II(level), alpha);

    //Compute Residual Error
    residual = calcResidual(phi, f, II(level), alpha);

    //Restrict
    residual.restrict();

    //V-Cycle Recursion
    if((residual.u.rows() < 5) || (residual.u.cols() < 5)) {
        for (size_t i = 0; i < coarsestSmooting; i++)
            rbgs(eps, residual, II(level + 1), alpha);
    }
    else {
        vCycle(eps, residual, II, alpha, (level + 1));
    }
    checkMultithreading(phi.u.rows(), phi.u.cols());

    //Prolongation and Correction
    eps.prolongateInPlace();
    phi += eps;

    //Post-Smoothing
    for (size_t i = 0; i < postSmooting; i++)
        rbgs(phi, f, II(level), alpha);
}

void wCycle(UV &phi, UV &f, const IStorage &II, float alpha, size_t level)
{
    //Pre-Smoothing
    checkMultithreading(phi.u.rows(), phi.u.cols());
    for (size_t i = 0; i < preSmooting; i++)
        rbgs(phi, f, II(level), alpha);

    //Compute Residual Error
    UV residual = calcResidual(phi, f, II(level), alpha);

    //Restrict
    residual.restrict();

    UV eps(residual.u.getShape(), 0.0, phi.u.getShape());

    //F-Cycle Recursion
    if((residual.u.rows() < 5) || (residual.u.cols() < 5)) {
        for (size_t i = 0; i < coarsestSmooting; i++)
            rbgs(eps, residual, II(level + 1), alpha);
    }
    else {
        wCycle(eps, residual, II, alpha, (level + 1));
    }
    checkMultithreading(phi.u.rows(), phi.u.cols());

    //Prolongation and Correction
    phi += eps.prolongate();

    //Re-Smoothing
    for (size_t i = 0; i < postSmooting; i++)
        rbgs(phi, f, II(level), alpha);

    //Compute Residual Error
    residual = calcResidual(phi, f, II(level), alpha);

    //Restrict
    residual.restrict();

    //V-Cycle Recursion
    if((residual.u.rows() < 5) || (residual.u.cols() < 5)) {
        for (size_t i = 0; i < coarsestSmooting; i++)
            rbgs(eps, residual, II(level + 1), alpha);
    }
    else {
        wCycle(eps, residual, II, alpha, (level + 1));
    }
    checkMultithreading(phi.u.rows(), phi.u.cols());

    //Prolongation and Correction
    eps.prolongateInPlace();
    phi += eps;

    //Post-Smoothing
    for (size_t i = 0; i < postSmooting; i++)
        rbgs(phi, f, II(level), alpha);
}