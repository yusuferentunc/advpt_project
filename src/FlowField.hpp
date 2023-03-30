#pragma once

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <memory>
#include <stdexcept>
#include <vector>
#include <cassert>
#include <utility>

#include "Matrix.hpp"

class UV
{
    public:

        //constructors
        UV() = delete;

        UV(const Matrix<float> &ou, const Matrix<float> &ov) :
            u(ou),
            v(ov)
            { }

        UV(const Matrix<float> &&ou, const Matrix<float> &&ov) :
            u(std::move(ou)),
            v(std::move(ov))
            { }

        UV(std::vector<size_t> shape, float init) :
            u(Matrix<float>(shape[0], shape[1], init)),
            v(Matrix<float>(shape[0], shape[1], init))
            { }

        UV(std::vector<size_t> shape, float init, std::vector<size_t> oldShape) :
            u(Matrix<float>(shape[0], shape[1], init, oldShape)),
            v(Matrix<float>(shape[0], shape[1], init, oldShape))
            { }

        UV& operator+=(const UV& rhs) {
            this->u += rhs.u;
            this->v += rhs.v;
            return *this; 
        }

        friend UV operator+(UV lhs, const UV& rhs) {
            lhs += rhs;
            return lhs;
        }

        Matrix<float> u;
        Matrix<float> v;


        inline void restrict() {
            this->u = std::move(u.restrict());
            this->v = std::move(v.restrict());
        }

        inline void prolongateInPlace() {
            this->u = std::move(u.prolongate());
            this->v = std::move(v.prolongate());
        }

        inline UV prolongate() {
            return UV(std::move(u.prolongate()), std::move(v.prolongate()));
        }

        inline void normalize() {

            float max = -1000;
            float min = 1000;

            float length;
            for(size_t j = 0; j < (u.cols() - 0); j++) {
                for(size_t i = 0; i < (u.rows() - 0); i++) {
                    length = sqrt((u(i, j) * u(i, j)) + (v(i, j) * v(i, j)));
                    if(max < length)
                        max = length;
                    if(min > length)
                        min = length;
                }
            }

            float range = max - min;
            
            for(size_t j = 0; j < (u.cols() - 0); j++) {
                for(size_t i = 0; i < (u.rows() - 0); i++) {
                    u(i, j) = u(i, j) / range;
                }
            }

            for(size_t j = 0; j < (u.cols() - 0); j++) {
                for(size_t i = 0; i < (u.rows() - 0); i++) {
                    v(i, j) = v(i, j) / range;
                }
            }
        }

        inline float compare(std::string path, bool l2 = true) {

            std::string pathRefU = path.substr(0, path.size() - 5) + "ref_u.bmp";
            std::string pathRefV = path.substr(0, path.size() - 5) + "ref_v.bmp";

            if (FILE *file = fopen(pathRefU.c_str(), "r")) {
                fclose(file);
            } else {
                std::cerr << "The file \"" << pathRefU << "\" used for reference comparison does not exist!\n";
                return -1.f;
            }

            if (FILE *file = fopen(pathRefV.c_str(), "r")) {
                fclose(file);
            } else {
                std::cerr << "The file \"" << pathRefV << "\" used for reference comparison does not exist!\n";
                return -1.f;
            }

            Matrix<float> uRef(pathRefU.c_str());
            Matrix<float> vRef(pathRefV.c_str());
            uRef = (uRef * 2.f) + (-1.f);
            vRef = (vRef * 2.f) + (-1.f);
            Matrix<float> uDiff = uRef + (u * (-1.f));
            Matrix<float> vDiff = vRef + (v * (-1.f));

            uDiff.writeToImage("uDiff.bmp");
            vDiff.writeToImage("vDiff.bmp");


            return  (   + cmp(u, uRef, l2) / uRef.norm(l2)
                        + cmp(v, vRef, l2) / vRef.norm(l2)  ) / 2.;
        }
    
        inline void writeToImage(std::string pathU, std::string pathV) {
            u.writeToImage(pathU);
            v.writeToImage(pathV);
        }
};