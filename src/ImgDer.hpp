#pragma once

#include <omp.h>
#include <vector>
#include "Matrix.hpp"

class I
{
    public:
        Matrix<float> x;
        Matrix<float> y;
        Matrix<float> t;

        I() = delete;

        I(const Matrix<float> &x, const Matrix<float> &y, const Matrix<float> &t) :
            x(x), y(y), t(t) {}

        I(const Matrix<float> &&x, const Matrix<float> &&y, const Matrix<float> &&t) :
            x(std::move(x)),
            y(std::move(y)),
            t(std::move(t))
            { }

        I(const Matrix<float> &a, const Matrix<float> &b) :
            x(Matrix<float>(a.rows(), a.cols(), 0.0)),
            y(Matrix<float>(a.rows(), a.cols(), 0.0)),
            t(Matrix<float>(a.rows(), a.cols(), 0.0))
        {
            float left, right;
            #pragma omp parallel for schedule(static) private(left, right)
            for(size_t y = 0; y < (a.cols() - 1); y++)
                for (size_t x = 0; x < (a.rows() - 1); x++)
                {
                    left =   a(x      , (y + 1))
                        - a(x      ,  y     )
                        + a((x + 1), (y + 1))
                        - a((x + 1),  y     );

                    right =  b(x      , (y + 1))
                        - b(x      ,  y     )
                        + b((x + 1), (y + 1))
                        - b((x + 1),  y     );
                    
                    this->x(x, y) = 0.25 * (left + right);
                }


            #pragma omp parallel for schedule(static) private(left, right)
            for(size_t y = 0; y < (a.cols() - 1); y++)
                for (size_t x = 0; x < (a.rows() - 1); x++)
                {
                    left =   a((x + 1),  y     )
                        - a(x      ,  y     )
                        + a((x + 1), (y + 1))
                        - a(x      , (y + 1));

                    right =  b((x + 1),  y     )
                        - b(x      ,  y     )
                        + b((x + 1), (y + 1))
                        - b(x      , (y + 1));
                    
                    this->y(x, y) = 0.25 * (left + right);
                }

            #pragma omp parallel for schedule(static) private(left, right)
            for(size_t y = 0; y < (a.cols() - 1); y++)
                for (size_t x = 0; x < (a.rows() - 1); x++)
                {
                    left =   a(x      ,  y     )
                        + a(x      , (y + 1))
                        + a((x + 1),  y     )
                        + a((x + 1), (y + 1));

                    right =  b(x      ,  y     )
                        + b(x      , (y + 1))
                        + b((x + 1),  y     )
                        + b((x + 1), (y + 1));
                    
                    this->t(x, y) = 0.25 * (-left + right);
                }
        }

        //TODO avoid copy!
        inline I restrict() {
            return I(   std::move(this->x.restrict()),
                        std::move(this->y.restrict()),
                        std::move(this->t.restrict())   );
        }

};

class IStorage
{
    public:
        IStorage() = delete;

        IStorage(const std::vector<I> &i) :
            is(i) {}

        IStorage(const std::vector<I> &&i) :
            is(std::move(i)) {}

        IStorage(const Matrix<float> &a, const Matrix<float> &b) :
            is(std::vector<I>())
        {
            I current(a, b);
            is.push_back(current);

            while((current.x.cols() >= 5) && (current.x.rows() >= 5)) {
                current = std::move(current.restrict());
                is.push_back(current);
            }
        }

        inline const I&
        operator()(size_t index) const {
            return is[index];
        }

        inline I&
        operator()(size_t index) {
            return is[index];
        }

    private:
        std::vector<I> is;

};