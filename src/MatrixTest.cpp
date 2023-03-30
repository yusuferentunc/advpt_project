
#include "Matrix.hpp"

void check(bool condition, const std::string& msg)
{
    if (!condition)
    {
        std::cout << "FAILED: " << msg << "\n";
    }
    else
    {
        std::cout << "PASSED: " << msg << "\n";
    }
}
/*
void test_matvec(std::vector< std::pair< bool, std::string > >& results)
{

    Matrix< double > A(5, 5, 1);
    Matrix< double > exA = A.prolongate();

    std::cout << A << std::endl;
    std::cout << exA << std::endl;

    // auto y_comp = matvec(A, x);

    results.push_back({1 == 1, "test_matvec: result equal to expected value"});
}*/

void test_readFromImage(std::vector< std::pair< bool, std::string > >& results);

void test_writeToImage(std::vector< std::pair< bool, std::string > >& results);

void test_prolongate(std::vector< std::pair< bool, std::string > >& results);

void test_computeWeightedSum(std::vector< std::pair< bool, std::string > >& results);

void test_restrict(std::vector< std::pair< bool, std::string > >& results);

void test_assembleI(std::vector< std::pair< bool, std::string > >& results);

// TODO: Decide on operator tests

int main()
{
    std::vector< std::pair< bool, std::string > > results;

    //test_matvec(results);

    size_t passed = 0;
    for (auto [condition, msg] : results)
    {
        check(condition, msg);
        if (condition)
        {
            passed++;
        }
    }

    std::cout << "--- " << passed << "/" << results.size() << " checks passed ---" << std::endl;

    return passed != results.size();
}
