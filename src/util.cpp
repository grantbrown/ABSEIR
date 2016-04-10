#include <iostream>
#include <sstream>
#include <random>
#include <math.h>
#include <Rmath.h>
#include "SEIRSimNodes.hpp"
#include "spatialSEIRModel.hpp"
#include <util.hpp>
#include <chrono>
#include <thread>
using namespace std;


compartment_tap::compartment_tap(int nrow, int ncol)
{
    int i,j;
    compartment = Eigen::MatrixXi(nrow, ncol);
    for (i = 0; i < compartment.rows(); i++)
    {
        for (j = 0; j < compartment.cols(); j++)
        {
            compartment(i,j) = -1;
        }
    }
    idx = 0;
    beenSet = std::vector<int>();
    nLags = compartment.rows();
    for (i = 0; i < nLags; i++)
    {
        beenSet.push_back(0);
    }
}

void compartment_tap::push(Eigen::VectorXi newComp)
{
    compartment.row(idx) = newComp.transpose();
    beenSet[idx] = 1;
    idx += 1;
    idx = (idx >= compartment.rows() ? 0 : idx);
}

Eigen::VectorXi compartment_tap::get(int lag)
{
    int i, j;
    int proposed = idx - lag - 1; 
    int itr = 0;
    while (proposed < 0 && itr < 1e6){
        proposed += nLags;
        itr++;
    }
    if (!beenSet[proposed])
    {
        // Error
    }
    return(compartment.row(proposed));
}
