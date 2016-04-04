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


compartment_tap::compartment_tap(Eigen::MatrixXi comp)
{
    compartment = comp;
    idx = 0;
    beenSet = std::vector<int>();
    nLags = comp.rows();
    for (int i = 0; i < nLags; i++)
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
    int proposed = idx - lag; 
    proposed = (proposed < 0 ? nLags + proposed : proposed);
    if (!beenSet[proposed])
    {
        // Error
    }
    return(compartment.row(proposed));
}
