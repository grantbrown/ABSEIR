#include <Rcpp.h>
#include <iostream>
#include <sstream>
#include <random>
#include "caf/all.hpp"
#include "SEIRSimNodes.hpp"

using namespace Rcpp;
using namespace std;
using namespace caf;

SEIR_sim_node::SEIR_sim_node(int w,
                             int sd,
                             std::vector<int> s,
                             std::vector<int> e,
                             std::vector<int> i,
                             std::vector<int> r,
                             Eigen::MatrixXi is,
                             std::vector<Eigen::MatrixXd> dmv,
                             Eigen::MatrixXd x,
                             Eigen::MatrixXd x_rs,
                             actor pr
                             ) : sim_width(w),
                                 random_seed(sd),
                                 S0(s),
                                 E0(e),
                                 I0(i),
                                 R0(r),
                                 I_star(is),
                                 DM_vec(dmv),
                                 X(x),
                                 X_rs(x_rs),
                                 parent(pr)
{
    //generator = new std::mt19936(sd);
    alive.assign(
        [=](sim_atom, unsigned int param_idx, std::vector<double> param_vals)
        {
            double result = simulate(param_vals);
            send(parent, param_idx, result);
        },
        [=](exit_atom)
        {
            quit();
        }
    );
}

behavior SEIR_sim_node::make_behavior(){
    send(this, wakeup_atom::value);
    return([=](wakeup_atom){become(alive);});
}

double SEIR_sim_node::simulate(std::vector<double> params)
{
    // Dummy workload. 
    return(params.size()*2.0);
}

SEIR_sim_node::~SEIR_sim_node()
{
    delete generator;
}

