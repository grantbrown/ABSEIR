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
                             Eigen::VectorXi s,
                             Eigen::VectorXi e,
                             Eigen::VectorXi i,
                             Eigen::VectorXi r,
                             Eigen::VectorXd offset,
                             Eigen::MatrixXi is,
                             std::vector<Eigen::MatrixXd> dmv,
                             Eigen::MatrixXd x,
                             Eigen::MatrixXd x_rs,
                             Eigen::VectorXd ei_prior,
                             Eigen::VectorXd ir_prior,
                             Eigen::VectorXd se_prec,
                             Eigen::VectorXd rs_prec,
                             Eigen::VectorXd se_mean,
                             Eigen::VectorXd rs_mean,
                             double ph,
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
                                 E_to_I_prior(ei_prior),
                                 I_to_R_prior(ir_prior),
                                 exposure_precision(se_prec),
                                 reinfection_precision(rs_prec),
                                 exposure_mean(se_mean),
                                 reinfection_mean(rs_mean),
                                 phi(ph),
                                 parent(pr)
{
    generator = new mt19937(sd);
    alive.assign(
        [=](sim_atom, unsigned int param_idx, Eigen::VectorXd param_vals)
        {
            aout(this) << "Message Received!\n";
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

double SEIR_sim_node::simulate(Eigen::VectorXd params)
{

    // Params is a vector made of:
    // [Beta, Beta_RS, rho, gamma_ei, gamma_ir]

    // Dummy workload. 
    return(params.size()*2.0);
}

SEIR_sim_node::~SEIR_sim_node()
{
    delete generator;
}

