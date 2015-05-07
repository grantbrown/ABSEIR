#include <Rcpp.h>
#include <iostream>
#include <sstream>
#include <random>
#include <math.h>
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
                             Eigen::VectorXd offs,
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
                                 offset(offs),
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
    has_reinfection = (reinfection_precision(0) > 0); 
    has_spatial = (I_star.cols() > 1);
    total_size = (has_reinfection && has_spatial ? X.cols() + X_rs.cols() + DM_vec.size() + 2 :
                    (has_reinfection ? X.cols() + X_rs.cols() + 2 : X.cols() + 2));

    aout(this) << "Node Created. has_reinfection: " << has_reinfection 
        << ", has_spatial: " << has_spatial 
        << ", total_size: " << total_size << "\n";

    alive.assign(
        [=](sim_atom, unsigned int param_idx, Eigen::VectorXd param_vals)
        {
            if (total_size != param_vals.size())
            {
                aout(this) << "Invalid parameter vector of length " << param_vals.size() <<  ", ignoring.\n";
                send(parent, param_idx, -2.0);

            }
            else 
            {
                aout(this) << "Valid params observed, simulating.\n";
                double result = simulate(param_vals);
                send(parent, param_idx, result); 
            }
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
    
    int time_idx, i, j;
    
    Eigen::VectorXd beta = params.segment(0, X.cols()); 
    Eigen::VectorXd beta_rs;
    if (has_reinfection) 
    {
        beta_rs = params.segment(X.cols(), X_rs.cols());
    }
    else 
    {
        beta_rs = Eigen::VectorXd(1);
    }

    Eigen::VectorXd rho;
    if (has_spatial && has_reinfection)
    {
        rho = params.segment(X.cols() + X_rs.cols(), DM_vec.size());
    }
    else if (has_spatial)
    {
        rho = params.segment(X.cols(), DM_vec.size());
    }
    else
    {
        rho = Eigen::VectorXd(1.0);
    }

    double gamma_ei = params(params.size() - 2);
    double gamma_ir = params(params.size() - 1);

    Eigen::MatrixXd tmp = (X*beta).unaryExpr([](double elem){return(std::exp(elem));});
    Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor> p_se_components(tmp.data(), 
            I_star.rows(), I_star.cols());

    Eigen::VectorXi N = (S0 + E0 + I0 + R0);

    Eigen::MatrixXi current_S(sim_width, S0.size());
    Eigen::MatrixXi current_E(sim_width, S0.size());
    Eigen::MatrixXi current_I(sim_width, S0.size());
    Eigen::MatrixXi current_R(sim_width, S0.size());

    Eigen::MatrixXi previous_S(sim_width, S0.size());
    Eigen::MatrixXi previous_E(sim_width, S0.size());
    Eigen::MatrixXi previous_I(sim_width, S0.size());
    Eigen::MatrixXi previous_R(sim_width, S0.size());

    Eigen::MatrixXi current_S_star(sim_width, S0.size());
    Eigen::MatrixXi current_E_star(sim_width, S0.size());
    Eigen::MatrixXi current_I_star(sim_width, S0.size());
    Eigen::MatrixXi current_R_star(sim_width, S0.size());

    Eigen::MatrixXi previous_S_star(sim_width, S0.size());
    Eigen::MatrixXi previous_E_star(sim_width, S0.size());
    Eigen::MatrixXi previous_I_star(sim_width, S0.size());
    Eigen::MatrixXi previous_R_star(sim_width, S0.size()); 

    // Initial Case
    time_idx = 0;

    for (i = 0; i < sim_width; i++)
    {
        previous_S.row(i) = S0;
        previous_E.row(i) = E0;
        previous_I.row(i) = I0;
        previous_R.row(i) = R0;
    }

    p_se_components.row(0) = (p_se_components.row(0)).array() * (previous_I.row(0)).cast<double>().array() / N.cast<double>().array();

    Eigen::VectorXd p_se = p_se_components.row(0);

    if (has_spatial)
    {
        for (i = 0; i < DM_vec.size(); i++)
        {
            p_se += rho[i]*(DM_vec[i] * p_se_components.row(0));
        }
    }

    p_se = (-1.0*p_se.array() * offset.array()).unaryExpr([](double e){return(1-std::exp(e));});

    Eigen::VectorXd p_ei = (-1.0*gamma_ei*offset)
                            .unaryExpr([](double e){return(1-std::exp(e));});
    
    Eigen::VectorXd p_ir = (-1.0*gamma_ir*offset)
                            .unaryExpr([](double e){return(1-std::exp(e));}); 

    Eigen::VectorXd p_rs;
    if (has_reinfection)
    {
        p_rs = ((((X_rs*beta_rs).unaryExpr([](double e){return(std::exp(e));})).array() 
                    * offset.array())).unaryExpr([](double e){return(1-std::exp(-e));});
    }
    else
    {
        p_rs = Eigen::VectorXd::Zero(I_star.rows());
    }

    // Initialize Random Draws
    
    std::binomial_distribution<int> S_star_gen(previous_S(0,0) ,p_rs(0));            
    std::binomial_distribution<int> E_star_gen(previous_E(0,0) ,p_se(0));            
    std::binomial_distribution<int> I_star_gen(previous_I(0,0) ,p_ei(0));            
    std::binomial_distribution<int> R_star_gen(previous_R(0,0) ,p_ir(0));            
    time_idx = 0;  
    for (i = 0; i < I_star.cols(); i++)
    {
        for (j = 0; j<sim_width; j++)
        {
            S_star_gen.param(std::binomial_distribution<>::param_type(previous_S(j,i) ,p_rs(0)));
            E_star_gen.param(std::binomial_distribution<>::param_type(previous_E(j,i) ,p_se(0)));
            I_star_gen.param(std::binomial_distribution<>::param_type(previous_I(j,i) ,p_ei(0)));
            R_star_gen.param(std::binomial_distribution<>::param_type(previous_R(j,i) ,p_ir(0)));

            previous_S_star(j,i) = S_star_gen(*generator);
            previous_E_star(j,i) = E_star_gen(*generator);
            previous_I_star(j,i) = I_star_gen(*generator);
            previous_R_star(j,i) = R_star_gen(*generator);
        }
    }

    current_S = previous_S + previous_S_star - previous_E_star;
    current_E = previous_E + previous_E_star - previous_I_star;
    current_I = previous_I + previous_I_star - previous_R_star;
    current_R = previous_R + previous_R_star - previous_S_star;

    // Evaluate Distribution @ previous_S,...,previous_R, 
    //                         previous_S_star, ..., previous_R_star

    previous_S = current_S;
    previous_E = current_E;
    previous_I = current_I;
    previous_R = current_R;

    for (time_idx = 1; time_idx < I_star.rows(); time_idx++)
    {
        for (i = 0; i < I_star.cols(); i++)
        {
            for (j = 0; j<sim_width; j++)
            {
                /*
                std::binomial_distribution<int> S_star_gen(previous_S(j,i) ,p_rs(0));            
                std::binomial_distribution<int> E_star_gen(previous_E(j,i) ,p_se(0));            
                std::binomial_distribution<int> I_star_gen(previous_I(j,i) ,p_ei(0));            
                std::binomial_distribution<int> R_star_gen(previous_R(j,i) ,p_ir(0));            
                */
                S_star_gen.param(std::binomial_distribution<>::param_type(previous_S(j,i) ,p_rs(0)));
                E_star_gen.param(std::binomial_distribution<>::param_type(previous_E(j,i) ,p_se(0)));
                I_star_gen.param(std::binomial_distribution<>::param_type(previous_I(j,i) ,p_ei(0)));
                R_star_gen.param(std::binomial_distribution<>::param_type(previous_R(j,i) ,p_ir(0)));

                previous_S_star(j,i) = S_star_gen(*generator);
                previous_E_star(j,i) = E_star_gen(*generator);
                previous_I_star(j,i) = I_star_gen(*generator);
                previous_R_star(j,i) = R_star_gen(*generator);
            }
        }

        current_S = previous_S + previous_S_star - previous_E_star;
        current_E = previous_E + previous_E_star - previous_I_star;
        current_I = previous_I + previous_I_star - previous_R_star;
        current_R = previous_R + previous_R_star - previous_S_star;

        // Evaluate Distribution @ previous_S,...,previous_R, 
        //                         previous_S_star, ..., previous_R_star

        previous_S = current_S;
        previous_E = current_E;
        previous_I = current_I;
        previous_R = current_R;
    }

    // Dummy workload. 
    return(params.size()*2.0);
}

SEIR_sim_node::~SEIR_sim_node()
{
    delete generator;
}

