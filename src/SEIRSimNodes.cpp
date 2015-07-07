#include <Rcpp.h>
#include <iostream>
#include <sstream>
#include <random>
#include <math.h>
#include <Rmath.h>
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
                             Eigen::VectorXd sp_prior,
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
                                 spatial_prior(sp_prior),
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
    total_size = (has_reinfection && has_spatial 
                  ? X.cols() + X_rs.cols() + DM_vec.size() + 2 
                  : (has_reinfection 
                      ? X.cols() + X_rs.cols() + 2 
                      : X.cols() + 2));

    alive.assign(
        [=](sim_atom, unsigned int param_idx, Eigen::VectorXd param_vals)
        {
            if (total_size != param_vals.size())
            {
                //aout(this) << "Invalid parameter vector of length " 
                // << param_vals.size() <<  ", ignoring.\n";
                send(parent, param_idx, -2.0);

            }
            else 
            {
                //aout(this) << "Valid params observed, simulating.\n";
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
    int time_idx, i, j, k;
    Eigen::VectorXd results = Eigen::VectorXd::Zero(sim_width); 
    for (i = 0; i < sim_width; i++)
    {
       results(i) = 0.0; 
    }
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

    Eigen::MatrixXd tmp = (X*beta).unaryExpr([](double elem)
            {
                return(std::exp(elem));
            });
    Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor> p_se_components(tmp.data(), 
            I_star.rows(), I_star.cols());

    Eigen::VectorXi N = (S0 + E0 + I0 + R0);
    
    Eigen::MatrixXd p_se_cache(sim_width, S0.size());

    Eigen::MatrixXi current_S(sim_width, S0.size());
    Eigen::MatrixXi current_E(sim_width, S0.size());
    Eigen::MatrixXi current_I(sim_width, S0.size());
    Eigen::MatrixXi current_R(sim_width, S0.size());

    Eigen::MatrixXi previous_S(sim_width, S0.size());
    Eigen::MatrixXi previous_E(sim_width, S0.size());
    Eigen::MatrixXi previous_I(sim_width, S0.size());
    Eigen::MatrixXi previous_R(sim_width, S0.size());

    Eigen::MatrixXi previous_S_star(sim_width, S0.size());
    Eigen::MatrixXi previous_E_star(sim_width, S0.size());
    Eigen::MatrixXi previous_I_star(sim_width, S0.size());
    Eigen::MatrixXi previous_R_star(sim_width, S0.size()); 

    // Begin Simulation: Initial Case and Setup
    time_idx = 0;
    for (i = 0; i < sim_width; i++)
    {
        previous_S.row(i) = S0;
        previous_E.row(i) = E0;
        previous_I.row(i) = I0;
        previous_R.row(i) = R0;

        p_se_cache.row(i) = (p_se_components.row(0)).array() * (previous_I.row(i).cast<double>().array()) / N.cast<double>().array();
    }
    // Calculate probabilities
    
    //p_se_components.row(0) = (p_se_components.row(0)).array() *
    //    (previous_I.row(0)).cast<double>().array() / N.cast<double>().array();


    Eigen::MatrixXd p_se = 1*p_se_cache; 
    if (has_spatial)
    {
        for (i = 0; i < DM_vec.size(); i++)
        {
            for (j = 0; j < sim_width; j++)
            {
                p_se.row(j) += rho[i]*(DM_vec[i] * p_se_cache.row(j));
            }
        }
    }

    for (i = 0; i < p_se.rows(); i++)
    {
        p_se.row(i) = (-1.0*p_se.row(i).array() * offset.array()).unaryExpr([](double e){return(1-std::exp(e));});
        //p_se = (-1.0*p_se.array() * offset.array()).unaryExpr([](double e)
    }

    Eigen::VectorXd p_ei = (-1.0*gamma_ei*offset)
                            .unaryExpr([](double e){return(1-std::exp(e));});
    
    Eigen::VectorXd p_ir = (-1.0*gamma_ir*offset)
                            .unaryExpr([](double e){return(1-std::exp(e));}); 

    Eigen::VectorXd p_rs;
    if (has_reinfection)
    {
        // untested
        p_rs = ((((X_rs*beta_rs).unaryExpr([](double e){return(std::exp(e));})).array() 
                    * offset.array())).unaryExpr([](double e){return(1-std::exp(-e));});
    }
    else
    {
        p_rs = Eigen::VectorXd::Zero(I_star.rows());
    }


    // Initialize Random Draws
    time_idx = 0;     
    std::binomial_distribution<int> S_star_gen(previous_S(0,0) ,p_rs(0));            
    std::binomial_distribution<int> E_star_gen(previous_E(0,0) ,p_se(0));            
    std::binomial_distribution<int> I_star_gen(previous_I(0,0) ,p_ei(0));            
    std::binomial_distribution<int> R_star_gen(previous_R(0,0) ,p_ir(0));            
    for (i = 0; i < I_star.cols(); i++)
    {
        for (j = 0; j<sim_width; j++)
        {
            S_star_gen.param(std::binomial_distribution<>::param_type(previous_R(j,i) ,p_rs(0)));
            E_star_gen.param(std::binomial_distribution<>::param_type(previous_S(j,i) ,p_se(j,i)));
            I_star_gen.param(std::binomial_distribution<>::param_type(previous_E(j,i) ,p_ei(0)));
            R_star_gen.param(std::binomial_distribution<>::param_type(previous_I(j,i) ,p_ir(0)));

            previous_S_star(j,i) = S_star_gen(*generator);
            previous_E_star(j,i) = E_star_gen(*generator);
            previous_I_star(j,i) = I_star_gen(*generator);
            previous_R_star(j,i) = R_star_gen(*generator);

            if (has_reinfection)
            {
                results(j) += R::dbinom(previous_S_star(j,i), previous_R(j,i), p_rs(0), 1);
            }
            results(j) += R::dbinom(previous_E_star(j,i), previous_S(j,i), p_se(i), 1);
            results(j) += R::dbinom(previous_I_star(j,i), previous_E(j,i), p_ei(0), 1);
            results(j) += R::dbinom(previous_R_star(j,i), previous_I(j,i), p_ir(0), 1);
            results(j) += R::dnorm(previous_I_star(j,i), I_star(0,i), 1.0/phi, 1);
        }
    }

    current_S = previous_S + previous_S_star - previous_E_star;
    current_E = previous_E + previous_E_star - previous_I_star;
    current_I = previous_I + previous_I_star - previous_R_star;
    current_R = previous_R + previous_R_star - previous_S_star;

    previous_S = current_S;
    previous_E = current_E;
    previous_I = current_I;
    previous_R = current_R;

    // Simulation: iterative case
    for (time_idx = 1; time_idx < I_star.rows(); time_idx++)
    {
        for (i = 0; i < sim_width; i++)
        {
            p_se_cache.row(i) = (p_se_components.row(time_idx)).array() * (previous_I.row(i).cast<double>().array()) / N.cast<double>().array();
        }
        
        p_se = 1*p_se_cache; 
        if (has_spatial)
        {
            for (i = 0; i < DM_vec.size(); i++)
            {
                for (j = 0; j < sim_width; j++)
                {
                    p_se.row(j) += rho[i]*(DM_vec[i] * p_se_cache.row(j));
                }
            }
        }

        for (i = 0; i < sim_width; i++)
        {
            p_se.row(i) = (-1.0*p_se.row(i).array() * offset.array()).unaryExpr([](double e){return(1-std::exp(e));});
        }
 
        for (i = 0; i < I_star.cols(); i++)
        {
            for (j = 0; j < sim_width; j++)
            {
                double oldResults = results(j);
                double newResult = 0.0;
                S_star_gen.param(std::binomial_distribution<>::param_type(previous_R(j,i) ,p_rs(time_idx)));
                E_star_gen.param(std::binomial_distribution<>::param_type(previous_S(j,i) ,p_se(j,i)));
                I_star_gen.param(std::binomial_distribution<>::param_type(previous_E(j,i) ,p_ei(time_idx)));
                R_star_gen.param(std::binomial_distribution<>::param_type(previous_I(j,i) ,p_ir(time_idx)));

                previous_S_star(j,i) = S_star_gen(*generator);
                previous_E_star(j,i) = E_star_gen(*generator);
                previous_I_star(j,i) = I_star_gen(*generator);
                previous_R_star(j,i) = R_star_gen(*generator);

                if (has_reinfection)
                {
                    results(j) += R::dbinom(previous_S_star(j,i), previous_R(j,i), p_rs(time_idx), 1);
                }
                results(j) += R::dbinom(previous_E_star(j,i), previous_S(j,i), p_se(j,i), 1);
                results(j) += R::dbinom(previous_I_star(j,i), previous_E(j,i), p_ei(time_idx), 1);
                results(j) += R::dbinom(previous_R_star(j,i), previous_I(j,i), p_ir(time_idx), 1);
                results(j) += R::dnorm(previous_I_star(j,i), I_star(time_idx,i), 1.0/phi, 1);
            }
        }

        current_S = previous_S + previous_S_star - previous_E_star;
        current_E = previous_E + previous_E_star - previous_I_star;
        current_I = previous_I + previous_I_star - previous_R_star;
        current_R = previous_R + previous_R_star - previous_S_star;

        previous_S = current_S;
        previous_E = current_E;
        previous_I = current_I;
        previous_R = current_R;
    }

    // Evaluate Priors
    double priorVal = 0.0;
    double resultVal = 0.0;

    for (j = 0; j < sim_width; j ++)
    {
        resultVal += (results(j) / sim_width);
    }

    for (i = 0; i < beta.size(); i++)
    {
        priorVal += R::dnorm(beta(i), exposure_mean(i), 1.0/exposure_precision(i), 1);
    }

    if (has_reinfection)
    {
        for (i = 0; i < beta_rs.size(); i++)
        {
            priorVal += R::dnorm(beta_rs(i), reinfection_mean(i), 1.0/reinfection_precision(i), 1);
        }
    }
    priorVal += R::dgamma(gamma_ei, E_to_I_prior(0), E_to_I_prior(1), 1);
    priorVal += R::dgamma(gamma_ir, I_to_R_prior(0), I_to_R_prior(1), 1);
    if (has_spatial)
    {
        for (i = 0; i < DM_vec.size(); i++)
        {
            priorVal += R::dbeta(rho(i), spatial_prior(0), spatial_prior(1), 1); 
        }
    }
    resultVal += priorVal;
    if (!std::isfinite(resultVal))
    {
        resultVal = -INFINITY;
    }
    return(resultVal);
}

SEIR_sim_node::~SEIR_sim_node()
{
    delete generator;
}

