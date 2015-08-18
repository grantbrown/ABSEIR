#include <iostream>
#include <sstream>
#include <random>
#include <math.h>
#include <Rmath.h>
#include "caf/all.hpp"
#include "SEIRSimNodes.hpp"
#include "spatialSEIRModel.hpp"

using namespace std;
using namespace caf;

SEIR_sim_node::SEIR_sim_node(int w,
                             int sd,
                             Eigen::VectorXi s,
                             Eigen::VectorXi e,
                             Eigen::VectorXi i,
                             Eigen::VectorXi r,
                             Eigen::VectorXd offs,
                             Eigen::MatrixXi y,
                             MatrixXb nm,
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
                             int dmc,
                             actor pr
                             ) : sim_width(w),
                                 random_seed(sd),
                                 S0(s),
                                 E0(e),
                                 I0(i),
                                 R0(r),
                                 offset(offs),
                                 Y(y),
                                 na_mask(nm),
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
                                 data_compartment(dmc),
                                 parent(pr)
{
    try
    {
        std::minstd_rand0 lc_generator(sd);
        std::uint_least32_t seed_data[std::mt19937::state_size];
        std::generate_n(seed_data, std::mt19937::state_size, std::ref(lc_generator));
        std::seed_seq q(std::begin(seed_data), std::end(seed_data));
        generator = new mt19937{q};   
    }
    catch (int e)
    {
        aout(this) << "Error in constructor: " << e << "\n";
    }
    has_reinfection = (reinfection_precision(0) > 0); 
    has_spatial = (Y.cols() > 1);
    total_size = (has_reinfection && has_spatial 
                      ? X.cols() + X_rs.cols() + DM_vec.size() + 2 
               : (has_reinfection 
                      ? X.cols() + X_rs.cols() + 2 
               : (has_spatial 
                      ? X.cols() + DM_vec.size() + 2 
               : X.cols() + 2)));

    alive.assign(
        [=](sim_atom, unsigned int param_idx, Eigen::VectorXd param_vals)
        {
            if (total_size != param_vals.size())
            {
                aout(this) << "Invalid parameter vector of length " 
                 << param_vals.size() <<  ", looking for: " << total_size << ", ignoring.\n"; 
                send(parent, param_idx, -2.0);

            }
            else 
            {
                //aout(this) << "Valid params observed, simulating.\n";
                double result = simulate(param_vals, false).result;
                send(parent, param_idx, result); 
            }
        },
        [=](sample_atom, unsigned int param_idx, Eigen::VectorXd param_vals)
        {
            if (total_size != param_vals.size())
            {
                aout(this) << "Invalid parameter vector of length " 
                 << param_vals.size() <<  ", looking for: " << total_size << ", ignoring.\n"; 
                send(parent, param_idx, -2.0);

            }
            else 
            {
                //aout(this) << "Valid params observed, simulating.\n";
                double result = simulate(param_vals, false).result;
                send(parent, param_idx, result); 
            }
        },
        [=](sim_result_atom, unsigned int param_idx, Eigen::VectorXd param_vals)
        {
            if (total_size != param_vals.size())
            {
                aout(this) << "Invalid parameter vector of length " 
                 << param_vals.size() <<  ", looking for: " << total_size << ", ignoring.\n"; 
                send(parent, param_idx, -2.0);

            }
            else 
            {
                //aout(this) << "Valid params observed, simulating result.\n";
                simulationResultSet result = simulate(param_vals, true);
                send(parent, param_idx, result); 
            }
        },

        [=](exit_atom)
        {
            //aout(this) << "Node quitting.\n";
            quit();
        }
    );
}

behavior SEIR_sim_node::make_behavior(){
    send(this, wakeup_atom::value);
    return([=](wakeup_atom){become(alive);});
}

simulationResultSet SEIR_sim_node::simulate(Eigen::VectorXd params, bool keepCompartments)
{
    // Params is a vector made of:
    // [Beta, Beta_RS, rho, gamma_ei, gamma_ir]    
    int oldWidth = this -> sim_width;
    if (keepCompartments)
    {
        // If we're saving debug information, each set of parameters
        // must only run one simulation. Duplicate sims in R if you
        // need to. 
        this -> sim_width = 1;
    }

    simulationResultSet compartmentResults;
    unsigned int idx; 
    int time_idx, i, j;
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

    Eigen::VectorXi N(S0.size());
    N = (S0 + E0 + I0 + R0);

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



    // Calculate probabilities

    // p_se calculation
    // Equivalent R expression: 
    // exp(matrix(X %*% beta, nrow = nrow(Y), ncol = ncol(Y)))
    Eigen::MatrixXd eta = ((X*beta).unaryExpr([](double elem){return(std::exp(elem));}));
    Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor> p_se_components(eta.data(), 
                Y.rows(), Y.cols());

    time_idx = 0;
    for (i = 0; i < sim_width; i++)
    {
        previous_S.row(i) = S0;
        previous_E.row(i) = E0;
        previous_I.row(i) = I0;
        previous_R.row(i) = R0;

        p_se_cache.row(i) = (((p_se_components.row(0)).array() 
                       * (previous_I.row(i).cast<double>().array())) 
                       * (((N.transpose().cast<double>().array()).unaryExpr([](double e)
                        {
                            return(1.0/e);
                        })).array())).matrix();
    }

    Eigen::MatrixXd p_se = 1*p_se_cache; 
    if (has_spatial)
    {
        for (idx = 0; idx < DM_vec.size(); idx++)
        {
            for (j = 0; j < sim_width; j++)
            {
                p_se.row(j) += rho[idx]*(DM_vec[idx] * (p_se_cache.row(j).transpose()));
            }
        }
    }

    for (i = 0; i < p_se.rows(); i++)
    {
        p_se.row(i) = (((-1.0*p_se.row(i).array()) * (offset(0)))).unaryExpr([](double e){return(1-std::exp(e));}); 
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
                    * offset(0)).matrix()).unaryExpr([](double e){return(1-std::exp(-e));});
    }
    else
    {
        p_rs = Eigen::VectorXd::Zero(Y.rows());
    }

    // Initialize debug info if applicable
    if (keepCompartments)
    {
        compartmentResults.result = 0.0;
        compartmentResults.S = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.E = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.I = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.R = Eigen::MatrixXi(Y.rows(), Y.cols());

        compartmentResults.S_star = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.E_star = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.I_star = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.R_star = Eigen::MatrixXi(Y.rows(), Y.cols());
        
        compartmentResults.X = Eigen::MatrixXd(X.rows(), X.cols());  
        compartmentResults.X = X; 
        compartmentResults.beta = Eigen::MatrixXd(1, beta.size());
        compartmentResults.beta = beta.transpose(); 
        // Todo: add reinfection component
        compartmentResults.p_se = Eigen::MatrixXd(Y.rows(), Y.cols());
        compartmentResults.p_se.row(0) = p_se.row(0);
        compartmentResults.p_ei = p_ei.transpose(); 
        compartmentResults.p_ir = p_ir.transpose(); 
        if (has_spatial)
        {
            compartmentResults.rho = rho.transpose();
        }
        else
        {
            compartmentResults.rho = Eigen::MatrixXd(1,1);
            compartmentResults.rho(0, 0) = 0.0; 
        }

    }


    // Initialize Random Draws
    time_idx = 0;     
    std::binomial_distribution<int> S_star_gen(previous_S(0,0) ,p_rs(0));            
    std::binomial_distribution<int> E_star_gen(previous_E(0,0) ,p_se(0));            
    std::binomial_distribution<int> I_star_gen(previous_I(0,0) ,p_ei(0));            
    std::binomial_distribution<int> R_star_gen(previous_R(0,0) ,p_ir(0));            
    for (i = 0; i < Y.cols(); i++)
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

            results(j) += (na_mask(0,i) ? 0 : 
                    pow((previous_I_star(j,i) - Y(0, i)), 2.0)); 
        }
    }

    if (keepCompartments)
    {
        for (i = 0; i < S0.size(); i++)
        {
            compartmentResults.S(time_idx, i) = S0(i);
            compartmentResults.E(time_idx, i) = E0(i);
            compartmentResults.I(time_idx, i) = I0(i);
            compartmentResults.R(time_idx, i) = R0(i);

            compartmentResults.S_star(time_idx, i) = previous_S_star(0,i);
            compartmentResults.E_star(time_idx, i) = previous_E_star(0,i);
            compartmentResults.I_star(time_idx, i) = previous_I_star(0,i);
            compartmentResults.R_star(time_idx, i) = previous_R_star(0,i);
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

    for (time_idx = 1; time_idx < Y.rows(); time_idx++)
    {
        for (i = 0; i < sim_width; i++)
        {
            p_se_cache.row(i) = (((p_se_components.row(time_idx)).array() 
                       * (previous_I.row(i).cast<double>().array())) 
                       * (((N.transpose().cast<double>().array()).unaryExpr([](double e)
                        {
                            return(1.0/e);
                        })).array())).matrix();
        }
        
        p_se = 1*p_se_cache; 
        if (has_spatial)
        {
            for (idx = 0; idx < DM_vec.size(); idx++)
            {
                for (j = 0; j < sim_width; j++)
                {
                    p_se.row(j) += rho[idx]*(DM_vec[idx] * p_se_cache.row(j).transpose());
                }
            }
        }
        for (i = 0; i < sim_width; i++)
        {
            p_se.row(i) = ((-1.0*p_se.row(i).array() * offset(time_idx)).matrix()).unaryExpr([](double e){return(1-std::exp(e));});
        }
 
        for (i = 0; i < Y.cols(); i++)
        {
            for (j = 0; j < sim_width; j++)
            {
                S_star_gen.param(std::binomial_distribution<>::param_type(previous_R(j,i) ,p_rs(time_idx)));
                E_star_gen.param(std::binomial_distribution<>::param_type(previous_S(j,i) ,p_se(j,i)));
                I_star_gen.param(std::binomial_distribution<>::param_type(previous_E(j,i) ,p_ei(time_idx)));
                R_star_gen.param(std::binomial_distribution<>::param_type(previous_I(j,i) ,p_ir(time_idx)));

                previous_S_star(j,i) = S_star_gen(*generator);
                previous_E_star(j,i) = E_star_gen(*generator);
                previous_I_star(j,i) = I_star_gen(*generator);
                previous_R_star(j,i) = R_star_gen(*generator);
                results(j) += (na_mask(time_idx, i) ? 0 : 
                        pow((previous_I_star(j,i) - Y(time_idx, i)), 2.0)); 
            }
        }

        if (keepCompartments)
        {
            for (i = 0; i < S0.size(); i++)
            {
                compartmentResults.S(time_idx, i) = previous_S(0,i);
                compartmentResults.E(time_idx, i) = previous_E(0,i);
                compartmentResults.I(time_idx, i) = previous_I(0,i);
                compartmentResults.R(time_idx, i) = previous_R(0,i);

                compartmentResults.S_star(time_idx, i) = previous_S_star(0,i);
                compartmentResults.E_star(time_idx, i) = previous_E_star(0,i);
                compartmentResults.I_star(time_idx, i) = previous_I_star(0,i);
                compartmentResults.R_star(time_idx, i) = previous_R_star(0,i);

                compartmentResults.p_se.row(time_idx) = p_se.row(0);
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

    double resultVal = 0.0;
    for (i = 0; i < sim_width; i++)
    {
       resultVal += std::sqrt(results(i)); 
    }
    resultVal /= sim_width;
    if (keepCompartments)
    {
        this -> sim_width = oldWidth;
    }
    compartmentResults.result = resultVal; 
    return(compartmentResults);
}

SEIR_sim_node::~SEIR_sim_node()
{
    delete generator;
}

