#include <iostream>
#include <sstream>
#include <random>
#include <math.h>
#include <Rmath.h>
#include "caf/all.hpp"
#include "SEIRSimNodes.hpp"
#include "spatialSEIRModel.hpp"
#include <chrono>
#include <thread>
using namespace std;
using namespace caf;

SEIR_sim_node::SEIR_sim_node(int sd,
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
                             std::string mode,
                             Eigen::MatrixXd ei_prior,
                             Eigen::MatrixXd ir_prior,
                             double avgH,
                             Eigen::VectorXd sp_prior,
                             Eigen::VectorXd se_prec,
                             Eigen::VectorXd rs_prec,
                             Eigen::VectorXd se_mean,
                             Eigen::VectorXd rs_mean,
                             double ph,
                             int dmc,
                             bool cmltv,
                             actor pr
                             ) : random_seed(sd),
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
                                 transitionMode(mode),
                                 E_to_I_prior(ei_prior),
                                 I_to_R_prior(ir_prior),
                                 avg_hazard(avgH),
                                 spatial_prior(sp_prior),
                                 exposure_precision(se_prec),
                                 reinfection_precision(rs_prec),
                                 exposure_mean(se_mean),
                                 reinfection_mean(rs_mean),
                                 phi(ph),
                                 data_compartment(dmc),
                                 cumulative(cmltv),
                                 parent(pr)
{
    try
    {
        std::minstd_rand0 lc_generator(sd);
        std::uint_least32_t seed_data[std::mt19937::state_size];
        std::generate_n(seed_data, std::mt19937::state_size, std::ref(lc_generator));
        std::seed_seq q(std::begin(seed_data), std::end(seed_data));
        generator = new mt19937{q};   

        if (transitionMode == "path_specific")
        {
            E_paths = Eigen::MatrixXi(E_to_I_prior.rows(),Y.cols());
            I_paths = Eigen::MatrixXi(I_to_R_prior.rows(),Y.cols());
        }
        else
        {
            E_paths = Eigen::MatrixXi(1,Y.cols());
            I_paths = Eigen::MatrixXi(1,Y.cols());
        }

        if (phi > 0)
        {   
            // We take the floor of the resulting continuous normal, so shift by 0.5
            overdispersion_distribution = std::normal_distribution<double>(0.5, 
                    1.0/phi);
        }
    }
    catch (int e)
    {
        aout(this) << "Error in constructor: " << e << "\n";
    }
    has_reinfection = (reinfection_precision(0) > 0); 
    has_spatial = (Y.cols() > 1);
    exp_transition = (transitionMode == "exponential");
    const int nRho = (has_spatial ? DM_vec.size() : 0);
    const int nReinf = (has_reinfection ? X_rs.cols() : 0);
    const int nBeta = X.cols();
    const int nTrans = (exp_transition ? 2 : 0);
    total_size = nRho + nReinf + nBeta + nTrans;
    
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
    int time_idx, i, j, k;   
    unsigned int idx; 

    simulationResultSet compartmentResults;
 
    Eigen::VectorXd results = Eigen::VectorXd(1); 
    results(0) = 0.0;

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
   
    double gamma_ei = (exp_transition ? params(params.size() - 2) : -1.0);
    double gamma_ir = (exp_transition ? params(params.size() - 1) : -1.0);

    if (!exp_transition)
    {
        // Clear out paths
        for (i = 0; i < E_paths.cols(); i++)
        {
            for(j = 0; j < E_paths.rows(); j++)
            {
                E_paths(j,i) = 0;
            }
        }
        for (i = 0; i < I_paths.cols(); i++)
        {
            for(j = 0; j < I_paths.rows(); j++)
            {
                I_paths(j,i) = 0;
            }
        }
    }

    Eigen::VectorXi N(S0.size());
    N = (S0 + E0 + I0 + R0);

    Eigen::VectorXd p_se_cache(S0.size());

    Eigen::VectorXi current_S(S0.size());
    Eigen::VectorXi current_E(S0.size());
    Eigen::VectorXi current_I(S0.size());
    Eigen::VectorXi current_R(S0.size());

    Eigen::VectorXi previous_S(S0.size());
    Eigen::VectorXi previous_E(S0.size());
    Eigen::VectorXi previous_I(S0.size());
    Eigen::VectorXi previous_R(S0.size());


    Eigen::VectorXi previous_S_star(S0.size());
    Eigen::VectorXi previous_E_star(S0.size());
    Eigen::VectorXi previous_I_star(S0.size());
    Eigen::VectorXi previous_R_star(S0.size()); 

    Eigen::VectorXi cumulative_compartment(S0.size());

    Eigen::VectorXi* comparison_compartment = (data_compartment == 0 ?
                                               &previous_I_star : 
                                              (data_compartment == 1 ? 
                                               &previous_R_star : 
                                              (data_compartment == 2 ? 
                                               &previous_I : &previous_I_star)));

    // Calculate probabilities
    // p_se calculation
    // Equivalent R expression: 
    // exp(matrix(X %*% beta, nrow = nrow(Y), ncol = ncol(Y)))
    Eigen::MatrixXd eta = (X*beta).unaryExpr([](double elem){return(
                std::exp(elem));
            });
    Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor> p_se_components(eta.data(), 
                Y.rows(), Y.cols());

    time_idx = 0;
    previous_S = S0;
    previous_E = E0;
    previous_I = I0;
    previous_R = R0;

    if (!exp_transition)
    {
        I_paths.row(0) = I0;
        E_paths.row(0) = E0;
    }


    p_se_cache = (((p_se_components.row(0).array().transpose())
                   * (previous_I.cast<double>().array())) 
                   * (((N.cast<double>().array()).unaryExpr([](double e)
                    {
                        return(1.0/e);
                    })).array()));
    Eigen::VectorXd p_se = 1*p_se_cache; 

    if (has_spatial)
    {
        for (idx = 0; idx < DM_vec.size(); idx++)
        {
            p_se += rho[idx]*(DM_vec[idx] * (p_se_cache));
        }
    }

    p_se = (((-1.0*p_se.array()) * (offset(0)))).unaryExpr([](double e){
            return(1-std::exp(e));
            }); 




    // Not used if !exp_transition
    Eigen::VectorXd p_ei = (-1.0*gamma_ei*offset)
                            .unaryExpr([](double e){return(1-std::exp(e));});
    
    Eigen::VectorXd p_ir = (-1.0*gamma_ir*offset)
                            .unaryExpr([](double e){return(1-std::exp(e));}); 

    Eigen::VectorXd p_rs;
    if (has_reinfection)
    {
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

        compartmentResults.rEA = Eigen::MatrixXd(Y.rows(), Y.cols());
        compartmentResults.r0t = Eigen::MatrixXd(Y.rows(), Y.cols());
        compartmentResults.effR0 = Eigen::MatrixXd(Y.rows(), Y.cols());
        compartmentResults.p_se = Eigen::MatrixXd(Y.rows(), Y.cols());
        compartmentResults.p_se.row(0) = p_se;
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

    int tmpDraw;
    for (i = 0; i < Y.cols(); i++)
    {
        S_star_gen.param(std::binomial_distribution<>::param_type(previous_R(i) ,p_rs(0)));
        E_star_gen.param(std::binomial_distribution<>::param_type(previous_S(i) ,p_se(i)));

        previous_S_star(i) = S_star_gen(*generator);
        previous_E_star(i) = E_star_gen(*generator);

        if (exp_transition)
        {
            I_star_gen.param(std::binomial_distribution<>::param_type(previous_E(i) ,p_ei(0)));
            R_star_gen.param(std::binomial_distribution<>::param_type(previous_I(i) ,p_ir(0)));
            previous_I_star(i) = I_star_gen(*generator);
            previous_R_star(i) = R_star_gen(*generator);
        }
        else
        {
            // Updating E_paths and I_paths is repetative - factor out into a function?
            // That might be costly, unless it's inlined...
            previous_I_star(i) = E_paths(E_paths.rows() - 1, i); // could take outside loop
            E_paths(E_paths.rows() -1, i) = 0;
            for (j = 0; j < offset(0); j++)
            {
                // TODO: stop early when possible
                // idea: cache previous max?
                for (k = E_paths.rows() - 2; 
                        k >= 0; k--)
                {
                    if (E_paths(k,i) > 0)
                    {
                        I_star_gen.param(std::binomial_distribution<>::param_type(
                                    E_paths(k,i),E_to_I_prior(k, 5)));
                        tmpDraw = I_star_gen(*generator);
                        previous_I_star(i) += tmpDraw;
                        E_paths(k,i) -= tmpDraw;
                        E_paths(k+1, i) = E_paths(k,i);
                        E_paths(k,i) = 0; // not needed?
                    }
                }
            }
            previous_R_star(i) = I_paths(I_paths.rows() - 1, i); // could take outside loop
            I_paths(I_paths.rows() -1, i) = 0;
            for (j = 0; j < offset(0); j++)
            {
                // TODO: stop early when possible
                // idea: cache previous max?
                for (k = I_paths.rows() - 2; 
                        k >= 0; k--)
                {
                    if (I_paths(k,i) > 0)
                    {
                        R_star_gen.param(std::binomial_distribution<>::param_type(
                                    I_paths(k,i),I_to_R_prior(k, 5)));
                        tmpDraw = R_star_gen(*generator);
                        previous_R_star(i) += tmpDraw;
                        I_paths(k,i) -= tmpDraw;
                        I_paths(k+1, i) = I_paths(k,i);
                        I_paths(k,i) = 0; // not needed?
                    }
                }
            }
        }

        if (cumulative)
        {
            cumulative_compartment(i) = (*comparison_compartment)(i);
        }

        results(0) += (na_mask(0,i) ? 0 : 
                pow(((*comparison_compartment)(i) + 
                        (phi > 0 ? 
                         std::floor(overdispersion_distribution(*generator)) : 0)
                        - Y(0, i)), 2.0)); 
    }

    if (keepCompartments)
    {
        for (i = 0; i < S0.size(); i++)
        {
            compartmentResults.S(time_idx, i) = S0(i);
            compartmentResults.E(time_idx, i) = E0(i);
            compartmentResults.I(time_idx, i) = I0(i);
            compartmentResults.R(time_idx, i) = R0(i);

            compartmentResults.S_star(time_idx, i) = previous_S_star(i);
            compartmentResults.E_star(time_idx, i) = previous_E_star(i);
            compartmentResults.I_star(time_idx, i) = previous_I_star(i);
            compartmentResults.R_star(time_idx, i) = previous_R_star(i);
        }
    }

    current_S = previous_S + previous_S_star - previous_E_star;
    current_E = previous_E + previous_E_star - previous_I_star;
    current_I = previous_I + previous_I_star - previous_R_star;
    current_R = previous_R + previous_R_star - previous_S_star;

    if (!exp_transition)
    {
        E_paths.row(0) = previous_E_star;
        I_paths.row(0) = previous_I_star;
    }

    previous_S = current_S;
    previous_E = current_E;
    previous_I = current_I;
    previous_R = current_R;

    // Simulation: iterative case

    for (time_idx = 1; time_idx < Y.rows(); time_idx++)
    {
        p_se_cache = (((p_se_components.row(time_idx).array().transpose())
                   * (previous_I.cast<double>().array())) 
                   * (((N.cast<double>().array()).unaryExpr([](double e)
                    {
                        return(1.0/e);
                    })).array()));
        
        p_se = 1*p_se_cache; 
        if (has_spatial)
        {
            for (idx = 0; idx < DM_vec.size(); idx++)
            {
                p_se += rho[idx]*(DM_vec[idx] * p_se_cache);
            }
        }

        p_se = ((-1.0*p_se.array() * offset(time_idx)).matrix()).unaryExpr([](double e){return(1-std::exp(e));});
 
        for (i = 0; i < Y.cols(); i++)
        {
            S_star_gen.param(std::binomial_distribution<>::param_type(previous_R(i) ,p_rs(time_idx)));
            E_star_gen.param(std::binomial_distribution<>::param_type(previous_S(i) ,p_se(i)));
            I_star_gen.param(std::binomial_distribution<>::param_type(previous_E(i) ,p_ei(time_idx)));
            R_star_gen.param(std::binomial_distribution<>::param_type(previous_I(i) ,p_ir(time_idx)));

            previous_S_star(i) = S_star_gen(*generator);
            previous_E_star(i) = E_star_gen(*generator);

            if (exp_transition)
            {
                I_star_gen.param(std::binomial_distribution<>::param_type(previous_E(i) ,p_ei(time_idx)));
                R_star_gen.param(std::binomial_distribution<>::param_type(previous_I(i) ,p_ir(time_idx)));
                previous_I_star(i) = I_star_gen(*generator);
                previous_R_star(i) = R_star_gen(*generator);
            }
            else
            {
                // Updating E_paths and I_paths is repetative - factor out into a function?
                // That might be costly, unless it's inlined...
                previous_I_star(i) = E_paths(E_paths.rows() - 1, i); // could take outside loop
                E_paths(E_paths.rows() -1, i) = 0;
                for (j = 0; j < offset(time_idx); j++)
                {
                    // TODO: stop early when possible
                    // idea: cache previous max?
                    for (k = E_paths.rows() - 2; 
                            k >= 0; k--)
                    {
                        if (E_paths(k,i) > 0)
                        {
                            I_star_gen.param(std::binomial_distribution<>::param_type(
                                        E_paths(k,i),E_to_I_prior(k, 5)));
                            tmpDraw = I_star_gen(*generator);
                            previous_I_star(i) += tmpDraw;
                            E_paths(k,i) -= tmpDraw;
                            E_paths(k+1, i) = E_paths(k,i);
                            E_paths(k,i) = 0; // not needed?
                        }
                    }
                }
                previous_R_star(i) = I_paths(I_paths.rows() - 1, i); // could take outside loop
                I_paths(I_paths.rows() -1, i) = 0;
                for (j = 0; j < offset(time_idx); j++)
                {
                    // TODO: stop early when possible
                    // idea: cache previous max?
                    for (k = I_paths.rows() - 2; 
                            k >= 0; k--)
                    {
                        if (I_paths(k,i) > 0)
                        {
                            R_star_gen.param(std::binomial_distribution<>::param_type(
                                        I_paths(k,i),I_to_R_prior(k, 5)));
                            tmpDraw = R_star_gen(*generator);
                            previous_R_star(i) += tmpDraw;
                            I_paths(k,i) -= tmpDraw;
                            I_paths(k+1, i) = I_paths(k,i);
                            I_paths(k,i) = 0; // not needed?
                        }
                    }
                }
            }
            if (cumulative)
            {
                cumulative_compartment(i) += (*comparison_compartment)(i); 

                results(0) += (na_mask(time_idx, i) ? 0 : 
                    pow(((cumulative_compartment)(i) 
                            + (phi > 0 ? 
                                std::floor(overdispersion_distribution(*generator)) 
                                : 0)
                            - Y(time_idx, i)), 2.0)); 

            }
            else
            {
                results(0) += (na_mask(time_idx, i) ? 0 : 
                        pow(((*comparison_compartment)(i) + (phi > 0 ? 
                                    std::floor(overdispersion_distribution(*generator))
                                    : 0)- Y(time_idx, i)), 2.0)); 
            }
        }

        if (keepCompartments)
        {
            for (i = 0; i < S0.size(); i++)
            {
                compartmentResults.S(time_idx, i) = previous_S(i);
                compartmentResults.E(time_idx, i) = previous_E(i);
                compartmentResults.I(time_idx, i) = previous_I(i);
                compartmentResults.R(time_idx, i) = previous_R(i);

                compartmentResults.S_star(time_idx, i) = previous_S_star(i);
                compartmentResults.E_star(time_idx, i) = previous_E_star(i);
                compartmentResults.I_star(time_idx, i) = previous_I_star(i);
                compartmentResults.R_star(time_idx, i) = previous_R_star(i);

                compartmentResults.p_se.row(time_idx) = p_se;
            }
        }

        current_S = previous_S + previous_S_star - previous_E_star;
        current_E = previous_E + previous_E_star - previous_I_star;
        current_I = previous_I + previous_I_star - previous_R_star;
        current_R = previous_R + previous_R_star - previous_S_star;

        if (!exp_transition)
        {
            E_paths.row(0) = previous_E_star;
            I_paths.row(0) = previous_I_star;
        }

        previous_S = current_S;
        previous_E = current_E;
        previous_I = current_I;
        previous_R = current_R;
    }


    double resultVal = std::sqrt(results(0)); 
    if (keepCompartments)
    {
        calculateReproductiveNumbers(&compartmentResults);
    }
    compartmentResults.result = resultVal; 
    return(compartmentResults);
}

void SEIR_sim_node::calculateReproductiveNumbers(simulationResultSet* results)
{
    int i, l, k;
    Eigen::VectorXi N(S0.size());
    N = S0 + E0 + I0 + R0;

    int time_idx;
    Eigen::VectorXd beta = (*results).beta; 

    Eigen::MatrixXd eta = (((*results).X*beta).unaryExpr([](double elem){return(std::exp(elem));}));
    Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor> p_se_components(eta.data(), 
                (*results).S.rows(), (*results).S.cols());
    Eigen::MatrixXd p_se_cache(1, (*results).S.cols());

    int nLoc = (*results).S.cols();
    int nTpt = (*results).S.rows();
    std::vector<Eigen::MatrixXd> GVector;
    double component1, component2;
    // Fill in R0(t) and eff R0(t)
    // I_to_R_prior(k, 4)
    for (time_idx = 0; time_idx < nTpt; time_idx++)
    {
        (*results).r0t.row(time_idx) = p_se_components.row(time_idx);
        for (l = 0; l < p_se_components.cols(); l++)
        {
            if (exp_transition)
            {
                (*results).r0t(time_idx, l) /= (
                        -std::log(1.0-(*results).p_ir(time_idx))/offset(time_idx));
                (*results).effR0(time_idx, l) = (*results).r0t(time_idx, l)
                                            *((*results).S(time_idx, l))/N(l);
            }
            else
            {
                (*results).r0t(time_idx, l) /= avg_hazard;
                (*results).effR0(time_idx, l) = (*results).r0t(time_idx, l)
                                            *((*results).S(time_idx, l))/N(l);
               
            }
        }
    }


    // Calculate next generation matrices (more vector logic could be used here)
    for (time_idx = 0; time_idx < nTpt; time_idx++)
    {
        //Create and zero out G(t)
        Eigen::MatrixXd G(nLoc, nLoc);
        for (i = 0; i < nLoc; i++){for (l = 0; l < nLoc; l++){G(i,l) = 0.0;}}

        // Out: rows
        for (i = 0; i < nLoc; i++) 
        {
            // Out: columns
            for (l = 0; l < nLoc; l++)
            { 
                component1 = ((*results).I(time_idx, l) 
                        * (p_se_components(time_idx, l)))/N(l);
                if (i != l)
                {
                    component2 = 0.0;
                    for (k = 0; k < DM_vec.size(); k++)
                    {
                        component2 += ((*results).rho(k))*((DM_vec[k](i,l))*component1);
                    }
                    G(i,l) = ((*results).I(time_idx, l) != 0 ? 
                              (*results).S(time_idx, l)/((*results).I(time_idx, l))
                              *(1-std::exp(-component2)) : 0.0);
                }
                else
                { 
                    G(i,l) = ((*results).I(time_idx, l) != 0 ? 
                              (*results).S(time_idx, l)/((*results).I(time_idx, l))
                              * (1-std::exp(-component1)) : 0.0);
                }
            }
        } 
        GVector.push_back(G);
    }    


    for (time_idx = 0; time_idx < nTpt; time_idx ++)
   {
        for (i = 0; i < nLoc; i++)
        {
            (*results).rEA(time_idx, i) = 0.0;
        }

        for (i = 0; i < nLoc; i++)
        {
            for (l = 0; l < nLoc; l++)
            {
                (*results).rEA(time_idx, l) += GVector[time_idx](l, i);        
            }
        }
    }



    double _1mpIR_cum = 1.0;
    int infTime = 0;
    if (exp_transition)
    {
        for (time_idx = 0; time_idx < nTpt; time_idx++)
        {
            _1mpIR_cum = (1-(*results).p_ir(time_idx, 0));
            for (l = time_idx+1; l < nTpt; l++)
            {
                (*results).rEA.row(time_idx) += 
                    _1mpIR_cum*(*results).rEA.row(l);
                _1mpIR_cum *= (1 - (*results).p_ir(l, 0)); 
            }
            while (_1mpIR_cum > 1e-8)
            {
                (*results).rEA.row(nTpt - 1) += _1mpIR_cum*(*results).rEA.row(nTpt - 1); 
                _1mpIR_cum*=((*results).p_ir(nTpt - 1,0));
            }
        }
    }
    else
    {
        for (time_idx = 0; time_idx < nTpt; time_idx++)
        {
            _1mpIR_cum = 1.0;
            infTime = 0;
            for (i = 0; i < offset(time_idx); i++)
            {
                _1mpIR_cum *= (1 - I_to_R_prior(infTime, 5)); 
                infTime ++;
            }
            for (l = time_idx+1; (l < nTpt && infTime < I_to_R_prior.rows()); l++)
            {
                (*results).rEA.row(time_idx) += 
                    _1mpIR_cum*(*results).rEA.row(l);
                for (i = 0; i < offset(time_idx); i++)
                {
                    _1mpIR_cum *= (1 - I_to_R_prior(infTime, 5)); 
                    infTime ++;
                }
            }
            while (_1mpIR_cum > 1e-8 && infTime < I_to_R_prior.rows())
            {
                (*results).rEA.row(nTpt - 1) += _1mpIR_cum*(*results).rEA.row(nTpt - 1); 
                _1mpIR_cum *= (1 - I_to_R_prior(infTime, 5)); 
                infTime ++;
            }
        }
    }
}

SEIR_sim_node::~SEIR_sim_node()
{
    delete generator;
}

