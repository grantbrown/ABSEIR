#include <Rcpp.h>
#include <Eigen/Core>
#include <RcppEigen.h>
#include <cmath>
#include <math.h>
#include <spatialSEIRModel.hpp>
#include <dataModel.hpp>
#include <exposureModel.hpp>
#include <reinfectionModel.hpp>
#include <distanceModel.hpp>
#include <transitionPriors.hpp>
#include <initialValueContainer.hpp>
#include <samplingControl.hpp>
#include <util.hpp>
#include <SEIRSimNodes.hpp>

Eigen::VectorXd calculate_weights(double cur_e,
                                  double prev_e, 
                                  Eigen::MatrixXd eps,
                                  Eigen::VectorXd prev_wts)
{
    Eigen::VectorXd out_wts = Eigen::VectorXd::Zero(eps.rows());
    int i, j;
    int M = eps.cols();
    double num, denom, tot;
    for (i = 0; i < eps.rows(); i++)
    {
        num = 0.0;
        denom = 0.0;
        for (j = 0; j < M; j++)
        {
            num += eps(i,j) < cur_e;
            denom += eps(i,j) < prev_e;
        } 
        out_wts(i) = (num/denom*prev_wts(i));
        tot += out_wts(i);
    }
    out_wts.array() /= tot;
    return(out_wts);
}

double ESS(Eigen::VectorXd wts)
{
    return(1.0/wts.squaredNorm()); 
}

double eps_f(double rhs,
             double cur_e,
             double prev_e, 
             Eigen::MatrixXd eps,
             Eigen::VectorXd prev_wts)
{
    return(std::pow(rhs - ESS(calculate_weights(cur_e, 
                                                prev_e, 
                                                eps, 
                                                prev_wts)), 2.0));
}

double solve_for_epsilon(double LB,
                       double UB,
                       double prev_e,
                       double alpha,
                       Eigen::MatrixXd eps,
                       Eigen::VectorXd prev_wts)
{

    double phi = (1.0 + std::sqrt(5))/2.0; 
    double a,b,c,d,fa,fb,fc,fd;
    double proposed_e = proposed_e; 
    double rhs = ESS(calculate_weights(proposed_e,
                                         prev_e, 
                                         eps,
                                         prev_wts))*alpha;

    a = LB;
    b = UB;
    fa = eps_f(rhs, a, prev_e, eps, prev_wts);
    fb = eps_f(rhs, b, prev_e, eps, prev_wts);
    bool mvLB = true;
    bool mvUB = true;
    int itrs = 0;
    double diff = 100.0;
    while (itrs < 10000 && diff > 1e-8)
    {
       if (mvLB)
       {
           c = b + (a - b)/phi; 
           mvLB = false;
       }
       if (mvUB)
       {
           d = a + (b - a)/phi;
           mvUB = false;
       }
       fc = eps_f(rhs, c, prev_e, eps, prev_wts);
       fd = eps_f(rhs, d, prev_e, eps, prev_wts);

       if (fc < fd)
       {
         b = d;
         fb = fd;
         d = c;
         fd = fc;
         mvLB = true;
       }
       else
       {
         a = c;
         fa = fc;
         c = d;
         fc = fd;
         mvUB = true;
       }
       diff = b-a;
    }
    if (diff > 1e-8)
    {
       Rcpp::Rcout << "Warning: optimization didn't converge.\n";
    }
    return(a);

}

Rcpp::List spatialSEIRModel::sample_DelMoral2012(int nSample, bool verbose, 
                                                 std::string sim_type_atom)
{
    int num_iterations = samplingControlInstance -> epochs;
    int N = nSample;
    const bool hasReinfection = (reinfectionModelInstance -> 
            betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;

    std::string transitionMode = transitionPriorsInstance -> mode;   

    double current_eps = std::numeric_limits<double>::infinity();
    double next_eps;

    int i,j;
    int iteration;

    auto U = std::uniform_real_distribution<double>(0,1);
    double drw;
    
    if (!is_initialized)
    {
        // Sample parameters from their prior
        Rcpp::Rcout << "N: " << N << "\n";
        param_matrix = generateParamsPrior(N);
        Rcpp::Rcout << "\n\n Generated!!\n\n";
        Rcpp::Rcout << "param_matrix: " << param_matrix.rows() << ", " << param_matrix.cols() << "\n";
        run_simulations(param_matrix, sim_atom, &results_double, &results_complete);
    }
    else
    {
        // The data in "param_matrix" is already accepted
        // To-do: finish this clause
    }
    // Step 0b: set weights to 1/N
    Eigen::VectorXd weights = Eigen::VectorXd::Zero(N).array() + 1.0/((double) N);
    Eigen::VectorXd next_weights = weights;
    Eigen::VectorXd cum_weights = next_weights;

    for (iteration = 0; iteration < num_iterations; iteration++)
    {   
        next_eps = solve_for_epsilon(results_double.minCoeff() + 1.0,
                                     results_double.maxCoeff(),
                                     current_eps,
                                     samplingControlInstance -> shrinkage,
                                     results_double,
                                     weights);
        next_weights = calculate_weights(next_eps,
                                         current_eps, 
                                         results_double,
                                         weights);
        if (ESS(next_weights) < N)
        {
           // Resample N particles 

           prev_results_double = results_double;
           prev_param_matrix = param_matrix;

           // Compute cumulative weights
           cum_weights(0) = weights(0);
           for (i = 1; i < weights.size(); i++)
           {
               cum_weights(i) = weights(i) + cum_weights(i-1);
           }
           // Fill in param_matrix with resamples
           for (i = 0; i < N; i++)
           {
               drw = U(*generator);
               for (j = 0; j < cum_weights.size(); j++)
               {
                    if (drw <= cum_weights(j))
                    {
                        param_matrix.row(i) = prev_param_matrix.row(j); 
                        results_double.row(i) = prev_results_double.row(j);
                    }
               }
           }
           // Reset weights
           for (i = 0; i < weights.size(); i++)
           {
               weights(i) = 1.0/((double) N);
           } 
        }
        else
        {
           prev_results_double = results_double;
           prev_param_matrix = param_matrix;
        }

        // Step 3: MCMC update
        Eigen::MatrixXd proposed_params = param_matrix;
        Eigen::VectorXd tau = ((param_matrix.rowwise() - 
                                param_matrix.colwise().mean()
                                ).colwise().norm()/std::sqrt((double) 
                                    param_matrix.rows()-1.0));
        // Propose new parameters
        for (j = 0; j < proposed_params.cols(); j++)
        {
            auto propDist = std::normal_distribution<double>(0.0, tau(j));
            for (i = 0; i < N; i++)
            {
                proposed_params(i,j) += propDist(*generator);
            }
        }
        run_simulations(proposed_params, sim_atom, &results_double, &results_complete);

        double acc_ratio, num, denom, pn, pd;
        bool accept;
        for (i = 0; i < N; i++)
        {
            pn = evalPrior(proposed_params.row(i));
            pd = evalPrior(param_matrix.row(i));
            num = 0.0;
            denom = 0.0;
            accept = false;
            if (std::isfinite(num) && num > 0.0)
            {
                for (j = 0; j < results_double.cols(); j++)
                {
                    num += (results_double(i,j) < current_eps);
                    denom += (prev_results_double(i,j) < current_eps);
                }
                num *= pn;
                denom *= pd;
                acc_ratio = num/denom;
                drw = U(*generator);
                if (!std::isfinite(acc_ratio))
                {
                    Rcpp::Rcout << "Potential problem: non-finite acceptance\n";
                }
                else if (drw < acc_ratio)
                {
                    param_matrix.row(i) = proposed_params.row(i);
                }
            }
        } 
    }

    // Todo: keep an eye on this object handling. It may have unreasonable
    // overhead, and is kind of complex.  
    Rcpp::List outList;
    if (sim_type_atom == sim_result_atom)
    {
        // keep_samples indicates a debug mode, so don't worry if we can't make
        // a regular data frame from the list.
        for (i = 0; i < results_complete.size(); i++)
        {
            Rcpp::List subList;
            subList["S"] = Rcpp::wrap(results_complete[i].S);
            subList["E"] = Rcpp::wrap(results_complete[i].E);
            subList["I"] = Rcpp::wrap(results_complete[i].I);
            subList["R"] = Rcpp::wrap(results_complete[i].R);

            subList["S_star"] = Rcpp::wrap(results_complete[i].S_star);
            subList["E_star"] = Rcpp::wrap(results_complete[i].E_star);
            subList["I_star"] = Rcpp::wrap(results_complete[i].I_star);
            subList["R_star"] = Rcpp::wrap(results_complete[i].R_star);
            subList["p_se"] = Rcpp::wrap(results_complete[i].p_se);
            // We p_ei and p_ir not generally defined in non-exponential case.  
            if (transitionMode == "exponential")
            {
                subList["p_ei"] = Rcpp::wrap(results_complete[i].p_ei);
                subList["p_ir"] = Rcpp::wrap(results_complete[i].p_ir);
            }
            subList["R_EA"] = Rcpp::wrap(results_complete[i].rEA);
            subList["R0t"] = Rcpp::wrap(results_complete[i].r0t);
            subList["effR0t"] = Rcpp::wrap(results_complete[i].effR0);
            if (hasSpatial)
            {
                subList["rho"] = Rcpp::wrap(results_complete[i].rho);
            }
            subList["beta"] = Rcpp::wrap(results_complete[i].beta);
            subList["X"] = Rcpp::wrap(results_complete[i].X);
            if (hasReinfection)
            {
                // TODO: output reinfection info
            }
            subList["result"] = Rcpp::wrap(results_complete[i].result);
            outList[std::to_string(i)] = subList;
        }
    }
    else if (sim_type_atom == sim_atom)
    {
        outList["result"] = Rcpp::wrap(results_double);
    }
    outList["params"] = Rcpp::wrap(param_matrix);
    outList["currentEps"] = current_eps;
    return(outList);
}
