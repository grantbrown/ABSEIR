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

Rcpp::List spatialSEIRModel::sample_Simulate(int nSample, 
                                             int verbose) 
{

    double eps = init_eps;
    if (!is_initialized)
    {
        Rcpp::stop("Simulation requires initialized parameters");
    }

    const bool hasReinfection = (reinfectionModelInstance -> 
            betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;

    std::string transitionMode = transitionPriorsInstance -> mode;   

    int Naccept = 0;
    int batch;
    int i;
    samplingControlInstance -> m = 1;
    results_double = Eigen::MatrixXd::Zero(param_matrix.rows(), 1);
    results_complete = std::vector<simulationResultSet>();

    std::vector<simulationResultSet> finalResults = std::vector<simulationResultSet>();
    for (batch = 0; batch < samplingControlInstance -> max_batches &&
            Naccept < nSample; batch ++)
    {
        run_simulations(param_matrix, 
                        sim_result_atom,
                        &results_double,
                        &results_complete); 

        for (i = 0; i < results_double.rows(); i++)
        {
            if (results_double(i,0) < eps)
            {
                finalResults.push_back(results_complete[i]);
                Naccept++;
            }
            if (Naccept == nSample)
            {
                break;
            }
        }
    }


    Rcpp::List outList;

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

    outList["params"] = Rcpp::wrap(param_matrix);
    outList["currentEps"] = eps;
    return(outList);
}
