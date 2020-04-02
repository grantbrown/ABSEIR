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


                                                                                
std::vector<size_t> sort_indexes(std::vector<int> inVec);
std::vector<size_t> sort_indexes_eigen(Eigen::MatrixXd inMat);
std::vector<size_t> sort_indexes_eigen_vec(Eigen::VectorXd inVec);

Rcpp::List spatialSEIRModel::sample_basic(int nSample, int vb, 
                                          std::string sim_type_atom)
{
    // This is set by constructor
    const int nParams = param_matrix.cols();

    // Accepted params/results are nSample size
    results_complete = std::vector<simulationResultSet>();
    proposed_results_complete = std::vector<simulationResultSet>();
    results_double = Eigen::MatrixXd::Zero(nSample, 
                                           samplingControlInstance -> m); 
    param_matrix = Eigen::MatrixXd::Zero(nSample, 
                                            nParams);
    prev_param_matrix = param_matrix;
    // Proposals matrices are batch size
    proposed_results_double = Eigen::MatrixXd::Zero(nSample,
                                               samplingControlInstance -> m); 
    proposed_param_matrix = Eigen::MatrixXd::Zero(nSample,
                                                  nParams);

    preproposal_params = Eigen::MatrixXd::Zero(samplingControlInstance -> init_batch_size, 
                                               nParams);
    preproposal_results = Eigen::MatrixXd::Zero(samplingControlInstance -> init_batch_size, 
                                                samplingControlInstance -> m); 
 




 
    const int verbose = vb;
    if (verbose > 1)
    {
        Rcpp::Rcout << "Starting sampler\n";
    }
    const bool hasReinfection = (reinfectionModelInstance -> 
               betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;

    std::string transitionMode = transitionPriorsInstance -> mode;   
    std::vector<size_t> reweight_idx;

    int i;

    if (verbose > 1)
    {
        dataModelInstance -> summary();
        exposureModelInstance -> summary();
        transitionPriorsInstance -> summary();
        reinfectionModelInstance -> summary();
        distanceModelInstance -> summary();
        initialValueContainerInstance -> summary();
        samplingControlInstance -> summary();
    }
    if (!is_initialized)
    {
        if (verbose > 1){Rcpp::Rcout << "Generating starting parameters from prior\n";}
        // Sample parameters from their prior

        preproposal_params = generateParamsPrior(samplingControlInstance -> init_batch_size);
        run_simulations(preproposal_params, sim_type_atom, &preproposal_results, &proposed_results_complete);

        std::vector<size_t> currentIndex = sort_indexes_eigen(preproposal_results); 
        std::vector<size_t> rcIdx = sort_indexes(result_idx); 
        for (i = 0; i < param_matrix.rows(); i++)
        {
            param_matrix.row(i) = preproposal_params.row(currentIndex[i]);
            results_double.row(i) = preproposal_results.row(currentIndex[i]); 
            if (sim_type_atom == sim_result_atom)
            {
                results_complete.push_back(proposed_results_complete[rcIdx[currentIndex[i]]]);
            } 
        }
    }
    Rcpp::List outList;

    // Todo: keep an eye on this object handling. It may have unreasonable
    // overhead, and is kind of complex.  
    if (sim_type_atom == sim_result_atom)
    {

        Rcpp::List simulationResults; 
        for (i = 0; i < (int) results_complete.size(); i++)
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
            simulationResults[std::to_string(i)] = subList;
        }
        outList["simulationResults"] = simulationResults;
    }
       
    outList["result"] = Rcpp::wrap(results_double);
    outList["params"] = Rcpp::wrap(param_matrix);
    outList["currentEps"] = results_double.maxCoeff() + 1.0;
;
    return(outList);
}
