#include <Rcpp.h>
#include <cmath>
#include <spatialSEIRModel.hpp>
#include <dataModel.hpp>
#include <exposureModel.hpp>
#include <reinfectionModel.hpp>
#include <distanceModel.hpp>
#include <transitionPriors.hpp>
#include <initialValueContainer.hpp>
#include <samplingControl.hpp>
#include <SEIRSimNodes.hpp>
#include "caf/all.hpp"

using namespace Rcpp;
using namespace caf;

Rcpp::IntegerMatrix createRcppIntFromEigen(Eigen::MatrixXi inMatrix)
{
    Rcpp::IntegerMatrix outMatrix(inMatrix.rows(), inMatrix.cols());
    int i, j;
    for (i = 0; i < inMatrix.cols(); i++)
    {
        for (j = 0; j < inMatrix.rows(); j++)
        {
            outMatrix(j,i) = inMatrix(j,i);
        }
    }
    return(outMatrix);
}

Rcpp::NumericMatrix createRcppNumericFromEigen(Eigen::MatrixXd inMatrix)
{
    Rcpp::NumericMatrix outMatrix(inMatrix.rows(), inMatrix.cols());
    int i, j;
    for (i = 0; i < inMatrix.cols(); i++)
    {
        for (j = 0; j < inMatrix.rows(); j++)
        {
            outMatrix(j,i) = inMatrix(j,i);
        }
    }
    return(outMatrix);
}



spatialSEIRModel::spatialSEIRModel(dataModel& dataModel_,
                                   exposureModel& exposureModel_,
                                   reinfectionModel& reinfectionModel_,
                                   distanceModel& distanceModel_,
                                   transitionPriors& transitionPriors_,
                                   initialValueContainer& initialValueContainer_,
                                   samplingControl& samplingControl_)
{
    // Make sure these pointers go to the real deal
    int err = (((dataModel_.getModelComponentType()) != LSS_DATA_MODEL_TYPE) ||
            ((exposureModel_.getModelComponentType()) != LSS_EXPOSURE_MODEL_TYPE) ||
            ((reinfectionModel_.getModelComponentType()) != LSS_REINFECTION_MODEL_TYPE) ||
            ((distanceModel_.getModelComponentType()) != LSS_DISTANCE_MODEL_TYPE) ||
            ((transitionPriors_.getModelComponentType()) != LSS_TRANSITION_MODEL_TYPE) ||
            ((initialValueContainer_.getModelComponentType()) != LSS_INIT_CONTAINER_TYPE) ||
            ((samplingControl_.getModelComponentType()) != LSS_SAMPLING_CONTROL_MODEL_TYPE));

    if (err != 0)
    { 
        ::Rf_error("Error: model components were not provided in the correct order. \n");
    }

    ncalls = 0;


    int i;
    dataModelInstance = &dataModel_;
    exposureModelInstance = &exposureModel_;
    reinfectionModelInstance = &reinfectionModel_;
    distanceModelInstance = &distanceModel_;
    transitionPriorsInstance = &transitionPriors_;
    initialValueContainerInstance = &initialValueContainer_;
    samplingControlInstance = &samplingControl_;

    dataModelInstance -> protect();
    exposureModelInstance -> protect();
    reinfectionModelInstance -> protect();
    distanceModelInstance -> protect();
    transitionPriorsInstance -> protect();
    initialValueContainerInstance -> protect();
    samplingControlInstance -> protect();

    if ((dataModelInstance -> nLoc) != (exposureModelInstance -> nLoc))
    { 
        ::Rf_error(("Exposure model and data model imply different number of locations: " 
                + std::to_string(dataModelInstance -> nLoc) + ", " 
                + std::to_string(exposureModelInstance -> nLoc) + ".\n").c_str());
    }
    if ((dataModelInstance -> nTpt) != (exposureModelInstance -> nTpt))
    { 
        ::Rf_error(("Exposure model and data model imply different number of time points:"
                    + std::to_string(dataModelInstance -> nTpt) + ", "
                    + std::to_string(exposureModelInstance -> nTpt) + ".\n").c_str());  
    }
    if ((dataModelInstance -> nLoc) != (distanceModelInstance -> numLocations))
    {       
        ::Rf_error(("Data model and distance model imply different number of locations:"
                    + std::to_string(dataModelInstance -> nLoc) + ", "
                    + std::to_string(distanceModelInstance -> numLocations) + ".\n").c_str()
                );
    }
    if ((dataModelInstance -> nLoc) != (initialValueContainerInstance -> S0.size())) 
    { 
        ::Rf_error("Data model and initial value container have different dimensions\n");
    }
    if ((reinfectionModelInstance -> reinfectionMode) == 3)
    {
        // No reinfection
    }
    else
    {
        if (((reinfectionModelInstance -> X_rs).rows()) != (dataModelInstance -> nTpt))
        { 
            ::Rf_error("Reinfection and data mode time points differ.\n");
        }
    }
    if ((reinfectionModelInstance -> reinfectionMode) > 2)
    {
        // pass
    }

}

Rcpp::List spatialSEIRModel::sample(SEXP nSamples, SEXP rejectionFraction)
{
    // yada, create param_matrix
    Eigen::MatrixXd param_matrix(10, 10);
    Rcpp::List outList = this -> simulate(param_matrix, sample_atom::value);
    return(outList);
}

Rcpp::List spatialSEIRModel::evaluate(SEXP inParams)
{
    Rcpp::NumericMatrix params(inParams);
    // Copy to Eigen matrix
    unsigned int i, j;
    Eigen::MatrixXd param_matrix(params.nrow(), params.ncol());
    for (i = 0; i < params.nrow(); i++)
    {
        for (j = 0; j < params.ncol(); j++)
        {
            param_matrix(i,j) = params(i,j);
        }
    }
    return(this -> simulate(param_matrix, sim_atom::value));
}

Rcpp::List spatialSEIRModel::simulate_given(SEXP inParams)
{
    Rcpp::NumericMatrix params(inParams);
    // Copy to Eigen matrix
    unsigned int i, j;
    Eigen::MatrixXd param_matrix(params.nrow(), params.ncol());
    for (i = 0; i < params.nrow(); i++)
    {
        for (j = 0; j < params.ncol(); j++)
        {
            param_matrix(i,j) = params(i,j);
        }
    }
    return(this -> simulate(param_matrix, sim_result_atom::value));
}

Rcpp::List spatialSEIRModel::simulate(Eigen::MatrixXd param_matrix, 
                                      caf::atom_value sim_type_atom)
{
    if (!(sim_type_atom == sim_atom::value || 
         sim_type_atom == sim_result_atom::value || 
         sim_type_atom == sample_atom::value))
    {
        Rcpp::Rcout << "Invalid simulation type requested\n";
        Rcpp::List outList;
        outList["error"] = true;
        return(outList);
    }

    self = new scoped_actor();
    ncalls += 1;    
    std::vector<caf::actor> workers;

    auto worker_pool = actor_pool::make(actor_pool::round_robin{});
    unsigned int ncore = (unsigned int) samplingControlInstance -> CPU_cores;
    unsigned int nrow =  (unsigned int) param_matrix.rows(); 
    unsigned int i;

    for (i = 0; i < ncore; i++)
    {
        workers.push_back((*self) -> spawn<SEIR_sim_node, monitored>(samplingControlInstance->simulation_width,
                                                                      samplingControlInstance->random_seed + 1000*i + ncalls,
                                                                      initialValueContainerInstance -> S0,
                                                                      initialValueContainerInstance -> E0,
                                                                      initialValueContainerInstance -> I0,
                                                                      initialValueContainerInstance -> R0,
                                                                      exposureModelInstance -> offset,
                                                                      dataModelInstance -> Y,
                                                                      distanceModelInstance -> dm_list,
                                                                      exposureModelInstance -> X,
                                                                      reinfectionModelInstance -> X_rs,
                                                                      transitionPriorsInstance -> gamma_ei_params,
                                                                      transitionPriorsInstance -> gamma_ir_params,
                                                                      distanceModelInstance -> spatial_prior,
                                                                      exposureModelInstance -> betaPriorPrecision,
                                                                      reinfectionModelInstance -> betaPriorPrecision, 
                                                                      exposureModelInstance -> betaPriorMean,
                                                                      reinfectionModelInstance -> betaPriorMean,
                                                                      dataModelInstance -> phi,
                                                                      (*self)));
        (*self) -> send(worker_pool, sys_atom::value, put_atom::value, workers[workers.size()-1]);
    }

    // Distribute jobs among workers
    unsigned int outIdx;
    Eigen::VectorXd outRow;

    // Send simulation orders
    for (i = 0; i < nrow; i++)
    {
        outRow = param_matrix.row(i);
        (*self) -> send(worker_pool, sim_type_atom, i, outRow); 
    }

    std::vector<int> result_idx;
    // Having two results containers is kinda ugly, better solution?
    std::vector<simulationResultSet> results_complete;
    std::vector<double> results_double;


    i = 0;
    if (sim_type_atom == sim_result_atom::value)
    {
        (*self)->receive_for(i, nrow)(
                         [&](unsigned int idx, simulationResultSet result) {
                              results_complete.push_back(result);
                              result_idx.push_back(idx);
                            });
    }
    else if (sim_type_atom == sim_atom::value || sim_type_atom == sample_atom::value)
    {
        (*self)->receive_for(i, nrow)(
                     [&](unsigned int idx, double result) {
                          results_double.push_back(result);
                          result_idx.push_back(idx);
                        });
    }

    (*self) -> send_exit(worker_pool, exit_reason::user_shutdown); 
    
    // Todo: keep an eye on this object handling. It may have unreasonable
    // overhead, and is kind of complex.  
    Rcpp::List outList;

    if (sim_type_atom == sim_result_atom::value)
    {
        // keep_samples indicates a debug mode, so don't worry if we can't make
        // a regular data frame from the list.
        for (i = 0; i < nrow; i++)
        {
            Rcpp::List subList;
            subList["S"] = createRcppIntFromEigen(results_complete[i].S);
            subList["E"] = createRcppIntFromEigen(results_complete[i].E);
            subList["I"] = createRcppIntFromEigen(results_complete[i].I);
            subList["R"] = createRcppIntFromEigen(results_complete[i].R);

            subList["S_star"] = createRcppIntFromEigen(results_complete[i].S_star);
            subList["E_star"] = createRcppIntFromEigen(results_complete[i].E_star);
            subList["I_star"] = createRcppIntFromEigen(results_complete[i].I_star);
            subList["R_star"] = createRcppIntFromEigen(results_complete[i].R_star);
            subList["p_se"] = createRcppNumericFromEigen(results_complete[i].p_se);
            subList["p_ei"] = createRcppNumericFromEigen(results_complete[i].p_ei);
            subList["p_ir"] = createRcppNumericFromEigen(results_complete[i].p_ir);
            subList["rho"] = createRcppNumericFromEigen(results_complete[i].rho);
            subList["beta"] = createRcppNumericFromEigen(results_complete[i].beta);
            subList["X"] = createRcppNumericFromEigen(results_complete[i].X);
            subList["result"] = results_complete[i].result;
            outList[std::to_string(i)] = subList;
        }
    }
    else if (sim_type_atom == sim_atom::value)
    {
        Rcpp::NumericVector outResults(nrow);
        for (i = 0; i < nrow; i++)
        {
            outResults(result_idx[i]) = results_double[i];
        }
        outList["result"] = outResults;
    }
    else if (sim_type_atom == sample_atom::value)
    {
        Rcpp::NumericVector outResults(nrow);
        for (i = 0; i < nrow; i++)
        {
            outResults(result_idx[i]) = results_double[i];
        }
        // Todo: add param_matrix results here. 
        outList["result"] = outResults;
    }
    delete self;
    return(outList);
}

spatialSEIRModel::~spatialSEIRModel()
{   
    shutdown();
    dataModelInstance -> unprotect();
    exposureModelInstance -> unprotect();
    reinfectionModelInstance -> unprotect();
    distanceModelInstance -> unprotect();
    transitionPriorsInstance -> unprotect();
    initialValueContainerInstance -> unprotect();
    samplingControlInstance -> unprotect();
}


RCPP_MODULE(mod_spatialSEIRModel)
{
    using namespace Rcpp;
    class_<spatialSEIRModel>( "spatialSEIRModel" )
    .constructor<dataModel&,
                 exposureModel&,
                 reinfectionModel&,
                 distanceModel&,
                 transitionPriors&,
                 initialValueContainer&,
                 samplingControl&>()
    .method("evaluate", &spatialSEIRModel::evaluate)
    .method("sample", &spatialSEIRModel::sample)
    .method("simulate", &spatialSEIRModel::simulate_given);

}

