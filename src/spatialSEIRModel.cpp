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


Rcpp::NumericVector spatialSEIRModel::marginalPosteriorEstimates(SEXP inParams)
{
    Rcpp::NumericMatrix params(inParams);
    self = new scoped_actor();
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
    
    std::vector<caf::actor> workers;

    auto worker_pool = actor_pool::make(actor_pool::round_robin{});
    unsigned int ncore = (unsigned int) samplingControlInstance -> CPU_cores;
    unsigned int nrow =  (unsigned int) params.nrow();
    
    /*
    auto tmp1 = samplingControlInstance -> simulation_width;
    auto tmp2 = samplingControlInstance->random_seed;
    auto tmp3 = initialValueContainerInstance -> S0;
    auto tmp4 = initialValueContainerInstance -> E0;
    auto tmp5 = initialValueContainerInstance -> I0;
    auto tmp6 = initialValueContainerInstance -> R0;
    auto tmp7 = exposureModelInstance -> offset;
    auto tmp8 = dataModelInstance -> Y;
    auto tmp9 = distanceModelInstance -> dm_list;
    auto tmp10 = exposureModelInstance -> X;
    auto tmp11 = reinfectionModelInstance -> X_rs;
    auto tmp12 = transitionPriorsInstance -> gamma_ei_params;
    auto tmp13 = transitionPriorsInstance -> gamma_ir_params;
    auto tmp14 = exposureModelInstance -> betaPriorPrecision;
    auto tmp15 = reinfectionModelInstance -> betaPriorPrecision; 
    auto tmp16 = exposureModelInstance -> betaPriorMean;
    auto tmp17 = reinfectionModelInstance -> betaPriorMean;
    auto tmp18 = dataModelInstance -> phi;
    */

    for (i = 0; i < ncore; i++)
    {
        workers.push_back((*self) -> spawn<SEIR_sim_node, monitored>(samplingControlInstance->simulation_width,
                                                                      samplingControlInstance->random_seed,
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

    for (i = 0; i < nrow; i++)
    {
        unsigned int outIdx = i;
        outRow = param_matrix.row(i);
        (*self) -> send(worker_pool, sim_atom::value, outIdx, outRow); 
    }

    std::vector<int> result_idx;
    std::vector<double> results;

    i = 0;
    (*self)->receive_for(i, nrow)(
                     [&](unsigned int idx, double result) {
                          results.push_back(result);
                          result_idx.push_back(idx);
                        });

    (*self) -> send_exit(worker_pool, exit_reason::user_shutdown);
    Rcpp::NumericVector out(nrow); 

    for (i = 0; i < nrow; i++)
    {
        out(result_idx[i]) = results[i];
    }
    delete self;
    shutdown();
    return(out);
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
    .method("marginalPosteriorEstimates", &spatialSEIRModel::marginalPosteriorEstimates);

}

