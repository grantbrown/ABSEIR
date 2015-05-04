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
        ::Rf_error("Exposure model and data model imply different number of locations\n");
    }
    if ((dataModelInstance -> nTpt) != (exposureModelInstance -> nTpt))
    { 
        ::Rf_error("Exposure model and data model imply different number of time points\n");  
    }
    if ((dataModelInstance -> nLoc) != (distanceModelInstance -> numLocations))
    {       
        ::Rf_error("Data model and distance model imply different number of locations\n");
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
    if (*(transitionPriorsInstance -> gamma_ei) < 0 || 
        *(transitionPriorsInstance -> gamma_ir) < 0)
    { 
        ::Rf_error("Transition priors haven't been populated (or have invalid values)\n");
    }
    if ((reinfectionModelInstance -> reinfectionMode) > 2)
    {
        // pass
    }

}


Rcpp::NumericVector spatialSEIRModel::marginalPosteriorEstimates(Rcpp::NumericMatrix params)
{
    self = new scoped_actor();

    std::vector<caf::actor> workers;
    size_t i;
    for (i = 0; i < samplingControlInstance -> CPU_cores; i++)
    {

        workers.push_back((*self) -> spawn<SEIR_sim_node, monitored>(samplingControlInstance->simulation_width,
                                                                      samplingControlInstance->random_seed,
                                                                      initialValueContainerInstance -> S0,
                                                                      initialValueContainerInstance -> E0,
                                                                      initialValueContainerInstance -> I0,
                                                                      initialValueContainerInstance -> R0,
                                                                      dataModelInstance -> Y,
                                                                      distanceModelInstance -> dm_list,
                                                                      exposureModelInstance -> X,
                                                                      reinfectionModelInstance -> X_rs,
                                                                      (*self)));
    }
    // Dummy output
    Rcpp::NumericVector out(0);

    delete self;
    return(out);
}

spatialSEIRModel::~spatialSEIRModel()
{   
    dataModelInstance -> unprotect();
    exposureModelInstance -> unprotect();
    reinfectionModelInstance -> unprotect();
    distanceModelInstance -> unprotect();
    transitionPriorsInstance -> unprotect();
    initialValueContainerInstance -> unprotect();
    samplingControlInstance -> unprotect();

    delete self;
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
    .method("runGridSearch", &spatialSEIRModel::marginalPosteriorEstimates);

}

