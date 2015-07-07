#ifndef SPATIALSEIR_MODEL_INC
#define SPATIALSEIR_MODEL_INC
#include <Rcpp.h>
#include "caf/all.hpp"
#include "./dataModel.hpp"
#include "./distanceModel.hpp"
#include "./exposureModel.hpp"
#include "./initialValueContainer.hpp"
#include "./reinfectionModel.hpp"
#include "./samplingControl.hpp"
#include "./spatialSEIRModel.hpp"
#include "./transitionPriors.hpp"


class dataModel;
class exposureModel;
class distanceModel;
class initialValueContainer;
class reinfectionModel;
class samplingControl;
class transitionPriors;

using namespace Rcpp;
using namespace caf;

class spatialSEIRModel
{
    public: 
        //Constructor
        spatialSEIRModel(dataModel& dataModel_,
                         exposureModel& exposureModel_,
                         reinfectionModel& reinfectionModel_,
                         distanceModel& distanceModel_,
                         transitionPriors& transitionPriors_,
                         initialValueContainer& initialValueContainer_,
                         samplingControl& samplingControl_);

        Rcpp::NumericVector marginalPosteriorEstimates(SEXP inParams);
        
        //Destructor
        ~spatialSEIRModel();
    private:
        int ncalls;
        dataModel* dataModelInstance;
        exposureModel* exposureModelInstance;
        reinfectionModel* reinfectionModelInstance;
        distanceModel* distanceModelInstance;
        transitionPriors* transitionPriorsInstance;
        initialValueContainer* initialValueContainerInstance;
        samplingControl* samplingControlInstance;
        scoped_actor* self;
};

#endif
