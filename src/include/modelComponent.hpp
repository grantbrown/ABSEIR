#ifndef SPATIALSEIR_MODEL_COMPONENT
#define SPATIALSEIR_MODEL_COMPONENT

#define LSS_INVALID_CONTAINER 0
#define LSS_DATA_MODEL_TYPE 1
#define LSS_EXPOSURE_MODEL_TYPE 2
#define LSS_REINFECTION_MODEL_TYPE 3
#define LSS_SAMPLING_CONTROL_MODEL_TYPE 4
#define LSS_DISTANCE_MODEL_TYPE 5
#define LSS_TRANSITION_MODEL_TYPE 6
#define LSS_INIT_CONTAINER_TYPE 7


#include <Rcpp.h>
#include <Eigen/Core>


using namespace Rcpp;

class modelComponent
{
    public:
        virtual ~modelComponent(){};
        virtual int getModelComponentType() = 0;
        virtual void summary() = 0;
        void protect(){prot=1;}
        void unprotect(){prot=0;}
        int prot=0;

};


#endif
