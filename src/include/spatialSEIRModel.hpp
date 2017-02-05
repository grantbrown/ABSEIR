#ifndef SPATIALSEIR_MODEL_INC
#define SPATIALSEIR_MODEL_INC

#include <Rcpp.h>
#include <memory>
#include "./dataModel.hpp"
#include "./distanceModel.hpp"
#include "./exposureModel.hpp"
#include "./initialValueContainer.hpp"
#include "./reinfectionModel.hpp"
#include "./samplingControl.hpp"
#include "./SEIRSimNodes.hpp"
#include "./transitionPriors.hpp"
#include "./transitionDistribution.hpp"

struct samplingResultSet
{
    Rcpp::NumericMatrix result;
    Rcpp::NumericMatrix params;
};

struct simulationResultSet
{
    Eigen::MatrixXi S;
    Eigen::MatrixXi E;
    Eigen::MatrixXi I;
    Eigen::MatrixXi R;
    Eigen::MatrixXi S_star;
    Eigen::MatrixXi E_star;
    Eigen::MatrixXi I_star;
    Eigen::MatrixXi R_star;
    Eigen::MatrixXd X;
    Eigen::MatrixXd p_se;
    Eigen::MatrixXd p_ei;
    Eigen::MatrixXd p_ir;
    Eigen::MatrixXd rho;
    Eigen::MatrixXd beta;
    Eigen::MatrixXd result; 
};

class dataModel;
class exposureModel;
class distanceModel;
class initialValueContainer;
class reinfectionModel;
class samplingControl;
class transitionPriors;
class NodePool;

using namespace Rcpp;

class spatialSEIRModel
{
    public: 
        /** Constructor */
        spatialSEIRModel(dataModel& dataModel_,
                         exposureModel& exposureModel_,
                         reinfectionModel& reinfectionModel_,
                         distanceModel& distanceModel_,
                         transitionPriors& transitionPriors_,
                         initialValueContainer& initialValueContainer_,
                         samplingControl& samplingControl_);
        /** The main ABC function - draw nSample samples from the approximated
         * posterior, and optionally set verbose to 0, 1, or 2 for 
         * different levels of output.*/ 
        Rcpp::List sample(SEXP nSample, SEXP returnComps, SEXP verbose);
        /** Evaluate the prior distribution of a particular set of parameters*/
        double evalPrior(Eigen::VectorXd param_values);
        /** Assign the parameter values manually */
        bool setParameters(Eigen::MatrixXd param_values, 
                           Eigen::VectorXd weights,
                           Eigen::MatrixXd results,
                           double eps);

        /** Public for caching purposes. Consider factoring accessors back 
         * into class methods*/
        /* Parameter sample covariance matrix*/
        Eigen::MatrixXd parameterCov;

        /* Parameter inverse sample covariance matrix*/
        Eigen::MatrixXd parameterICov;

        /* Parameter sample covariance matrix*/
        double parameterICovDet;

        /* Parameter sample covariance matrix cholesky factor*/
        Eigen::MatrixXd parameterL;

        /** Destructor */
        ~spatialSEIRModel();

    private:
        /** Set parameters from prior distribution*/
        Eigen::MatrixXd generateParamsPrior(int N);

        /** Simulate epidemics based on parameters*/
        void run_simulations(Eigen::MatrixXd params, 
                             std::string sim_type_atom,
                             Eigen::MatrixXd* result_recip,
                             std::vector<simulationResultSet>* result_c_recip);

        /** Run simulation using basic ABC algorithm */
        Rcpp::List sample_basic(int nSample, int verbose, 
                                std::string sim_type_atom);

        /** Run simulation using Beaumont 2009 algorithm */
        Rcpp::List sample_Beaumont2009(int nSample, int verbose, 
                                std::string sim_type_atom);

        /** Run simulation using Del Moral 2012 algorithm */
        Rcpp::List sample_DelMoral2012(int nSample, int verbose, 
                                std::string sim_type_atom);

        /** Use current parameters to simulate epidemics*/
        Rcpp::List sample_Simulate(int nSample, int enforceEps, int verbose);

        /** Flag for whether params have been initialized*/
        bool is_initialized;

        /** If simulation is re-started, need initial epsilon stored*/
        double init_eps;

        /** If simulation is re-started, need initial weights stored*/
        Eigen::VectorXd init_weights;

        /** If simulation is re-started, need initial params stored */
        Eigen::MatrixXd init_param_matrix;

        /** If simulation is re-started, need initial results stored */
        Eigen::MatrixXd init_results_double;

        /** General E to I transition Distribution*/
        std::unique_ptr<transitionDistribution> EI_transition_dist;

        /** General I to R transition Distribution*/
        std::unique_ptr<transitionDistribution> IR_transition_dist;

        /** Matrix of parameters */
        Eigen::MatrixXd proposed_param_matrix;

        /** Matrix of parameters */
        Eigen::MatrixXd param_matrix;

        /** Matrix of parameters */
        Eigen::MatrixXd prev_param_matrix;  

        /** Result index vector */
        std::vector<int> result_idx;

        /** Results vector*/
        Eigen::MatrixXd results_double;

        /** Results vector*/
        Eigen::MatrixXd prev_results_double;

        /** Results vector*/
        Eigen::MatrixXd proposed_results_double;

        /** particles - cache*/
        Eigen::MatrixXd proposal_cache;

        /** particles - cache*/
        Eigen::MatrixXd preproposal_params;

        /** particles - cache*/
        Eigen::MatrixXd preproposal_results;

        /** Complete results vector */
        std::vector<simulationResultSet> results_complete;

        /** Complete results vector */
        std::vector<simulationResultSet> proposed_results_complete;

        /** Pointer to a dataModel object*/
        dataModel* dataModelInstance;

        /** Pointer to an exposureModel object*/
        exposureModel* exposureModelInstance;

        /** Pointer to a reinfectionModel object*/
        reinfectionModel* reinfectionModelInstance;

        /** Pointer to a distanceModel object*/
        distanceModel* distanceModelInstance;

        /** Pointer to a transitionPriors object*/
        transitionPriors* transitionPriorsInstance;

        /** Pointer to initialValueContainer object*/
        initialValueContainer* initialValueContainerInstance;

        /** Pointer to a samplingControl object.*/
        samplingControl* samplingControlInstance;

        /** Thread pool */
        std::unique_ptr<NodePool> worker_pool; 

        /** A persistant pointer to a properly initialized random 
         * number generator.*/
        std::mt19937* generator;
};

#endif
