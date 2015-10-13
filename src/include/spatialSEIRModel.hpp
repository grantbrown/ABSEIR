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

struct samplingResultSet
{
    Rcpp::NumericVector result;
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
    Eigen::MatrixXd rEA;
    //Eigen::MatrixXd r0t;
    //Eigen::MatrixXd effR0;
    Eigen::MatrixXd p_se;
    Eigen::MatrixXd p_ei;
    Eigen::MatrixXd p_ir;
    Eigen::MatrixXd rho;
    Eigen::MatrixXd beta;
    double result;
};

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
        /** Constructor */
        spatialSEIRModel(dataModel& dataModel_,
                         exposureModel& exposureModel_,
                         reinfectionModel& reinfectionModel_,
                         distanceModel& distanceModel_,
                         transitionPriors& transitionPriors_,
                         initialValueContainer& initialValueContainer_,
                         samplingControl& samplingControl_);
        /** The main ABC function - draw nSample samples from the approximated
         * posterior, and optionally set verbose to a nonzero integer for additional
         * output.*/ 
        Rcpp::List sample(SEXP nSample, SEXP verbose);
        /** Accept a matrix of parameters, simulate epidemics, and return the 
         * euclidean distance from the outcome. */
        Rcpp::List evaluate(SEXP inParams);
        /** Simulate epdiemics based on the matrix of parameters inParams,
         * return the simulated epidemic quantities.  
         */
        Rcpp::List simulate_given(SEXP inParams);
        /** Continue sampling from the approximated posterior after an inital 
         * call to sample. 
         */
        Rcpp::List update(SEXP nSample, SEXP inParams, SEXP inEps, 
                SEXP inWeights, SEXP inCurEps, SEXP verbose);
        double evalPrior(Rcpp::NumericVector param_values);
        /** Destructor */
        ~spatialSEIRModel();

    private:
        /** Simulated epidemics are processed in batches for parallelizeability. 
         * This applies both to the BasicABC and Beaumont2009 algorithms.*/
        int batchNum;
        /** An indicator as to whether or not reweighting is required 
         * (Beaumont2009)*/
        int reweight;
        /** A counter of how many times this object has been asked to simulate
         * data. This is not currently used much, as we destroy objects after
         * use and recreate them as needed.*/
        int ncalls;
        /** A flag to indicate that the Beaumont2009 sampler should make 
         * perterbation proposals from a multivariate normal distribution
         * rather than independant normals.
         */
        bool multivariatePerturbation;

        /** A persistant reference to the number of requested samples*/
        int numSamples;
        /** The fraction of the current estimates which was replaced by new 
         * draws in the latest batch (BasicABC)*/
        double updateFraction;
        /** The minimum distance estimate currently accepted*/
        double minEps;
        /** The maximum distance estimate currently accepted*/
        double maxEps;
        /** The currently enforced uppder bound on distance*/
        double currentEps;
        /** A vector of current parameter means*/
        Eigen::RowVectorXd parameterMeans;
        /** Centered values of current parameters */
        Eigen::MatrixXd parameterCentered;
        /** Current parameter covariance matrix*/
        Eigen::MatrixXd parameterCov;
        /** Current parameter inverse covariance matrix*/
        Eigen::MatrixXd parameterICov; 
        /** Current parameter inverse covariance matrix determinant*/
        double parameterICovDet;
        /** internal implementation of the main sampling function*/
        Rcpp::List sample_internal(int nSample, bool verbose, bool init);
        /** Generic function to propose new parameters. */
        void updateParams();
        /** Function to propose new parameters from the prior distribution*/
        void updateParams_prior();
        /** Function to propose new parameters accordin to the weights from
         * Beaumont 2009*/
        void updateParams_SMC();
        /** Function to compute new weights from Beaumont 2009*/
        void updateWeights();
        /** Storage vector for accepted parameters during repeated sampling.*/
        std::vector<Eigen::VectorXd> currentAccepted;
        /** Storage vector for accepted parameter distances 
         * during repeated sampling.*/
        std::vector<double> currentAcceptedResult;
        /** Matrix from which parameters are sent to the workers. */
        Eigen::MatrixXd param_matrix;
        /** Vector of Gaussian standard deviations used in forward SMC kernel*/
        Eigen::VectorXd tau;
        /** The current resampling weights (Beaumont 2009) */
        Eigen::VectorXd weights;
        /** A data structure containing current epoch parameter values and 
         * distances.*/
        samplingResultSet currentSamples;
        /** A data structure containing previous epoch parameter values and 
         * distances, used in weight computation.*/
        samplingResultSet previousSamples;

        /** Main simulation function. */
        Rcpp::List simulate(Eigen::MatrixXd params, caf::atom_value sim_type);

        /** General function to take new samples and use them to update the 
         * currently accepted ones. */
        samplingResultSet combineResults(Rcpp::NumericVector currentResults, 
                                  Rcpp::NumericMatrix currentParams,
                                  Rcpp::NumericVector newResults,
                                  Eigen::MatrixXd newParams);

        /** Parameter update function for BasicABC rejection algorithm */
        samplingResultSet combineResults_basic(Rcpp::NumericVector currentResults, 
                                  Rcpp::NumericMatrix currentParams,
                                  Rcpp::NumericVector newResults,
                                  Eigen::MatrixXd newParams);
        /** Parameter update function for Beaumont2009 algorithm */
        samplingResultSet combineResults_SMC(
                                  Rcpp::NumericVector newResults,
                                  Eigen::MatrixXd newParams);

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
        /** Pointer to a scoped_actor object - can't we delete this?*/
        scoped_actor* self;
        /** A persistant pointer to a properly initialized random 
         * number generator.*/
        std::mt19937* generator;
};

#endif
