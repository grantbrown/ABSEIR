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

using namespace Rcpp;

double rbeta(double a, double b, std::mt19937* generator){
    double x = std::gamma_distribution<double>(a,1)(*generator);
    double y = std::gamma_distribution<double>(b,1)(*generator);
    return(x/(x+y));
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
        Rcpp::stop("Error: model components were not provided in the correct order. \n");
    }

    // Store model components
    dataModelInstance = &dataModel_;
    exposureModelInstance = &exposureModel_;
    reinfectionModelInstance = &reinfectionModel_;
    distanceModelInstance = &distanceModel_;
    transitionPriorsInstance = &transitionPriors_;
    initialValueContainerInstance = &initialValueContainer_;
    samplingControlInstance = &samplingControl_;

    // Check for model component compatibility
    if ((dataModelInstance -> nLoc) != (exposureModelInstance -> nLoc))
    { 
        Rcpp::stop(("Exposure model and data model imply different number of locations: " 
                + std::to_string(dataModelInstance -> nLoc) + ", " 
                + std::to_string(exposureModelInstance -> nLoc) + ".\n").c_str());
    }
    if ((dataModelInstance -> nTpt) != (exposureModelInstance -> nTpt))
    { 
        Rcpp::stop(("Exposure model and data model imply different number of time points:"
                    + std::to_string(dataModelInstance -> nTpt) + ", "
                    + std::to_string(exposureModelInstance -> nTpt) + ".\n").c_str());  
    }
    if ((dataModelInstance -> nLoc) != (distanceModelInstance -> numLocations))
    {       
        Rcpp::stop(("Data model and distance model imply different number of locations:"
                    + std::to_string(dataModelInstance -> nLoc) + ", "
                    + std::to_string(distanceModelInstance -> numLocations) + ".\n").c_str()
                );
    }
    if ((int) (distanceModelInstance -> tdm_list).size() != (dataModelInstance -> nTpt))
    {
        Rcpp::stop("TDistance model and data model imply a different number of time points.\n");
    }
    int sz1 = (distanceModelInstance -> tdm_list)[0].size();
    for (int i = 0; i < (int) (distanceModelInstance -> tdm_list).size(); i++)
    {
        if ((int) distanceModelInstance -> tdm_list[i].size() != sz1)
        {
            Rcpp::stop("Differing number of lagged contact matrices across time points.\n");
        }
    }
    if ((dataModelInstance -> nLoc) != (initialValueContainerInstance -> S0.size())) 
    { 
        Rcpp::stop("Data model and initial value container have different dimensions\n");
    }
    if ((reinfectionModelInstance -> reinfectionMode) == 3)
    {
        // No reinfection
    }
    else
    {
        if (((reinfectionModelInstance -> X_rs).rows()) != (dataModelInstance -> nTpt))
        { 
            Rcpp::stop("Reinfection and data mode time points differ.\n");
        }
    }
    if ((reinfectionModelInstance -> reinfectionMode) > 2)
    {
        // pass
    }

    if (transitionPriorsInstance -> mode != "exponential" && 
            transitionPriorsInstance -> mode != "path_specific" &&
            transitionPriorsInstance -> mode != "weibull")
    {
        Rcpp::stop("Invalid transition mode: " + 
                (transitionPriorsInstance -> mode));
    }

    // Optionally, set up transition distribution
    if (transitionPriorsInstance -> mode == "weibull")
    {
        EI_transition_dist = std::unique_ptr<weibullTransitionDistribution>(
                new weibullTransitionDistribution(
                (transitionPriorsInstance -> E_to_I_params).col(0))); 
        IR_transition_dist = std::unique_ptr<weibullTransitionDistribution>(
                new weibullTransitionDistribution(
                (transitionPriorsInstance -> I_to_R_params).col(0))); 
    }
    else
    {
        Eigen::VectorXd DummyParams(4);
        DummyParams(0) = 1.0;
        DummyParams(1) = 1.0;
        DummyParams(2) = 1.0;
        DummyParams(3) = 1.0;
        EI_transition_dist = std::unique_ptr<weibullTransitionDistribution>(new 
            weibullTransitionDistribution(DummyParams));
        IR_transition_dist = std::unique_ptr<weibullTransitionDistribution>(new 
            weibullTransitionDistribution(DummyParams));
    }
    // Set up param matrix
    const bool hasReinfection = (reinfectionModelInstance -> 
            betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;
    std::string transitionMode = transitionPriorsInstance -> mode;

    const int nBeta = (exposureModelInstance -> X).cols();
    const int nBetaRS = (reinfectionModelInstance -> X_rs).cols()*hasReinfection;
    const int nRho = ((distanceModelInstance -> dm_list).size() + 
                      (distanceModelInstance -> tdm_list)[0].size())*hasSpatial;
    const int nTrans = (transitionMode == "exponential" ? 2 :
                       (transitionMode == "weibull" ? 4 : 0));
    const int nReport = (dataModelInstance -> dataModelType == 2 ? 1 : 0);

    const int nParams = nBeta + nBetaRS + nRho + nTrans + nReport;

    // Set up random number provider 
    std::minstd_rand0 lc_generator(samplingControlInstance -> random_seed + 1);
    std::uint_least32_t seed_data[std::mt19937::state_size];
    std::generate_n(seed_data, std::mt19937::state_size, std::ref(lc_generator));
    std::seed_seq q(std::begin(seed_data), std::end(seed_data));
    generator = new std::mt19937{q};   

    // Parameters are not initialized
    is_initialized = false;

    results_complete = std::vector<simulationResultSet>();
    results_double = Eigen::MatrixXd::Zero(samplingControlInstance -> init_batch_size, 
                                               samplingControlInstance -> m); 
    param_matrix = Eigen::MatrixXd::Zero(samplingControlInstance -> init_batch_size, 
                                            nParams);

    // Parameter covariance matrix and related items are not initialized

    parameterCov = Eigen::MatrixXd::Zero(nParams, nParams);
    parameterICov = Eigen::MatrixXd::Zero(nParams, nParams);
    parameterL = Eigen::MatrixXd::Zero(nParams, nParams);
    parameterICovDet = 0.0;

    result_idx = std::vector<int>();

    // Create the worker pool
    worker_pool = std::unique_ptr<NodePool>(
                new NodePool(&results_double,
                     &results_complete,
                     &result_idx,
                     (unsigned int) samplingControlInstance -> CPU_cores,
                     samplingControlInstance->random_seed,
                     initialValueContainerInstance -> S0,
                     initialValueContainerInstance -> E0,
                     initialValueContainerInstance -> I0,
                     initialValueContainerInstance -> R0,
                     exposureModelInstance -> offset,
                     dataModelInstance -> Y,
                     dataModelInstance -> na_mask,
                     dataModelInstance -> dataModelType,
                     distanceModelInstance -> dm_list,
                     distanceModelInstance -> tdm_list,
                     distanceModelInstance -> tdm_empty,
                     exposureModelInstance -> X,
                     reinfectionModelInstance -> X_rs,
                     transitionPriorsInstance -> mode,
                     transitionPriorsInstance -> E_to_I_params,
                     transitionPriorsInstance -> I_to_R_params,
                     transitionPriorsInstance -> inf_mean,
                     distanceModelInstance -> spatial_prior,
                     exposureModelInstance -> betaPriorPrecision,
                     reinfectionModelInstance -> betaPriorPrecision, 
                     exposureModelInstance -> betaPriorMean,
                     reinfectionModelInstance -> betaPriorMean,
                     dataModelInstance -> phi,
                     dataModelInstance -> dataModelCompartment,
                     dataModelInstance -> cumulative,
                     samplingControlInstance -> m
                ));
}

Eigen::MatrixXd spatialSEIRModel::generateParamsPrior(int nParticles)
{
    const bool hasReinfection = (reinfectionModelInstance -> 
            betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;
    std::string transitionMode = transitionPriorsInstance -> mode;

    const int nBeta = (exposureModelInstance -> X).cols();
    const int nBetaRS = (reinfectionModelInstance -> X_rs).cols()*hasReinfection;
    const int nRho = ((distanceModelInstance -> dm_list).size() + 
                      (distanceModelInstance -> tdm_list)[0].size())*hasSpatial;
    const int nTrans = (transitionMode == "exponential" ? 2 :
                       (transitionMode == "weibull" ? 4 : 0));
    const int nReport = (dataModelInstance -> dataModelType == 2 ? 1 : 0);
    const int nParams = nBeta + nBetaRS + nRho + nTrans + nReport;


    int N = nParticles;

    int i, j;

    double rf_alpha = (dataModelInstance -> dataModelType == 2  ?
                       (dataModelInstance -> report_fraction)*(dataModelInstance -> report_fraction_ess) :
                       -1.0);
    double rf_beta = (dataModelInstance -> dataModelType == 2 ?
                       (1.0 - dataModelInstance -> report_fraction)*(dataModelInstance -> report_fraction_ess) :
                       -1.0);



    Eigen::MatrixXd outParams = Eigen::MatrixXd::Zero(N, nParams);

    // Set up random samplers 
    // beta, beta_RS
    std::normal_distribution<double> standardNormal(0,1); 
    // rho  
    std::gamma_distribution<> rhoDist(
            (distanceModelInstance -> spatial_prior)(0),
        1.0/(distanceModelInstance -> spatial_prior)(1));
    // Hyperprior distributions for E to I and I to R transitions
    std::vector<std::gamma_distribution<> > gammaEIDist;
    std::vector<std::gamma_distribution<> > gammaIRDist;

    if (transitionMode == "exponential")
    {
        gammaEIDist.push_back(std::gamma_distribution<>(
                (transitionPriorsInstance -> E_to_I_params)(0,0),
            1.0/(transitionPriorsInstance -> E_to_I_params)(1,0)));
        gammaIRDist.push_back(std::gamma_distribution<>(
                (transitionPriorsInstance -> I_to_R_params)(0,0),
            1.0/(transitionPriorsInstance -> I_to_R_params)(1,0)));
    }   
    else if (transitionMode == "weibull")
    {
        gammaEIDist.push_back(std::gamma_distribution<>(
                (transitionPriorsInstance -> E_to_I_params)(0,0),
            1.0/(transitionPriorsInstance -> E_to_I_params)(1,0)));
        gammaEIDist.push_back(std::gamma_distribution<>(
                (transitionPriorsInstance -> E_to_I_params)(2,0),
            1.0/(transitionPriorsInstance -> E_to_I_params)(3,0)));

        gammaIRDist.push_back(std::gamma_distribution<>(
                (transitionPriorsInstance -> I_to_R_params)(0,0),
            1.0/(transitionPriorsInstance -> I_to_R_params)(1,0)));
        gammaIRDist.push_back(std::gamma_distribution<>(
                (transitionPriorsInstance -> I_to_R_params)(2,0),
            1.0/(transitionPriorsInstance -> I_to_R_params)(3,0)));
    }
    // If this is too slow, consider column-wise operations
    double rhoTot = 0.0;
    int rhoItrs = 0;
    for (i = 0; i < nParticles; i++)
    {
        // Draw beta
        for (j = 0; j < nBeta; j++)
        {
            outParams(i, j) = (exposureModelInstance -> betaPriorMean(j)) + 
                                 standardNormal(*generator) /
                                 (exposureModelInstance -> betaPriorPrecision(j));
        }
    }
    // draw gammaEI, gammaIR
    if (transitionMode == "exponential")
    {
        for (i = 0; i < nParticles; i++)
        {
            // Draw gamma_ei
            outParams(i, nBeta + nBetaRS + nRho) = 
                gammaEIDist[0](*generator);
            // Draw gamma_ir
            outParams(i, nBeta + nBetaRS + nRho + 1) = 
                gammaIRDist[0](*generator);
        }
    }
    else if (transitionMode == "weibull")
    {
        for (i = 0; i < nParticles; i++)
        {
            outParams(i, nBeta + nBetaRS + nRho) = gammaEIDist[0](*generator);
            outParams(i, nBeta + nBetaRS + nRho + 1) = gammaEIDist[1](*generator);
            outParams(i, nBeta + nBetaRS + nRho + 2) = gammaIRDist[0](*generator);
            outParams(i, nBeta + nBetaRS + nRho + 3) = gammaIRDist[1](*generator);
        }
    }
    // Draw reinfection parameters
    if (hasReinfection)
    {
        for (i = 0; i < nParticles; i++)
        {
            for (j = nBeta; j < nBeta + nBetaRS; j++)
            {
                outParams(i, j) = 
                    (reinfectionModelInstance -> betaPriorMean(j-nBeta)) + 
                     standardNormal(*generator) /
                    (reinfectionModelInstance -> betaPriorPrecision(j-nBeta));           
            }
        }
    }
    // Draw rho
    if (hasSpatial)
    {
        for (i = 0; i < nParticles; i++)
        {
            rhoTot = 2.0; 
            rhoItrs = 0;
            while (rhoTot > 1.0 && rhoItrs < 100)
            {
                rhoTot = 0.0;
                for (j = nBeta + nBetaRS; j < nBeta + nBetaRS + nRho; j++)
                {
                   outParams(i, j) = rhoDist(*generator); 
                   rhoTot += outParams(i,j);
                }
                rhoItrs++;
            }
            if (rhoTot > 1.0)
            {
                Rcpp::Rcout << "Error, valid rho value not obtained\n";
            }
        }
    }
    // Draw report fraction
    if (dataModelInstance -> dataModelType == 2)
    {
        int lastcol = outParams.cols() -1;
        for (i = 0; i < nParticles; i++)
        {
            outParams(i,lastcol) = rbeta(rf_alpha, rf_beta, generator);
        }
    }
    return(outParams);
}

Rcpp::List spatialSEIRModel::sample(SEXP nSample, SEXP returnComps, SEXP verbose)
{
    Rcpp::IntegerVector n(nSample);
    Rcpp::IntegerVector r(returnComps);
    Rcpp::IntegerVector v(verbose);

    int N = n(0);
    bool R = r(0) > 0;
    int V = v(0);

    if (R && V)
    {
        Rcpp::Rcout << "return compartments: true\n";
    }

    std::string sim_type_atom = (R ? sim_result_atom : sim_atom);
    
    if (samplingControlInstance -> algorithm == ALG_BasicABC)
    {
        return(sample_basic(N, V, sim_type_atom));
    }
    else if (samplingControlInstance -> algorithm == ALG_ModifiedBeaumont2009)
    {
        return(sample_Beaumont2009(N, V, sim_type_atom));
    }
    else if (samplingControlInstance -> algorithm == ALG_DelMoral2012)
    {
        return(sample_DelMoral2012(N, V, sim_type_atom));
    }
    else 
    {
	if (samplingControlInstance -> algorithm != ALG_Simulate)
	{
            Rcpp::stop("Unknown algorithm.");
	}
        if (!is_initialized)
        {
            Rcpp::stop("Model must be initialized before simulating.");
        }
        return(sample_Simulate(N, 0, V));
    }
}

bool spatialSEIRModel::setParameters(Eigen::MatrixXd params, 
        Eigen::VectorXd weights, Eigen::MatrixXd results, double eps)
{
    const bool hasReinfection = (reinfectionModelInstance ->                    
                        betaPriorPrecision)(0) > 0; 
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;                
    std::string transitionMode = transitionPriorsInstance -> mode;
    const int nBeta = (exposureModelInstance -> X).cols();
    const int nBetaRS = (reinfectionModelInstance -> X_rs).cols()*hasReinfection;
    const int nRho = ((distanceModelInstance -> dm_list).size() + 
                      (distanceModelInstance -> tdm_list)[0].size())*hasSpatial;

    const int nTrans = (transitionMode == "exponential" ? 2 :
                       (transitionMode == "weibull" ? 4 : 0));
    const int nReport = (dataModelInstance -> dataModelType == 2 ? 1 : 0);

    const int nParams = nBeta + nBetaRS + nRho + nTrans + nReport;
    
    if (params.cols() != nParams)
    {
        Rcpp::stop("Number of supplied parameters does not match model specification.\n");
    }
    if (params.rows() != weights.size())
    {
        Rcpp::stop("Number of weights not equal to number of particles.\n");
    }

    init_eps = eps;

    init_param_matrix = params; 
    param_matrix = params;

    init_results_double = results;
    results_double = results;

    init_weights = weights;
    is_initialized = true;
    return(true);
}

double spatialSEIRModel::evalPrior(Eigen::VectorXd param_vector)
{
    double outPrior = 0.0;
    double constr = 0.0;
    const bool hasReinfection = (reinfectionModelInstance -> betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;
    std::string transitionMode = transitionPriorsInstance -> mode;
    const int nBeta = (exposureModelInstance -> X).cols();
    const int nBetaRS = (reinfectionModelInstance -> X_rs).cols()*hasReinfection;
    const int nRho = ((distanceModelInstance -> dm_list).size() + 
                      (distanceModelInstance -> tdm_list)[0].size())*hasSpatial;
    //const int nReport = (dataModelInstance -> dataModelType == 2 ? 1 : 0);
    double rf_alpha = (dataModelInstance -> dataModelType == 2 ?
                       (dataModelInstance -> report_fraction)*(dataModelInstance -> report_fraction_ess) :
                       -1.0);
    double rf_beta = (dataModelInstance -> dataModelType == 2 ?
                       (1.0 - dataModelInstance -> report_fraction)*(dataModelInstance -> report_fraction_ess) :
                       -1.0);



    int i;
    int paramIdx = 0;
    for (i = 0; i < nBeta; i++)
    {
        outPrior += R::dnorm(param_vector(paramIdx), 
                (exposureModelInstance -> betaPriorMean)(i), 
                1.0/((exposureModelInstance -> betaPriorPrecision)(i)), 1);
        paramIdx++;
    }

    if (nBetaRS > 0)
    {
        for (i = 0; i < nBetaRS; i++)
        {
            outPrior += R::dnorm(param_vector(paramIdx), 
                    (reinfectionModelInstance -> betaPriorMean)(i), 
                    1.0/((reinfectionModelInstance -> betaPriorPrecision)(i)), 1);
            paramIdx++;
        }
    }

    if (nRho > 0)
    {
        for (i = 0; i < nRho; i++)
        {
             constr += param_vector(paramIdx);
             outPrior += R::dbeta(param_vector(paramIdx), 
                         (distanceModelInstance -> spatial_prior)(0),
                         (distanceModelInstance -> spatial_prior)(1), 1);
             paramIdx++;
        }
        if (constr > 1){
            outPrior = -std::numeric_limits<double>::infinity();
        }
    }

    if (transitionMode == "exponential")
    {
        outPrior += R::dgamma(param_vector(paramIdx), 
                (transitionPriorsInstance -> E_to_I_params)(0,0),
                1.0/(transitionPriorsInstance -> E_to_I_params)(1,0), 1);
        paramIdx++;

        outPrior += R::dgamma(param_vector(paramIdx), 
                (transitionPriorsInstance -> I_to_R_params)(0,0),
                1.0/(transitionPriorsInstance -> I_to_R_params)(1,0), 1);
    }
    else if (transitionMode == "weibull")
    {
        outPrior += EI_transition_dist -> evalParamPrior(
                param_vector.segment(paramIdx, 2)); 
        paramIdx += 2;
        outPrior += IR_transition_dist -> evalParamPrior(
                param_vector.segment(paramIdx, 2)); 
        paramIdx += 2;
    }
    if (dataModelInstance -> dataModelType == 2)
    {
        outPrior += R::dbeta(param_vector(param_vector.size()-1), 
                            rf_alpha,
                            rf_beta,
                            1);
    }
    return(std::exp(outPrior));
}

void spatialSEIRModel::run_simulations(Eigen::MatrixXd params, 
                                       std::string sim_type_atom,
                                       Eigen::MatrixXd* results_dest,
                                       std::vector<simulationResultSet>* results_c_dest)
{
    int i;
    worker_pool -> setResultsDest(results_dest, 
                                  results_c_dest);
    for (i = 0; i < params.rows(); i++)
    {
        worker_pool -> enqueue(sim_type_atom, i, params.row(i));
    }
    worker_pool -> awaitFinished();
}

spatialSEIRModel::~spatialSEIRModel()
{   
    delete generator;
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
    .method("sample", &spatialSEIRModel::sample)
    .method("setParameters", &spatialSEIRModel::setParameters);
}

