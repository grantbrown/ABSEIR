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
#include <SEIRSimNodes.hpp>

using namespace Rcpp;

std::vector<size_t> sort_indexes(Rcpp::NumericVector inVec)
{
    vector<size_t> idx(inVec.size());
    for (size_t i = 0; i < idx.size(); i++)
    {
        idx[i] = i;
    }
    std::sort(idx.begin(), idx.end(),
         [&inVec](size_t i1, size_t i2){return(inVec(i1) < inVec(i2));});
    return(idx);    
}

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
        Rcpp::stop("Error: model components were not provided in the correct order. \n");
    }

    ncalls = 0;
    updateFraction = 0;
    minEps = 0.0;
    maxEps = 0.0;
    numSamples = 0;

    dataModelInstance = &dataModel_;
    exposureModelInstance = &exposureModel_;
    reinfectionModelInstance = &reinfectionModel_;
    distanceModelInstance = &distanceModel_;
    transitionPriorsInstance = &transitionPriors_;
    initialValueContainerInstance = &initialValueContainer_;
    samplingControlInstance = &samplingControl_;

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
    currentEps = std::numeric_limits<double>::infinity();

    // Set up random number provider 
    std::minstd_rand0 lc_generator(samplingControlInstance -> random_seed + 1);
    std::uint_least32_t seed_data[std::mt19937::state_size];
    std::generate_n(seed_data, std::mt19937::state_size, std::ref(lc_generator));
    std::seed_seq q(std::begin(seed_data), std::end(seed_data));
    generator = new std::mt19937{q};   


    // TODO: there should be some better way to synchronize this 
    // so that we don't need to do all sorts of locking-pushing-sorting
    result_idx = std::vector<int>();
    // Having two results containers is kinda ugly, better solution?
    results_complete = std::vector<simulationResultSet>();
    results_double = std::vector<double>();
    // Prefill results_double
    for (int idx = 0; idx < samplingControlInstance -> batch_size; idx++)
    {
        results_double.push_back(0.0);
    }

    std::vector<int>* index_pointer = &result_idx;

    unsigned int ncore = (unsigned int) samplingControlInstance -> CPU_cores;

    worker_pool = std::unique_ptr<NodePool>(
                new NodePool(&results_double,
                     &results_complete,
                     index_pointer,
                     ncore,
                     samplingControlInstance->random_seed + ncalls,
                     initialValueContainerInstance -> S0,
                     initialValueContainerInstance -> E0,
                     initialValueContainerInstance -> I0,
                     initialValueContainerInstance -> R0,
                     exposureModelInstance -> offset,
                     dataModelInstance -> Y,
                     dataModelInstance -> na_mask,
                     distanceModelInstance -> dm_list,
                     distanceModelInstance -> tdm_list,
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
                     dataModelInstance -> cumulative
                ));
}
    

samplingResultSet spatialSEIRModel::combineResults(
                                            Rcpp::NumericVector currentResults, 
                                            Rcpp::NumericMatrix currentParams,
                                            Rcpp::NumericVector newResults,
                                            Eigen::MatrixXd newParams)
{
    if ((batchNum == 0) || samplingControlInstance -> algorithm == ALG_BasicABC)
    {
        return(combineResults_basic(currentResults, currentParams, newResults, 
                                    newParams));
    }
    return(combineResults_SMC(newResults, newParams));
}

samplingResultSet spatialSEIRModel::combineResults_basic(
                                            Rcpp::NumericVector currentResults, 
                                            Rcpp::NumericMatrix currentParams,
                                            Rcpp::NumericVector newResults,
                                            Eigen::MatrixXd newParams)
{
    // If we need to expand result sets, this is where it happens. 
    Rcpp::NumericVector outResults = Rcpp::NumericVector(numSamples);
    Rcpp::NumericMatrix outParams = Rcpp::NumericMatrix(numSamples, 
                                        currentParams.ncol());

    std::vector<size_t> currentIndex = sort_indexes(currentResults);
    std::vector<size_t> newIndex = sort_indexes(newResults);
    size_t idx1 = 0;
    size_t idx2 = 0;
    size_t err = 0;
    size_t i, j;


    // Zipper merge
    for (i = 0; i < (size_t) numSamples; i++)
    {
        // Make sure we don't run out of samples
        if (idx1 >= (size_t) currentResults.size() || 
                currentResults(currentIndex[idx1]) > newResults(newIndex[idx2]))
        {
            outResults(i) = newResults(newIndex[idx2]);
            for (j = 0; j < (size_t) currentParams.ncol(); j++)
            {
                outParams(i,j) = newParams(newIndex[idx2], j);
            } 
            idx2++;
        }
        else if (currentResults(currentIndex[idx1]) <= newResults(newIndex[idx2]))
        {
            outResults(i) = currentResults(currentIndex[idx1]);
            for (j = 0; j < (size_t) currentParams.ncol(); j++)
            {
                outParams(i,j) = currentParams(currentIndex[idx1], j);
            }
            idx1++;
        }
        else
        {   
            err += 1;
        }
    }
    minEps = outResults(0);
    maxEps = outResults(outResults.size()-1);
    currentEps = maxEps;
    updateFraction = (idx2*1.0)/(idx1+idx2+err);
    samplingResultSet output;
    output.result = outResults;
    output.params = outParams;
    return(output);
}

samplingResultSet spatialSEIRModel::combineResults_SMC(
                                            Rcpp::NumericVector newResults,
                                            Eigen::MatrixXd newParams)
{
    int idx = 0;
    int i,j;
    int nrow = newResults.size();
    int N = numSamples;
    while (idx < nrow && currentAccepted.size() < (size_t) N) 
    {
        if (newResults(idx) < currentEps)
        {
            currentAccepted.push_back(newParams.row(idx));
            currentAcceptedResult.push_back(newResults(idx));
        }
        idx++;
    }

    double minRslt = std::numeric_limits<double>::infinity();
    double maxRslt = 0;
    if (currentAccepted.size() == (size_t) N)
    {
        reweight = 1;
        samplingResultSet output;
        output.result = Rcpp::NumericVector(N);
        output.params = Rcpp::NumericMatrix(N, newParams.cols());
        for (i = 0; i < output.params.nrow(); i++)
        {
            for (j = 0; j < output.params.ncol(); j++)
            {
                output.params(i,j) = currentAccepted[i](j);
            }
            output.result(i) = currentAcceptedResult[i]; 
            minRslt = std::min(output.result(i), minRslt);
            maxRslt = std::max(output.result(i), maxRslt);
        }
        minEps = minRslt;
        maxEps = maxRslt;
        currentEps = (samplingControlInstance -> shrinkage)*currentEps;
        currentAccepted.clear();
        currentAcceptedResult.clear();
        return(output);
    }
    reweight = 0;
    return(currentSamples);
}


void spatialSEIRModel::updateParams()
{
    if (batchNum == 0 || ((samplingControlInstance -> algorithm) == ALG_BasicABC))
    {
        updateParams_prior();
    }
    else if ((samplingControlInstance -> algorithm) == ALG_ModifiedBeaumont2009)
    {
        updateParams_SMC();
    }
    else
    {
        Rcpp::stop("Unknown algorithm number");
    }
}

void spatialSEIRModel::updateParams_SMC()
{
    int i,j;
    int bs = param_matrix.rows();
    int csSize = currentSamples.params.nrow();
    int nParams = currentSamples.params.ncol();
    // Back up current results for later weight calculation.
    previousSamples.params = Rcpp::clone(currentSamples.params);
    previousSamples.result = Rcpp::clone(currentSamples.result);

    Eigen::MatrixXd L;
    Eigen::VectorXd Z;
    std::uniform_real_distribution<double> runif(0.0,1.0);
    std::normal_distribution<double> rnorm(0.0, 1.0);
    std::vector<std::normal_distribution<double> > K;



    if (samplingControlInstance -> multivariatePerturbation)
    {
        // TODO: Integrate with/replace above
        Eigen::Map<Eigen::MatrixXd> curSamp(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(currentSamples.params));
        Eigen::RowVectorXd parameterMeans = curSamp.colwise().mean();
        parameterCentered = (curSamp.rowwise() - parameterMeans);
        Eigen::MatrixXd parameterCov = (parameterCentered).transpose() 
                * (parameterCentered)/(curSamp.rows()-1);
        Eigen::LLT<Eigen::MatrixXd> paramLLT = parameterCov.llt();
        L = Eigen::MatrixXd(paramLLT.matrixL());
        Z = Eigen::VectorXd(parameterCov.cols());
        parameterICov = paramLLT.solve(Eigen::MatrixXd::Identity(L.rows(), L.cols()));
        parameterICovDet = 1.0/std::pow(L.diagonal().prod(), 2.0); 
    }
    else
    {
        for (i = 0; i < nParams; i++)
        {
            tau(i) = Rcpp::sd(currentSamples.params( _, i));
            K.push_back(std::normal_distribution<double>(0.0, tau(i)));
        }
    }

    // Calculate cumulative weights for resampling. 
    std::vector<double> cumulativeWeights(csSize);
    cumulativeWeights[0] = weights(0);
    for (i = 1; i < csSize; i++)
    {
        cumulativeWeights[i] = cumulativeWeights[i-1] + weights(i);
    }
    if (std::abs(cumulativeWeights[csSize-1]-1) > 1e-6)
    {
        Rcpp::Rcout << "Warning: weights don't add up to 1: "
            << cumulativeWeights[csSize-1] << "\n";
    }
    cumulativeWeights[csSize-1] = 1.0;

    int up, itrs;
    bool hasValidProposal;
    Rcpp::NumericVector proposal(nParams);
    for (i = 0; i < bs; i++)
    {
        itrs = 0;
        up = std::upper_bound(cumulativeWeights.begin(), 
                              cumulativeWeights.end(), runif(*generator)) 
            - cumulativeWeights.begin();
        hasValidProposal = false;
        if (samplingControlInstance -> multivariatePerturbation)
        {
            while (!hasValidProposal && itrs < 10000)
            {
                // Fill up standard normal vector
                for (j = 0; j < param_matrix.cols(); j++)
                {
                    Z(j) = rnorm(*generator);
                }
                // Generate multivariate normal
                proposal = Rcpp::wrap(L * Z);
                // Add back appropriate mean
                proposal += currentSamples.params(up, _);
                hasValidProposal = (evalPrior(proposal) != 0);  
                itrs++;
            }
        }
        else
        {
            while (!hasValidProposal && itrs < 10000)
            {
                for (j = 0; j < param_matrix.cols(); j ++)
                {
                    proposal(j) = currentSamples.params(up, j) + K[j](*generator);
                }
                // proposals with impossible values receive 0 weight anyway.
                hasValidProposal = (evalPrior(proposal) != 0);
                itrs++;
            }
        }
        // Check that we got a valid proposal, otherwise panic
        if (!hasValidProposal)
        {
            Rcpp::Rcout << "Proposal Iterations: " << itrs << "\n";
            for (j = 0; j < proposal.size(); j++)
            {
                Rcpp::Rcout << proposal(j) << ", ";
            }
            Rcpp::Rcout << "\n";
            Rcpp::stop("Unable to generate a valid proposal.\n");
        }
        // Assign accepted proposal to param_matrix
        for (j = 0; j < nParams; j++)
        {
            param_matrix(i, j) = proposal(j);
        }
    }
}
    
void spatialSEIRModel::updateParams_prior()
{
    const bool hasReinfection = (reinfectionModelInstance -> 
            betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;
    std::string transitionMode = transitionPriorsInstance -> mode;

    const int nBeta = (exposureModelInstance -> X).cols();
    const int nBetaRS = (reinfectionModelInstance -> X_rs).cols()*hasReinfection;
    const int nRho = ((distanceModelInstance -> dm_list).size() + 
                      (distanceModelInstance -> tdm_list).size() > 0 ? 
                      (distanceModelInstance -> tdm_list)[0].size() : 0)*hasSpatial;

    const int bs = samplingControlInstance -> batch_size;
    int i, j;

    // Set up random samplers 
    // beta, beta_RS
    std::normal_distribution<> standardNormal(0,1); 
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
    for (i = 0; i < bs; i++)
    {
        // Draw beta
        for (j = 0; j < nBeta; j++)
        {
            param_matrix(i, j) = (exposureModelInstance -> betaPriorMean(j)) + 
                                 standardNormal(*generator) /
                                 (exposureModelInstance -> betaPriorPrecision(j));
        }
    }
    // draw gammaEI, gammaIR
    if (transitionMode == "exponential")
    {
        for (i = 0; i < bs; i++)
        {
            // Draw gamma_ei
            param_matrix(i, nBeta + nBetaRS + nRho) = 
                gammaEIDist[0](*generator);
            // Draw gamma_ir
            param_matrix(i, nBeta + nBetaRS + nRho + 1) = 
                gammaIRDist[0](*generator);
        }
    }
    else if (transitionMode == "weibull")
    {
        for (i = 0; i < bs; i++)
        {
            param_matrix(i, nBeta + nBetaRS + nRho) = gammaEIDist[0](*generator);
            param_matrix(i, nBeta + nBetaRS + nRho + 1) = gammaEIDist[1](*generator);
            param_matrix(i, nBeta + nBetaRS + nRho + 2) = gammaIRDist[0](*generator);
            param_matrix(i, nBeta + nBetaRS + nRho + 3) = gammaIRDist[1](*generator);
        }
    }
    // Draw reinfection parameters
    if (hasReinfection)
    {
        for (i = 0; i < bs; i++)
        {
            for (j = nBeta; j < nBeta + nBetaRS; j++)
            {
                param_matrix(i, j) = 
                    (reinfectionModelInstance -> betaPriorMean(j)) + 
                     standardNormal(*generator) /
                    (reinfectionModelInstance -> betaPriorPrecision(j));           
            }
        }
    }
    // Draw rho
    if (hasSpatial)
    {
        for (i = 0; i < bs; i++)
        {
            rhoTot = 2.0; 
            rhoItrs = 0;
            while (rhoTot > 1.0 && rhoItrs < 100)
            {
                rhoTot = 0.0;
                for (j = nBeta + nBetaRS; j < nBeta + nBetaRS + nRho; j++)
                {
                   param_matrix(i, j) = rhoDist(*generator); 
                   rhoTot += param_matrix(i,j);
                }
                rhoItrs++;
            }
            if (rhoTot > 1.0)
            {
                Rcpp::Rcout << "Error, valid rho value not obtained\n";
            }
        }
    }

}

double spatialSEIRModel::evalPrior(Rcpp::NumericVector param_vector)
{
    double outPrior = 1.0;
    double constr = 0.0;
    const bool hasReinfection = (reinfectionModelInstance -> betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;
    std::string transitionMode = transitionPriorsInstance -> mode;
    const int nBeta = (exposureModelInstance -> X).cols();
    const int nBetaRS = (reinfectionModelInstance -> X_rs).cols()*hasReinfection;
    const int nRho = ((distanceModelInstance -> dm_list).size() + 
                      (distanceModelInstance -> tdm_list).size() > 0 ? 
                      (distanceModelInstance -> tdm_list)[0].size() : 0)*hasSpatial;

    int i;

    int paramIdx = 0;
    for (i = 0; i < nBeta; i++)
    {
        outPrior *= R::dnorm(param_vector(paramIdx), 
                (exposureModelInstance -> betaPriorMean)(i), 
                1.0/((exposureModelInstance -> betaPriorPrecision)(i)), 0);
        paramIdx++;
    }

    if (nBetaRS > 0)
    {
        for (i = 0; i < nBetaRS; i++)
        {
            outPrior *= R::dnorm(param_vector(paramIdx), 
                    (reinfectionModelInstance -> betaPriorMean)(i), 
                    1.0/((reinfectionModelInstance -> betaPriorPrecision)(i)), 0);
            paramIdx++;
        }
    }

    if (nRho > 0)
    {
        for (i = 0; i < nRho; i++)
        {
             constr += param_vector(paramIdx);
             outPrior *= R::dbeta(param_vector(paramIdx), 
                         (distanceModelInstance -> spatial_prior)(0),
                         (distanceModelInstance -> spatial_prior)(1), 0);
             paramIdx++;
        }
        outPrior *= (constr <= 1);
    }

    if (transitionMode == "exponential")
    {
        outPrior *= R::dgamma(param_vector(paramIdx), 
                (transitionPriorsInstance -> E_to_I_params)(0,0),
                1.0/(transitionPriorsInstance -> E_to_I_params)(1,0), 0);
        paramIdx++;

        outPrior *= R::dgamma(param_vector(paramIdx), 
                (transitionPriorsInstance -> I_to_R_params)(0,0),
                1.0/(transitionPriorsInstance -> I_to_R_params)(1,0), 0);
    }
    else if (transitionMode == "weibull")
    {
        // The following may be slow. Consider 
        // 1. override of evalParamPrior for Weibull case
        // 2. rewrite of evalParamPrior to accept Rcpp object
        Eigen::Map<Eigen::VectorXd> eigen_param_vec( 
            Rcpp::as<Eigen::Map<Eigen::VectorXd> >(param_vector));
        outPrior *= EI_transition_dist -> evalParamPrior(
                eigen_param_vec.segment(paramIdx, 2)); 
        paramIdx += 2;
        outPrior *= IR_transition_dist -> evalParamPrior(
                eigen_param_vec.segment(paramIdx, 2)); 
        paramIdx += 2;
    }
    return(outPrior);
}

void spatialSEIRModel::updateWeights()
{
    if ((samplingControlInstance -> algorithm) 
            == ALG_ModifiedBeaumont2009 
            && reweight == 1)
    {
        reweight = 0;
        int i,j,k;
        int N = currentSamples.params.nrow();
        // If we need to expand the weight vector, this is where it happens
        int oldN = weights.size(); 
        int nParams = currentSamples.params.ncol();
        double tmpWeightComp;
        double totalWeight;
        double tmpWeight;
        Eigen::VectorXd newWeights(N);
        Eigen::VectorXd tmpParams;
        if (samplingControlInstance -> multivariatePerturbation)
        {
            tmpWeightComp = (-nParams/2.0)*std::log(2.0*3.14159) 
                - 0.5*std::log((parameterICovDet));
            totalWeight = 0.0;
            for (i = 0; i < N; i++)
            {
                newWeights(i) = 0.0;
                for (j = 0; j < oldN; j++) 
                {
                    tmpWeight = tmpWeightComp;
                    tmpWeight -= (0.5*
                            (parameterCentered.row(i)
                            *parameterICov*parameterCentered.transpose()
                                .col(i))(0));
                    newWeights(i) += weights(j)*std::exp(tmpWeight);
                }
                newWeights(i) = evalPrior(currentSamples.params(i, _))
                    /newWeights(i);
                totalWeight += newWeights(i);
            }
        }
        else
        {
            tmpWeightComp = 0.0;
            for (k = 0; k < nParams; k++)
            {
                tmpWeightComp -= (0.5*std::log(2*3.14159) + std::log(tau[k]));  
            }
            totalWeight = 0.0;
            for (i = 0; i < N; i++)
            {
                newWeights(i) = 0.0;
                for (j = 0; j < oldN; j++) 
                {
                    tmpWeight = tmpWeightComp;
                    for (k = 0; k < nParams; k++)
                    {
                        tmpWeight -= 1.0*std::pow(currentSamples.params(i,k)
                                                - previousSamples.params(j,k), 
                                                2)/(2.0*std::pow(tau[k], 2));
                    }

                    newWeights(i) += weights(j)*std::exp(tmpWeight);
                }
                newWeights(i) = evalPrior(currentSamples.params(i, _))/newWeights(i);
                if (!std::isfinite(newWeights(i))){
                    Rcpp::stop("newWeights(i) isn't finite 2");
                }

                totalWeight += newWeights(i);
            }
        }
        if (!std::isfinite(totalWeight)){
            Rcpp::stop("Non-finite weights were generated. ");
        }

        weights = newWeights;
        for (i = 0; i < N; i++)
        {
            weights(i) = weights(i)/totalWeight;
        }
    }

}

Rcpp::List spatialSEIRModel::sample_internal(int N, bool verbose, bool init)
{
    numSamples = N;
    const bool hasReinfection = (reinfectionModelInstance -> betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;
    std::string transitionMode = transitionPriorsInstance -> mode;

    const double target_eps = samplingControlInstance -> target_eps;
    const double r = samplingControlInstance -> accept_fraction;
    const int bs = samplingControlInstance -> batch_size;
    int i;
    const int nBeta = (exposureModelInstance -> X).cols();
    const int nBetaRS = (reinfectionModelInstance -> X_rs).cols()*hasReinfection;
    const int nRho = ((distanceModelInstance -> dm_list).size() + 
                      (distanceModelInstance -> tdm_list).size() > 0 ? 
                      (distanceModelInstance -> tdm_list)[0].size() : 0)*hasSpatial;
    const int nTrans = (transitionMode == "exponential" ? 2 :
                       (transitionMode == "weibull" ? 4 : 0));
    const int nParams = nBeta + nBetaRS + nRho + nTrans;


    const int nBatches = (samplingControlInstance -> algorithm == ALG_BasicABC ? 
                        std::ceil(((1.0*N)/r)/bs) :  
                        (samplingControlInstance -> epochs));

    if (bs < N)
    {
        Rcpp::stop("Simulation batch size must be at least as large as final sample size\n");
    }

    tau = Eigen::VectorXd(nParams);

    batchNum = (init ? 1 : 0); // If we're bringing in existing samples/weights, 
                               // then don't do the usual batch 0 stuff. 
    reweight = batchNum;

    // Initialize weights if we don't already have them. 
    // Only used for Beaumont 2009 algorithm, but set and returned for 
    // consistency. 
    if (!init)
    {
        weights = Eigen::VectorXd(N);
        for (i = 0; i < N; i++)
        {
            weights(i) = 1.0/N;
        }
    }

    // Declare output containers
    Rcpp::List tmpList, outList;

    Rcpp::NumericVector outputValues(N);
    for (i = 0; i < N; i++)
    {
        outputValues(i) = std::numeric_limits<double>::infinity(); 
    }

    Rcpp::NumericMatrix outputParams(N, nParams);

    if (!init)
    {
        currentSamples.result = outputValues;
        currentSamples.params = outputParams;
    }

    int incompleteBatches = 0;
    param_matrix = Eigen::MatrixXd(bs, nParams);
    while (batchNum < nBatches && incompleteBatches < 
            (samplingControlInstance -> max_batches) && 
            currentEps/(samplingControlInstance -> shrinkage) > target_eps)
    {
        updateParams();
        tmpList = this -> simulate(param_matrix, sim_atom);
        currentSamples = combineResults(currentSamples.result, 
                                        currentSamples.params,
                                        as<NumericVector>(tmpList["result"]),
                                        param_matrix);

        if (((samplingControlInstance -> algorithm) 
                == ALG_ModifiedBeaumont2009) && reweight == 0
                && batchNum != 0)
        {
            incompleteBatches++;
        }
        else 
        {
            incompleteBatches = 0;
            batchNum++;
        }

        if (verbose && (samplingControlInstance -> algorithm) 
                == ALG_ModifiedBeaumont2009)
        {
            if (reweight != 0 || (!init && batchNum == 1 
                        && incompleteBatches == 0))
            {
                Rcpp::Rcout << "Completed batch " << batchNum << " of " << 
                    nBatches << ". Eps: [" << 
                    minEps << ", " << maxEps << "] < " << 
                    currentEps/(samplingControlInstance -> shrinkage) << "\n";
            }
            else
            {
                Rcpp::Rcout << "Incomplete batch number: " 
                    << incompleteBatches
                    << " of max " << (samplingControlInstance -> max_batches)
                    << ". Current size: " 
                    << currentAccepted.size() << " of " << N << ".\n"; 
            }
        }
        else if (verbose)
        {
            Rcpp::Rcout << "Completed batch " << batchNum << " of " << 
                nBatches << ". Eps: [" << 
                minEps << ", " << maxEps << "]\n";
        }
        Rcpp::checkUserInterrupt();
        updateWeights();
    }

    Rcpp::NumericVector outputWeights(N);
    for (i = 0; i < N; i++)
    {
        outputWeights(i) = weights(i); 
    }

    outList["result"] = currentSamples.result;
    outList["params"] = currentSamples.params;
    outList["weights"] = outputWeights;
    outList["currentEps"] = (incompleteBatches == 
            (samplingControlInstance -> max_batches) 
            ? currentEps/(samplingControlInstance 
                                -> shrinkage) 
            : currentEps); // If we didn't get a complete batch in, prev eps
    outList["completedEpochs"] = batchNum;
    return(outList);
}

Rcpp::List spatialSEIRModel::sample(SEXP nSamples, SEXP vb)
{
    Rcpp::IntegerVector nSamp(nSamples);
    Rcpp::IntegerVector vbs(vb);

    bool verbose = (vbs(0) != 0);
    int N = nSamp(0);
    return(sample_internal(N, verbose, false));
}

Rcpp::List spatialSEIRModel::evaluate(SEXP inParams)
{
    Rcpp::NumericMatrix params(inParams);
    // Copy to Eigen matrix
    int i, j;
    Eigen::MatrixXd param_matrix(params.nrow(), params.ncol());
    for (i = 0; i < params.nrow(); i++)
    {
        for (j = 0; j < params.ncol(); j++)
        {
            param_matrix(i,j) = params(i,j);
        }
    }
    return(this -> simulate(param_matrix, sim_atom));
}

Rcpp::List spatialSEIRModel::simulate_given(SEXP inParams)
{
    Rcpp::NumericMatrix params(inParams);
    // Copy to Eigen matrix
    int i, j;
    Eigen::MatrixXd param_matrix(params.nrow(), params.ncol());
    for (i = 0; i < params.nrow(); i++)
    {
        for (j = 0; j < params.ncol(); j++)
        {
            param_matrix(i,j) = params(i,j);
        }
    }
    return(this -> simulate(param_matrix, sim_result_atom));
}

Rcpp::List spatialSEIRModel::update(SEXP nSample, SEXP inParams,
        SEXP inEps, SEXP inWeights, SEXP inCurEps, SEXP inVerbose)
{
    Rcpp::IntegerVector nSamp(nSample);
    Rcpp::IntegerVector vbs(inVerbose);

    bool verbose = (vbs(0) != 0);
    int N = nSamp(0);

    Rcpp::NumericMatrix params(inParams);
    Rcpp::NumericVector eps(inEps);
    Rcpp::NumericVector wts(inWeights);
    Rcpp::NumericVector curEps(inCurEps);

    currentSamples.result = eps;
    currentSamples.params = params;

    int i;
    // Need for expansion.
    // Need for initialization 
    minEps = 0.0;
    maxEps = 0.0;
    currentEps = curEps(0);
    
    weights = Eigen::VectorXd(wts.size());
    for (i = 0; i < wts.size(); i++)
    {
        if (!std::isfinite(wts(i)))
        {
            Rcpp::stop("Non-finite weights encountered.\n");
        }
        weights(i) = wts(i);
    }
    return(sample_internal(N, verbose, true));
}

Rcpp::List spatialSEIRModel::simulate(Eigen::MatrixXd param_matrix, 
                                      std::string sim_type_atom)
{
    const bool hasReinfection = (reinfectionModelInstance -> 
            betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;
    std::string transitionMode = transitionPriorsInstance -> mode;

    if (!(sim_type_atom == sim_atom || 
         sim_type_atom == sim_result_atom))
    {
        Rcpp::Rcout << "Invalid simulation type requested\n";
        Rcpp::List outList;
        outList["error"] = true;
        return(outList);
    }

    ncalls += 1;    

    unsigned int i;
    unsigned int nrow = (unsigned int) param_matrix.rows(); 

    // TODO: there should be some better way to synchronize this 
    // so that we don't need to do all sorts of locking-pushing-sorting
    result_idx.clear();
    results_complete.clear();
    // Don't clear results_double

    Eigen::VectorXd outRow;
    // Send simulation orders
    for (i = 0; i < nrow; i++)
    {
        outRow = param_matrix.row(i);
        worker_pool -> enqueue(sim_type_atom, i, outRow);
    }

    worker_pool -> awaitFinished();
   
    // Todo: keep an eye on this object handling. It may have unreasonable
    // overhead, and is kind of complex.  
    Rcpp::List outList;

    if (sim_type_atom == sim_result_atom)
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
            // We p_ei and p_ir not generally defined in non-exponential case.  
            if (transitionMode == "exponential")
            {
                subList["p_ei"] = createRcppNumericFromEigen(results_complete[i].p_ei);
                subList["p_ir"] = createRcppNumericFromEigen(results_complete[i].p_ir);
            }
            subList["R_EA"] = createRcppNumericFromEigen(results_complete[i].rEA);
            subList["R0t"] = createRcppNumericFromEigen(results_complete[i].r0t);
            subList["effR0t"] = createRcppNumericFromEigen(results_complete[i].effR0);
            if (hasSpatial)
            {
                subList["rho"] = createRcppNumericFromEigen(results_complete[i].rho);
            }
            subList["beta"] = createRcppNumericFromEigen(results_complete[i].beta);
            subList["X"] = createRcppNumericFromEigen(results_complete[i].X);
            if (hasReinfection)
            {
                // TODO: output reinfection info
            }
            subList["result"] = results_complete[i].result;
            outList[std::to_string(i)] = subList;
        }
    }
    else if (sim_type_atom == sim_atom)
    {
        Rcpp::NumericVector outResults(nrow);
        for (i = 0; i < nrow; i++)
        {
            // No resorting
            outResults(i) = results_double[i];
        }
        outList["result"] = outResults;
    }
    return(outList);
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
    .method("evaluate", &spatialSEIRModel::evaluate)
    .method("sample", &spatialSEIRModel::sample)
    .method("update", &spatialSEIRModel::update)
    .method("simulate", &spatialSEIRModel::simulate_given);
}

