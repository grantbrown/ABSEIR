#include <Rcpp.h>
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
#include "caf/all.hpp"

using namespace Rcpp;
using namespace caf;


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

    // Set up random number provider 
    std::minstd_rand0 lc_generator(samplingControlInstance -> random_seed + 1);
    std::uint_least32_t seed_data[std::mt19937::state_size];
    std::generate_n(seed_data, std::mt19937::state_size, std::ref(lc_generator));
    std::seed_seq q(std::begin(seed_data), std::end(seed_data));
    generator = new std::mt19937{q};   
}

samplingResultSet spatialSEIRModel::combineResults(Rcpp::NumericVector currentResults, 
                                            Rcpp::NumericMatrix currentParams,
                                            Rcpp::NumericVector newResults,
                                            Eigen::MatrixXd newParams)
{
    Rcpp::NumericVector outResults = Rcpp::clone(currentResults);
    Rcpp::NumericMatrix outParams = Rcpp::clone(currentParams);
    std::vector<size_t> currentIndex = sort_indexes(currentResults);
    std::vector<size_t> newIndex = sort_indexes(newResults);
    size_t idx1 = 0;
    size_t idx2 = 0;
    size_t i, j;

    // Zipper merge
    for (i = 0; i < currentIndex.size(); i++)
    {
        if (currentResults(currentIndex[idx1]) > newResults(newIndex[idx2]))
        {
            outResults(i) = newResults(newIndex[idx2]);
            for (j = 0; j < (size_t) currentParams.ncol(); j++)
            {
                outParams(i,j) = newParams(newIndex[idx2], j);
            } 
            idx2++;
        }
        else if (currentResults(currentIndex[idx1]) >= newResults(newIndex[idx2]))
        {
            outResults(i) = currentResults(currentIndex[idx1]);
            for (j = 0; j < (size_t) currentParams.ncol(); j++)
            {
                outParams(i,j) = currentResults(currentIndex[idx1], j);
            }
            idx1++;
        }
    }
    samplingResultSet output;
    output.result = outResults;
    output.params = outParams;
    return(output);
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
    // Back up current results for later weight calculation.
    previousSamples.params = Rcpp::clone(currentSamples.params);
    previousSamples.result = Rcpp::clone(currentSamples.result);

    Eigen::VectorXd newValues(currentSamples.params.nrow());
    std::uniform_real_distribution<double> runif(0.0,1.0);
    std::vector<std::normal_distribution<double>> K;

    for (i = 0; i < currentSamples.params.ncol(); i++)
    {
        tau(i) = 2.0*Rcpp::sd(currentSamples.params( _, i));
        K.push_back(std::normal_distribution<double>(0.0, tau(i)));
    }

    // Calculate cumulative weights for resampling. 
    std::vector<double> cumulativeWeights(currentSamples.params.nrow());
    cumulativeWeights[0] = weights(0);
    for (i = 1; i < currentSamples.params.rows(); i++)
    {
        cumulativeWeights[i] = cumulativeWeights[i-1] + weights(i);
    }

    int up;
    for (i = 0; i < currentSamples.params.nrow(); i++)
    {
        up = std::upper_bound(cumulativeWeights.begin(), 
                              cumulativeWeights.end(), runif(*generator)) 
            - cumulativeWeights.begin();

        newValues(i) = currentSamples.result(up);
        for (j = 0; j < param_matrix.cols(); j ++)
        {
            param_matrix(i, j) = currentSamples.params(up, j) + 
                K[j](*generator);
        }
    }
}
    
void spatialSEIRModel::updateParams_prior()
{
    const bool hasReinfection = (reinfectionModelInstance -> 
            betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;

    int nBeta = (exposureModelInstance -> X).cols();
    int nBetaRS = (reinfectionModelInstance -> X_rs).cols()*hasReinfection;
    int nRho = (distanceModelInstance -> dm_list).size()*hasSpatial;
    int bs = samplingControlInstance -> batch_size;
    int i, j;

    // Set up random samplers 
    // beta, beta_RS
    std::normal_distribution<> standardNormal(0,1); 
    // rho  
    std::gamma_distribution<> rhoDist(
            (distanceModelInstance -> spatial_prior)(0),
        1.0/(distanceModelInstance -> spatial_prior)(1));
    // gamma_ei 
    std::gamma_distribution<> gammaEIDist(
            (transitionPriorsInstance -> gamma_ei_params)(0),
        1.0/(transitionPriorsInstance -> gamma_ei_params)(1)); 
    // gamma_ir
    std::gamma_distribution<> gammaIRDist(
            (transitionPriorsInstance -> gamma_ir_params)(0),
        1.0/(transitionPriorsInstance -> gamma_ir_params)(1)); 

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
        // Draw gamma_ei
        param_matrix(i, nBeta + nBetaRS + nRho) = gammaEIDist(*generator);
        // Draw gamma_ir
        param_matrix(i, nBeta + nBetaRS + nRho + 1) = gammaIRDist(*generator);
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

void spatialSEIRModel::updateWeights()
{
    if ((samplingControlInstance -> algorithm) 
            == ALG_ModifiedBeaumont2009 && batchNum > 0)
    {
        int i,j,k;
        int N = currentSamples.params.nrow();
        int nParams = currentSamples.params.ncol();
        double tmpWeight;
        double totalWeight = 0.0;
        Eigen::VectorXd newWeights(N);
        for (i = 0; i < N; i++)
        {
            newWeights(i) = 0.0;
            for (j = 0; j < N; j++) 
            {
                tmpWeight = 0.0;
                for (k = 0; k < nParams; k++)
                {
                    tmpWeight += pow(currentSamples.params(i, k) 
                            - currentSamples.params(j, k), 2)/(2.0*tau[k]);
                }
                newWeights(i) += weights(j)*std::exp(-1.0*tmpWeight);
            }
            totalWeight += newWeights(i);
        }

        for (i = 0; i < N; i++)
        {
            newWeights(i) /= totalWeight;
            weights(i) = newWeights(i);
        }
    }

}

Rcpp::List spatialSEIRModel::sample(SEXP nSamples)
{
    Rcpp::Rcout << "Sampling\n";
    Rcpp::IntegerVector nSamp(nSamples);

    bool hasReinfection = (reinfectionModelInstance -> betaPriorPrecision)(0) > 0;
    bool hasSpatial = (dataModelInstance -> Y).cols() > 1;

    int N = nSamp(0);
    double r = samplingControlInstance -> accept_fraction;
    int bs = samplingControlInstance -> batch_size;
    int i;
    int nBeta = (exposureModelInstance -> X).cols();
    int nBetaRS = (reinfectionModelInstance -> X_rs).cols()*hasReinfection;
    int nRho = (distanceModelInstance -> dm_list).size()*hasSpatial;
    int nParams = nBeta + nBetaRS + nRho + 2;

    int nBatches = (samplingControlInstance -> algorithm == ALG_BasicABC ? 
                        std::ceil(((1.0*N)/r)/bs) :  
                        (samplingControlInstance -> epochs));

    tau = Eigen::VectorXd(nParams);

    batchNum = 0;

    if (samplingControlInstance -> algorithm == ALG_ModifiedBeaumont2009)
    {
        weights = Eigen::VectorXd(N);
        for (i = 0; i < N; i++)
        {
            weights(i) = 1.0/N;
        }
    }


    Rcpp::List tmpList, outList;

    Rcpp::NumericVector outputValues(N);
    for (i = 0; i < N; i++)
    {
        outputValues(i) = std::numeric_limits<double>::infinity(); 
    }

    Rcpp::NumericMatrix outputParams(nSamp(0), nParams);

    currentSamples.result = outputValues;
    currentSamples.params = outputParams;

    param_matrix = Eigen::MatrixXd(bs, nParams);
    for (batchNum = 0; batchNum < nBatches; batchNum ++)
    {
        updateParams();
        tmpList = this -> simulate(param_matrix, sample_atom::value);
        currentSamples = combineResults(currentSamples.result, 
                                        currentSamples.params,
                                        as<NumericVector>(tmpList["result"]),
                                        param_matrix);
        updateWeights();
    }
    outList["result"] = currentSamples.result;
    outList["params"] = currentSamples.params;
    return(outList);
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
    return(this -> simulate(param_matrix, sim_atom::value));
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

    auto worker_pool = actor_pool::make(actor_pool::round_robin());
    unsigned int ncore = (unsigned int) samplingControlInstance -> CPU_cores;
    unsigned int nrow =  (unsigned int) param_matrix.rows(); 
    unsigned int i;

    for (i = 0; i < ncore; i++)
    {
        workers.push_back((*self) -> spawn<SEIR_sim_node, monitored>(samplingControlInstance->simulation_width,
                                                                      samplingControlInstance->random_seed + 1000*(i + 1) + ncalls,
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
        outList["result"] = outResults;
    }
    delete self;
    return(outList);
}

spatialSEIRModel::~spatialSEIRModel()
{   
    shutdown();
    delete generator;
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

