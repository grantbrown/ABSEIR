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


                                                                                
std::vector<size_t> sort_indexes_eigen(Eigen::MatrixXd inMat)                     
{                                                                               
        vector<size_t> idx(inMat.rows());                                           
        for (size_t i = 0; i < idx.size(); i++)                                     
        {                                                                           
            idx[i] = i;                                                             
        }                                                                           
        std::sort(idx.begin(), idx.end(),                                           
                         [&inMat](size_t i1, size_t i2){return(inMat(i1,0) < inMat(i2,0));});       
            return(idx);                                                                
}                                                                               
   

void proposeParams_beaumont_multivariate(Eigen::MatrixXd* params,
                                         Eigen::VectorXd* tau,
                                         std::mt19937* generator,
                                         spatialSEIRModel* model)
{
    int i, j;
    int N = params -> rows();
    Eigen::MatrixXd curSamp = *params;
    Eigen::RowVectorXd parameterMeans = curSamp.colwise().mean();
    auto parameterCentered = (curSamp.rowwise() - parameterMeans);
    Eigen::MatrixXd parameterCov = (parameterCentered).transpose() 
            * (parameterCentered)/(N-1);
    // Add a SD to overdisperse particles:
    for (i = 0; i < parameterCov.cols(); i++) 
    {
        parameterCov(i,i) += 0.5*(*tau)(i);
    }
    Eigen::LLT<Eigen::MatrixXd> paramLLT = parameterCov.llt();
    auto L = Eigen::MatrixXd(paramLLT.matrixL());
    auto parameterICov = paramLLT.solve(Eigen::MatrixXd::Identity(L.rows(), L.cols()));
    auto parameterICovDet = 1.0/std::pow(L.diagonal().prod(), 2.0); 
    auto Z = Eigen::VectorXd(parameterCov.cols());
    //auto parameterICov = paramLLT.solve(Eigen::MatrixXd::Identity(L.rows(), L.cols()));
    //auto parameterICovDet = 1.0/std::pow(L.diagonal().prod(), 2.0); 
    Eigen::VectorXd proposal = (params -> row(0))*0.0; 

    bool hasValid;
    int itrs = 0;
    Eigen::VectorXd param_cache = params -> row(0);
    for (i = 0; i < N; i++)
    {
        itrs = 0;
        hasValid = false;
        param_cache = params -> row(i);

        while (!hasValid && itrs < 1000)
        {
            params -> row(i) = param_cache;
            // Fill up standard normal vector
            for (j = 0; j < params -> cols(); j++)
            {
                Z(j) = std::normal_distribution<double>(0.0,1.0)(*generator);
            }
            // Generate multivariate normal
            proposal = L * Z;
            // Add back appropriate mean
            params -> row(i) += proposal;
            hasValid = (model -> evalPrior(params -> row(i)) > 0);  
            itrs++;
        }
        if (!hasValid)
        {
            Rcpp::Rcout << "Unable to generate parameters with nonzero probability.\n";
            Rcpp::Rcout << "  Param " << i << " of " << params -> rows() << "\n"; 
            Rcpp::Rcout << "  Pror prob: " << (model -> evalPrior(params -> row(i)));
            Rcpp::Rcout << "  Param: \n" << params -> row(i) << "\n";
            Rcpp::stop("No parameters");
        }
    }
    // Store important components for later reuse (optimization)
    model -> parameterCov = parameterCov;
    model -> parameterICov = parameterICov;
    model -> parameterICovDet = parameterICovDet;
    model -> parameterL = L;
}


void proposeParams_beaumont(Eigen::MatrixXd* params,
                            Eigen::VectorXd* tau,
                            std::mt19937* generator,
                            spatialSEIRModel* model)
{
    int i,j;
    // Propose new parameters
    int p = params -> cols();
    int N = params -> rows();
    bool hasValid;
    int itrs = 0;
    Eigen::VectorXd param_cache = params -> row(0);
    for (i = 0; i < N; i++)
    {
        itrs = 0;
        hasValid = false;
        param_cache = params -> row(i);
        for (itrs = 0; itrs < 1000 &&  !hasValid; itrs++)
        { 
            params -> row(i) = param_cache;
            for (j = 0; j < p; j++)
            {
                auto propDist = std::normal_distribution<double>(0.0, (*tau)(j));
                (*params)(i,j) += propDist(*generator);
            }
            hasValid = (model -> evalPrior(params -> row(i))) > 0;
        }
        if (!hasValid)
        {
            Rcpp::Rcout << "Unable to generate parameters with nonzero probability.\n";
            Rcpp::Rcout << "  Param " << i << " of " << params -> rows() << "\n"; 
            Rcpp::Rcout << "  Pror prob: " << (model -> evalPrior(params -> row(i)));
            Rcpp::Rcout << "  Param: \n" << params -> row(i) << "\n";
            Rcpp::stop("No parameters");
        }
    }
}

Rcpp::List spatialSEIRModel::sample_Beaumont2009(int nSample, int vb, 
                                                 std::string sim_type_atom)
{
    // This is set by constructor
    const int nParams = param_matrix.cols();

    // Accepted params/results are nSample size
    results_complete = std::vector<simulationResultSet>();
    results_double = Eigen::MatrixXd::Zero(nSample, 
                                           samplingControlInstance -> m); 
    param_matrix = Eigen::MatrixXd::Zero(nSample, 
                                            nParams);
    prev_param_matrix = param_matrix;
    // Proposals matrices are batch size
    proposed_results_double = Eigen::MatrixXd::Zero(nSample,
                                               samplingControlInstance -> m); 
    proposed_param_matrix = Eigen::MatrixXd::Zero(nSample,
                                                  nParams);
    proposal_cache = Eigen::MatrixXd::Zero(samplingControlInstance -> batch_size, 
                                           nParams);
    preproposal_params = Eigen::MatrixXd::Zero(samplingControlInstance -> batch_size, 
                                               nParams);
    preproposal_results = Eigen::MatrixXd::Zero(samplingControlInstance -> batch_size, 
                                                samplingControlInstance -> m); 
 
    bool terminate = false;
    const int verbose = vb;
    if (verbose > 1)
    {
        Rcpp::Rcout << "Starting sampler\n";
    }
    const int num_iterations = samplingControlInstance -> epochs;
    const int Nsim = samplingControlInstance -> batch_size;
    const int Npart = nSample;

    const int maxBatches= samplingControlInstance -> max_batches;
    const bool hasReinfection = (reinfectionModelInstance -> 
               betaPriorPrecision)(0) > 0;
    const bool hasSpatial = (dataModelInstance -> Y).cols() > 1;

    std::string transitionMode = transitionPriorsInstance -> mode;   

    double e0 = std::numeric_limits<double>::infinity();
    double e1 = std::numeric_limits<double>::infinity();

    int i,j,k;
    int iteration;

    auto U = std::uniform_real_distribution<double>(0,1);
    double drw;

    if (verbose > 1)
    {
        Rcpp::Rcout << "Number of iterations requested: " 
                    << num_iterations << "\n";

        dataModelInstance -> summary();
        exposureModelInstance -> summary();
        reinfectionModelInstance -> summary();
        distanceModelInstance -> summary();
        initialValueContainerInstance -> summary();
        samplingControlInstance -> summary();

    }
    // Step 0b: set weights to 1/N
    Eigen::VectorXd w0 = Eigen::VectorXd::Zero(Npart).array() + 1.0/((double) Npart);
    Eigen::VectorXd w1 = Eigen::VectorXd::Zero(Npart).array() + 1.0/((double) Npart);
    Eigen::VectorXd cum_weights = Eigen::VectorXd::Zero(Npart).array() + 1.0/((double) Npart);

    
    if (!is_initialized)
    {
        if (verbose > 1){Rcpp::Rcout << "Generating starting parameters from prior\n";}
        // Sample parameters from their prior
        preproposal_params = generateParamsPrior(Nsim);
        param_matrix = Eigen::MatrixXd::Zero(Npart, preproposal_params.cols());
        run_simulations(preproposal_params, sim_atom, &preproposal_results, &results_complete);

        std::vector<size_t> currentIndex = sort_indexes_eigen(preproposal_results); 
        for (i = 0; i < param_matrix.rows(); i++)
        {
            param_matrix.row(i) = preproposal_params.row(currentIndex[i]);
            results_double.row(i) = preproposal_results.row(currentIndex[i]); 
        }
        e0 = results_double.maxCoeff() + 1.0;
    }
    else
    {
        if (verbose > 1){Rcpp::Rcout << "Starting parameters provided\n";}
        e0 = init_eps;
        w0 = init_weights;
        param_matrix = init_param_matrix;
        results_double = init_results_double;
    }
    // Calculate current parameter SD's
    Eigen::VectorXd tau = (param_matrix.rowwise() - 
                      (param_matrix.colwise()).mean()
                      ).colwise().norm()/std::sqrt((double) 
                        (param_matrix.rows())-1.0);


    for (iteration = 0; iteration < num_iterations && !terminate; iteration++)
    {   
        // Todo: figure out how to return results even if user interrupt
        Rcpp::checkUserInterrupt();       
        if (verbose > 0){Rcpp::Rcout << "Iteration " << iteration << " [" << 
            results_double.minCoeff() << ", " << results_double.maxCoeff() <<
                "]  eps: " << e0 << "\n";
        }

        // Calculating tau

        tau = 2.0*(param_matrix.rowwise() - 
              (param_matrix.colwise()).mean()
                ).colwise().norm()/std::sqrt((double) 
                        (param_matrix.rows())-1.0);

        if (verbose > 2)
        {
            Rcpp::Rcout << "tau inverse: \n" << tau << "\n";
        }


        e1 = (samplingControlInstance -> shrinkage)*e0;

        if (verbose > 2)
        {
            Rcpp::Rcout << "   e1 = " << e1 << "\n";
            Rcpp::Rcout << "   w0, 1-10: ";
            for (i = 0; i < std::min(10, (int) w1.size()); i++)
            {
                Rcpp::Rcout << w0(i) << ", ";
            }
            Rcpp::Rcout << "\n";
        }

        // Resample Npart particles 
        // Compute cumulative weights
        cum_weights(0) = w1(0);
        for (i = 1; i < w1.size(); i++)
        {
            cum_weights(i) = w1(i) + cum_weights(i-1);
        }
        if (std::abs(cum_weights.maxCoeff() - 1) > 1e-10)
        {
            Rcpp::Rcout << "cumulative weight: " << cum_weights.maxCoeff() << "\n";
            Rcpp::stop("particle weights do not sum to one\n");
        }

        // Propose params and run simulations
        int currentIdx = 0;
        int nBatches = 0;
        while (currentIdx < Npart && 
               nBatches < maxBatches)
        {

            // Fill in param_matrix with resamples
            for (i = 0; i < Nsim; i++)
            {
               drw = U(*generator);
               for (j = 0; j < Npart; j++)
               {
                    if (drw <= cum_weights(j)) 
                    {
                        preproposal_params.row(i) = param_matrix.row(j); 
                        preproposal_results.row(i) = results_double.row(j);
                        break;
                    }
               }
            }

            // perturb parameters
            if (samplingControlInstance -> multivariatePerturbation)
            {
                proposeParams_beaumont_multivariate(&preproposal_params, 
                              &tau,
                              generator,
                              this);
            }
            else
            {
                proposeParams_beaumont(&preproposal_params, 
                              &tau,
                              generator,
                              this);     
            }

            // run simulations
            run_simulations(preproposal_params,
                            sim_atom,
                            &preproposal_results, 
                            &results_complete);

           for (i = 0; i < Nsim && currentIdx < Npart; i++)
           {
               if (preproposal_results(i,0) < e1)
               {
                   proposed_param_matrix.row(currentIdx) = 
                       preproposal_params.row(i);
                   proposed_results_double.row(currentIdx) = 
                       preproposal_results.row(i);
                   currentIdx++;
               }
           }
           if (currentIdx < Npart && verbose > 1)
           {
                Rcpp::Rcout << "  batch " << nBatches << ", " << currentIdx << 
                    "/" << Npart << " accepted\n";
           }
           nBatches ++;
        }

        w0 = w1;
        e0 = e1;
        double wtTot = 0.0;
        if (currentIdx + 1 < Npart)
        {
            if (verbose > 1)
            {
                Rcpp::Rcout << "\n";
                Rcpp::Rcout << "Maximum batches exceeded: " << currentIdx + 1 << "/" 
                    << Npart << " acceptances in " << nBatches << " batches of max " <<
                    maxBatches << "\n";
                Rcpp::Rcout << "Returning last params\n";
            }
            terminate = 1;
        }
        else
        {
            wtTot = 0.0;
            double newWt, tmpWt;
            double tmpWeightComp;
            double tmpWeight;
            Eigen::VectorXd tmpParams;
            if (samplingControlInstance -> multivariatePerturbation)
            {
                Eigen::RowVectorXd parameterMeans = param_matrix.colwise().mean();
                //auto parameterCentered = (param_matrix.rowwise() - parameterMeans);
                auto proposedParameterCentered = (proposed_param_matrix.rowwise() - parameterMeans);
                //Eigen::MatrixXd parameterCov = (parameterCentered).transpose() 
                //        * (parameterCentered)/(param_matrix.rows()-1);
                // Add a SD to overdisperse particles:
                //for (i = 0; i < parameterCov.cols(); i++) 
                //{
                //    parameterCov(i,i) += 0.5*tau(i);
                //}
                //Eigen::LLT<Eigen::MatrixXd> paramLLT = parameterCov.llt();
                //auto L = Eigen::MatrixXd(paramLLT.matrixL());
                //auto parameterICov = paramLLT.solve(Eigen::MatrixXd::Identity(L.rows(), L.cols()));
                //auto parameterICovDet = 1.0/std::pow(L.diagonal().prod(), 2.0); 

                tmpWeightComp = (-nParams/2.0)*std::log(2.0*3.14159) 
                    - 0.5*std::log((parameterICovDet));
                for (i = 0; i < w1.size(); i++)
                {
                    w1(i) = 0.0;
                    tmpWeight = tmpWeightComp;
                    tmpWeight -= (0.5*
                            (proposedParameterCentered.row(i)
                            *parameterICov*proposedParameterCentered.transpose()
                                .col(i))(0));
                    for (j = 0; j < proposedParameterCentered.rows(); j++) 
                    {
                        w1(i) += w0(j)*std::exp(tmpWeight);
                    }
                    w1(i) = evalPrior(proposed_param_matrix.row(i))
                        /w1(i);
                    wtTot += w1(i);
                }
                w1.array() /= wtTot;
            }
            else
            {
                for (i = 0; i < w1.size(); i++)
                {
                   newWt = 0.0; 
                   tmpWt = 0.0;
                   for (j = 0; j < param_matrix.rows(); j++)
                   {
                       tmpWt = 0.0;
                       for (k = 0; k < proposed_param_matrix.cols(); k++)
                       {
                          // TODO: need multivariate perturbation
                          tmpWt += R::dnorm(proposed_param_matrix(i,k),
                                            param_matrix(j,k),
                                            tau(k), 1);
                       }
                       newWt += w0(j)*std::exp(tmpWt);
                   }

                   w1(i) = evalPrior(proposed_param_matrix.row(i))/newWt;    
                   if (std::isnan(w1(i)))
                   {
                       Rcpp::stop("nan weights encountered.");
                   }
                   wtTot += w1(i);
                }
                w1.array() /= wtTot;
            }
        }

        param_matrix = proposed_param_matrix;
        results_double = proposed_results_double;
    }

    // Todo: keep an eye on this object handling. It may have unreasonable
    // overhead, and is kind of complex.  
    Rcpp::List outList;
    if (sim_type_atom == sim_result_atom)
    {
        // keep_samples indicates a debug mode, so don't worry if we can't make
        // a regular data frame from the list.
        for (i = 0; i < (int) results_complete.size(); i++)
        {
            Rcpp::List subList;
            subList["S"] = Rcpp::wrap(results_complete[i].S);
            subList["E"] = Rcpp::wrap(results_complete[i].E);
            subList["I"] = Rcpp::wrap(results_complete[i].I);
            subList["R"] = Rcpp::wrap(results_complete[i].R);

            subList["S_star"] = Rcpp::wrap(results_complete[i].S_star);
            subList["E_star"] = Rcpp::wrap(results_complete[i].E_star);
            subList["I_star"] = Rcpp::wrap(results_complete[i].I_star);
            subList["R_star"] = Rcpp::wrap(results_complete[i].R_star);
            subList["p_se"] = Rcpp::wrap(results_complete[i].p_se);
            // We p_ei and p_ir not generally defined in non-exponential case.  
            if (transitionMode == "exponential")
            {
                subList["p_ei"] = Rcpp::wrap(results_complete[i].p_ei);
                subList["p_ir"] = Rcpp::wrap(results_complete[i].p_ir);
            }
            if (hasSpatial)
            {
                subList["rho"] = Rcpp::wrap(results_complete[i].rho);
            }
            subList["beta"] = Rcpp::wrap(results_complete[i].beta);
            subList["X"] = Rcpp::wrap(results_complete[i].X);
            if (hasReinfection)
            {
                // TODO: output reinfection info
            }
            subList["result"] = Rcpp::wrap(results_complete[i].result);
            outList[std::to_string(i)] = subList;
        }
    }
    else if (sim_type_atom == sim_atom)
    {
        outList["result"] = Rcpp::wrap(results_double);
    }
    outList["params"] = Rcpp::wrap(param_matrix);
    outList["completedEpochs"] = iteration;
    outList["weights"] = Rcpp::wrap(w1);
    outList["currentEps"] = e1;
    return(outList);
}
