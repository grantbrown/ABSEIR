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


                                                                                
std::vector<size_t> sort_indexes(std::vector<int> inVec)
{
    vector<size_t> idx(inVec.size());
    for (size_t i = 0; i < idx.size(); i++)
    {
        idx[i] = i;
    }
    std::sort(idx.begin(), idx.end(),
         [&inVec](size_t i1, size_t i2){return(inVec[i1] < inVec[i2]);});
    return(idx);    
}

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

std::vector<size_t> sort_indexes_eigen_vec(Eigen::VectorXd inVec)                     
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
 
   

void proposeParams_beaumont_multivariate(Eigen::MatrixXd* outParams,
                                         Eigen::MatrixXd* inParams,
                                         Eigen::VectorXd* cum_weights,
                                         Eigen::VectorXd* tau,
                                         std::mt19937* generator,
                                         spatialSEIRModel* model)
{
    Rcpp::stop("Multivariate proposals are depricated")
}


void proposeParams_beaumont(Eigen::MatrixXd* outParams,
                            Eigen::MatrixXd* inParams,
                            Eigen::VectorXd* cum_weights,
                            Eigen::VectorXd* tau,
                            Eigen::VectorXi fixed,
                            std::mt19937* generator,
                            spatialSEIRModel* model)
{
    Rcpp::Rcout << "PROPOSING BEAUMONT\n";
    Rcpp::Rcout << "Dims: (" << outParams -> rows() << ", " << outParams -> cols() << ")\n"; 
    // Check if tau is valid
    int i,j;
    double drw = 0.0;
    auto U = std::uniform_real_distribution<double>(0,1);
    // Propose new parameters
    int p = (outParams -> cols());
    int N = outParams -> rows();
    // Check whether any tau's have gotten degenerate
    bool hasValid;
    int itrs = 0;
    for (i = 0; i < N; i++)
    {
        itrs = 0;
        hasValid = false;
        for (itrs = 0; itrs < 1000 &&  !hasValid; itrs++)
        { 
            drw = U(*generator);
            for (j = 0; j < N; j++)
            {
                if (drw <= (*cum_weights)(j)) 
                {
                    //Rcpp::Rcout << "Selecting " << (*inParams)(j,2) << " at " 
                    //<< j << ", " << drw << ", " << (*cum_weights)(j) << "\n";
                    outParams -> row(i) = inParams -> row(j); 
                    break;
                }
            }
            for (j = 0; j < p; j++)
            {
                if (!fixed(j)){
                    (*outParams)(i,j) += std::normal_distribution<double>(0.0, 
                        (*tau)(j))(*generator);
                }
            }
            hasValid = (model -> evalPrior(outParams -> row(i))) > 0;
        }
        if (!hasValid)
        {
            Rcpp::Rcout << "Unable to generate parameters with nonzero probability.\n";
            Rcpp::Rcout << "  Param " << i << " of " << outParams -> rows() << "\n"; 
            Rcpp::Rcout << "  Pror prob: " << (model -> evalPrior(outParams -> row(i)));
            Rcpp::Rcout << "  Param: \n" << outParams -> row(i) << "\n";
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
    proposed_results_complete = std::vector<simulationResultSet>();
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
    preproposal_params = Eigen::MatrixXd::Zero(samplingControlInstance -> init_batch_size, 
                                               nParams);
    preproposal_results = Eigen::MatrixXd::Zero(samplingControlInstance -> init_batch_size, 
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
    std::vector<size_t> reweight_idx;

    int i,j,k;
    int iteration;

    if (verbose > 1)
    {
        Rcpp::Rcout << "Number of iterations requested: " 
                    << num_iterations << "\n";

        dataModelInstance -> summary();
        exposureModelInstance -> summary();
        transitionPriorsInstance -> summary();
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
        /*
        preproposal_params = generateParamsPrior(Nsim);
        param_matrix = Eigen::MatrixXd::Zero(Npart, preproposal_params.cols());
        */
        preproposal_params = generateParamsPrior(samplingControlInstance -> init_batch_size);
        param_matrix = Eigen::MatrixXd::Zero(Npart, preproposal_params.cols());

        run_simulations(preproposal_params, sim_atom, &preproposal_results, &results_complete);

        std::vector<size_t> currentIndex = sort_indexes_eigen(preproposal_results); 
        for (i = 0; i < param_matrix.rows(); i++)
        {
            param_matrix.row(i) = preproposal_params.row(currentIndex[i]);
            results_double.row(i) = preproposal_results.row(currentIndex[i]); 
        }
        preproposal_params = Eigen::MatrixXd::Zero(Nsim, preproposal_params.cols());
        preproposal_results = Eigen::MatrixXd::Zero(Nsim, 
                                                samplingControlInstance -> m); 
 

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
    Eigen::VectorXd tau = 2*(param_matrix.rowwise() - 
                      (param_matrix.colwise()).mean()
                      ).colwise().norm()/std::sqrt((double) 
                        (param_matrix.rows())-1.0);

    // Determine fixed parameters
    Eigen::VectorXi fixed = Eigen::VectorXi::Zero(tau.size());
    int sz = (initialValueContainerInstance -> S0).size();
    if (initialValueContainerInstance -> type == 1){
        // In this case, we have some constant parameters. 
        for (i = tau.size()-1; i >= (tau.size()-sz*4); i--){
            fixed(i) = 1;
        }
    }

    for (iteration = 0; iteration < num_iterations && !terminate; iteration++)
    {   
        // Todo: figure out how to return results even if user interrupt
        Rcpp::checkUserInterrupt();       
        if (verbose > 0){Rcpp::Rcout << "Iteration " << iteration << " [" << 
            results_double.minCoeff() << ", " << results_double.maxCoeff() <<
                "]  eps: " << e0 << "\n";
        }

        // Calculating tau

        tau = 2*(param_matrix.rowwise() - 
              (param_matrix.colwise()).mean()
                ).colwise().norm()/std::sqrt((double) 
                        (param_matrix.rows())-1.0);
        for (i = 0; i < tau.size(); i++){
            if ((tau)(i) == 0){
                // Is this expected?
                if (!(i >= tau.size()- (initialValueContainerInstance -> S0.size())*4)){
                    Rcpp::warning("Degenerate particles detected!");
                    (tau)(i) = 0.1;
                } else if (initialValueContainerInstance -> type == 1){
                    // No pasa nada, tudo bem
                } else {
                    // Not worth warning about, discrete parameters be like that
                    (tau)(i) = 1; 
                }
            }
        }
 

        //for (i = 0; i < tau.size(); i++){
        //    tau(i) *= 1.414214;
        //}
        //Rcpp::Rcout << "Tau: " << tau << "\n";

        if (verbose > 2)
        {
            Rcpp::Rcout << "tau inverse: \n" << tau << "\n";
            Rcpp::Rcout << "...tau.size(): \n" << tau.size() << "\n"; 
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



        // Reorder parameters by weight
        reweight_idx = sort_indexes_eigen_vec(w0); 
        for (i = w0.size()-1; i >= 0; i--){
            Rcpp::Rcout << "i=" << i <<",reweight_idx[i]=" << reweight_idx[i] << "\n";
            w1(i) = w0(reweight_idx[i]);
            Rcpp::Rcout << "w1(i)=" << w1(i) << "\n";
            Rcpp::Rcout << "preproposal_params.row(i): " << preproposal_params.row(i) << "\n"; 
            Rcpp::Rcout << "param_matrix.row(reweight_idx[i]): " <<  param_matrix.row(reweight_idx[i]) << "\n";
            preproposal_params.row(i) = param_matrix.row(reweight_idx[i]);
        }

        for (i = 0; i < param_matrix.rows(); i++)
        {
            w0(i) = w1(i);
            param_matrix.row(i) = preproposal_params.row(i);
        }


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

            // perturb parameters
            if (samplingControlInstance -> multivariatePerturbation)
            {
                Rcpp::stop("Multivariate proposals are depricated")
            }
            else
            {
                proposeParams_beaumont(&preproposal_params, 
                              &param_matrix,
                              &cum_weights,
                              &tau,
                              fixed,
                              generator,
                              this);     
            }


            // run simulations
            run_simulations(preproposal_params,
                            sim_atom,
                            &preproposal_results, 
                            &results_complete);


           //std::vector<size_t> preproposal_order = sort_indexes_eigen(preproposal_results); 
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

        e0 = e1;
        w0 = w1;
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
            Eigen::VectorXd tmpParams;
           
            for (i = 0; i < proposed_param_matrix.rows(); i++)
            {
               newWt = 0.0; 
               tmpWt = 0.0;
               for (j = 0; j < param_matrix.rows(); j++)
               {
                   tmpWt = 0.0;
                   for (k = 0; k < proposed_param_matrix.cols(); k++)
                   {
                      if (!fixed(k)){
                          tmpWt += R::dnorm(proposed_param_matrix(i,k),
                                        param_matrix(j,k),
                                        tau(k), 1);
                      }
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


        w0 = w1;
        param_matrix = proposed_param_matrix;
        results_double = proposed_results_double;
    }

    Rcpp::List outList;
    if (sim_type_atom == sim_result_atom)
    {
        // Need to do an extra iteration to generate compartment data.  
        
        // NOTE: this code is largely copied from the main loop. 
        // If revisions need to be made, make sure to check both spots
        //
        // Todo: think about a way to refactor this
        
        Rcpp::checkUserInterrupt();       
        if (verbose > 0){
            Rcpp::Rcout << "Running additional iteration to capture compartments";
        }

        // Calculating tau
        tau = 2*(param_matrix.rowwise() - 
              (param_matrix.colwise()).mean()
                ).colwise().norm()/std::sqrt((double) 
                        (param_matrix.rows())-1.0);

        for (i = 0; i < tau.size(); i++){
            if ((tau)(i) == 0){
                // Is this expected?
                if (!(i >= tau.size()- (initialValueContainerInstance -> S0.size())*4)){
                    Rcpp::warning("Degenerate particles detected!");
                    (tau)(i) = 0.1;
                } else if (initialValueContainerInstance -> type == 1){
                    // No pasa nada, tudo bem
                } else {
                    // Not worth warning about, discrete parameters be like that
                    Rcpp::Rcout << "error!\n";
                    (tau)(i) = 1; 
                }
            }
        }
    



        if (verbose > 2)
        {
            Rcpp::Rcout << "tau inverse: \n" << tau << "\n";
            Rcpp::Rcout << "tau.size(): \n" << tau.size() << "\n"; 
        }


        // Back off the shrinkage
        e1 = e0/std::pow((samplingControlInstance -> shrinkage), 2);
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

        // Reorder parameters by weight
        reweight_idx = sort_indexes_eigen_vec(w0); 
        for (i = w0.size()-1; i >= 0; i--){
            w1(i) = w0(reweight_idx[i]);
            preproposal_params.row(i) = param_matrix.row(reweight_idx[i]);
        }
        for (i = 0; i < param_matrix.rows(); i++)
        {
            w0(i) = w1(i);
            param_matrix.row(i) = preproposal_params.row(i);
        }

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
        results_complete.clear();
        while (currentIdx < Npart)
        {
            // perturb parameters
            if (samplingControlInstance -> multivariatePerturbation)
            {
                proposeParams_beaumont_multivariate(&preproposal_params, 
                              &param_matrix,
                              &cum_weights,
                              &tau,
                              generator,
                              this);
            }
            else
            {
                proposeParams_beaumont(&preproposal_params, 
                              &param_matrix,
                              &cum_weights,
                              &tau,
                              fixed,
                              generator,
                              this);     
            }

           proposed_results_complete.clear();
           result_idx.clear();

            // run simulations
            run_simulations(preproposal_params,
                            sim_type_atom,
                            &preproposal_results, 
                            &proposed_results_complete);

           std::vector<size_t> result_order = sort_indexes(result_idx); 
           for (i = 0; i < Nsim && currentIdx < Npart; i++)
           {
               if (preproposal_results(i,0) < e1)
               {
                   proposed_param_matrix.row(currentIdx) = 
                       preproposal_params.row(i);
                   proposed_results_double.row(currentIdx) = 
                       preproposal_results.row(i);
                   results_complete.push_back(
                           proposed_results_complete[result_order[i]]);
                   currentIdx++;
               }
           }
           if (currentIdx < Npart && verbose > 1)
           {
                Rcpp::Rcout << "  batch " << nBatches << ", " << currentIdx << 
                    "/" << Npart << " accepted\n";
           }
           nBatches ++;
           // Need to use indexes for results complete
        }

        e0 = e1;
        w0 = w1;
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
            Eigen::VectorXd tmpParams;
            if (samplingControlInstance -> multivariatePerturbation)
            {
                // Column vector of difference
                tmpWeightComp = (-nParams/2.0)*std::log(2.0*3.14159) 
                    - 0.5*std::log((parameterICovDet));
                Eigen::MatrixXd a = Eigen::MatrixXd::Zero(param_matrix.cols(), 1);
                for (i = 0; i < proposed_param_matrix.rows(); i++)
                {
                   newWt = 0.0; 
                   for (j = 0; j < param_matrix.rows(); j++)
                   {
                       a.col(0) = (proposed_param_matrix.row(i) - 
                                   param_matrix.row(j));
                       newWt += w0(j)*std::exp(tmpWeightComp - 
                               0.5*((a.transpose()*parameterICov*a)(0)));
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
            else
            {
                for (i = 0; i < proposed_param_matrix.rows(); i++)
                {
                   newWt = 0.0; 
                   tmpWt = 0.0;
                   for (j = 0; j < param_matrix.rows(); j++)
                   {
                       tmpWt = 0.0;
                       for (k = 0; k < proposed_param_matrix.cols(); k++)
                       {
                          if (!fixed(k)){
                              tmpWt += R::dnorm(proposed_param_matrix(i,k),
                                            param_matrix(j,k),
                                            tau(k), 1);
                          }
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

        w0 = w1;
        param_matrix = proposed_param_matrix;
        results_double = proposed_results_double;
    
        // Todo: keep an eye on this object handling. It may have unreasonable
        // overhead, and is kind of complex.  
        Rcpp::List simulationResults; 
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
            simulationResults[std::to_string(i)] = subList;
        }
        outList["simulationResults"] = simulationResults;
    }
    outList["result"] = Rcpp::wrap(results_double);
    outList["params"] = Rcpp::wrap(param_matrix);
    outList["completedEpochs"] = iteration;
    outList["weights"] = Rcpp::wrap(w1);
    outList["currentEps"] = e1;
    return(outList);
}
