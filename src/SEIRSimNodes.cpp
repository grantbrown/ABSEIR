#include <iostream>
#include <sstream>
#include <random>
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <util.hpp>
#include "SEIRSimNodes.hpp"
#include "spatialSEIRModel.hpp"
#include <chrono>
#include <thread>
using namespace std;

#ifdef SPATIALSEIR_SINGLETHREAD

void printDMatrix(Eigen::MatrixXd inMat, std::string name)
{
    Rcpp::Rcout << "Matrix (" << inMat.rows() << ", " << inMat.cols() << "): " << name << "\n";
    Rcpp::Rcout << "[";
    int i,j;
    for (i = 0; i < inMat.rows(); i++)
    {
        Rcpp::Rcout << "[";
        for (j = 0; j < inMat.cols(); j++)
        {
            Rcpp::Rcout << inMat(i,j);
            if (j + 1 != inMat.cols())
            {
                Rcpp::Rcout << ", ";
            }
            else
            {
                Rcpp::Rcout << "]";
            }
        }
        if (i + 1 != inMat.rows())
        {
            Rcpp::Rcout << "\n";
        }
        else Rcpp::Rcout << "]\n";
    } 
}

void printDVector(Eigen::VectorXd inVec, std::string name)
{
    Rcpp::Rcout << "Vector (" << inVec.size() << "): " << name << "\n";
    int i;
    Rcpp::Rcout << "[";
    for (i = 0; i < inVec.size(); i++)
    {
        Rcpp::Rcout << inVec(i);
        if (i + 1 != inVec.size())
        {
            Rcpp::Rcout << ", ";
        }
        else Rcpp::Rcout << "]\n";
    }

}

void printIMatrix(Eigen::MatrixXi inMat, std::string name)
{
    Rcpp::Rcout << "Matrix (" << inMat.rows() << ", " << inMat.cols() << "): " << name << "\n";
    Rcpp::Rcout << "[";
    int i,j;
    for (i = 0; i < inMat.rows(); i++)
    {
        Rcpp::Rcout << "[";
        for (j = 0; j < inMat.cols(); j++)
        {
            Rcpp::Rcout << inMat(i,j);
            if (j + 1 != inMat.cols())
            {
                Rcpp::Rcout << ", ";
            }
            else
            {
                Rcpp::Rcout << "]";
            }
        }
        if (i + 1 != inMat.rows())
        {
            Rcpp::Rcout << "\n";
        }
        else Rcpp::Rcout << "]\n";
    } 
}

void printIVector(Eigen::VectorXi inVec, std::string name)
{
    Rcpp::Rcout << "Vector (" << inVec.size() << "): " << name << "\n";
    int i;
    Rcpp::Rcout << "[";
    for (i = 0; i < inVec.size(); i++)
    {
        Rcpp::Rcout << inVec(i);
        if (i + 1 != inVec.size())
        {
            Rcpp::Rcout << ", ";
        }
        else Rcpp::Rcout << "]\n";
    }

}
#endif
#ifndef SPATIALSEIR_SINGLETHREAD
// Dummy functions
void printDMatrix(Eigen::MatrixXd inMat, std::string name)
{
}

void printDVector(Eigen::VectorXd inVec, std::string name)
{
}

void printIMatrix(Eigen::MatrixXi inMat, std::string name)
{
}

void printIVector(Eigen::VectorXi inVec, std::string name)
{
}

#endif


NodeWorker::NodeWorker(NodePool* pl,
                       int sd,
                       Eigen::VectorXi s,
                       Eigen::VectorXi e,
                       Eigen::VectorXi i,
                       Eigen::VectorXi r,
                       Eigen::VectorXd offs,
                       Eigen::MatrixXi y,
                       MatrixXb nm,
                       std::vector<Eigen::MatrixXd> dmv,
                       std::vector<std::vector<Eigen::MatrixXd> > tdmv,
                       std::vector<int> tdme,
                       Eigen::MatrixXd x,
                       Eigen::MatrixXd x_rs,
                       std::string mode,
                       Eigen::MatrixXd ei_prior,
                       Eigen::MatrixXd ir_prior,
                       double avgI,
                       Eigen::VectorXd sp_prior,
                       Eigen::VectorXd se_prec,
                       Eigen::VectorXd rs_prec,
                       Eigen::VectorXd se_mean,
                       Eigen::VectorXd rs_mean,
                       double ph,
                       int dmc,
                       bool cmltv,
                       int m)
{
    pool = pl;
    node = std::unique_ptr<SEIR_sim_node>(new SEIR_sim_node(this, sd,s,e,i,
                         r,offs,y,nm,dmv,tdmv,tdme,x,x_rs,mode,ei_prior,ir_prior,avgI,
                         sp_prior,se_prec,rs_prec,se_mean,rs_mean, ph,dmc,cmltv, m));
}

void NodeWorker::operator()()
{
    instruction task;
#ifdef SPATIALSEIR_SINGLETHREAD
    while ((pool -> tasks).size() > 0)
    {
        task = (pool -> tasks).front();
        (pool -> nBusy)++;
        (pool -> tasks).pop_front();
        if (task.action_type == sim_atom)
        {
            Eigen::VectorXd result = node -> simulate(task.params, false).result;
            (*(pool -> result_pointer)).row(task.param_idx) = result; 
        }
        else if (task.action_type == sim_result_atom)
        {
            // Do these need to be re-sorted?
            simulationResultSet result = node -> simulate(task.params, true);
            pool -> result_complete_pointer -> push_back(result);
            pool -> index_pointer -> push_back(task.param_idx);
        }
        (pool -> nBusy)--;
    }
#else
    while(true)
    {
        {
            std::unique_lock<std::mutex> lock(pool -> queue_mutex);

            while (!pool -> exit && (pool -> tasks).empty())
            {
                (pool -> condition).wait(lock);
            }
            if (pool -> exit)
                return;
            task = (pool -> tasks).front();
            (pool -> nBusy)++;
            (pool -> tasks).pop_front();
        }
        if (task.action_type == sim_atom)
        {
            Eigen::VectorXd result = node -> simulate(task.params, false).result;
            {
                std::unique_lock<std::mutex> lock(pool -> result_mutex);
                (*(pool -> result_pointer)).row(task.param_idx) = result; 
                while (!(node -> messages).empty()) 
                {
                    (pool -> messages).push_back((node -> messages).front()); 
                    (node -> messages).pop_front();
                }
            }
        }
        else if (task.action_type == sim_result_atom)
        {
            simulationResultSet result = node -> simulate(task.params, true);
            {
                std::unique_lock<std::mutex> lock(pool -> result_mutex);
                pool -> index_pointer -> push_back(task.param_idx);
                pool -> result_complete_pointer -> push_back(result);
                while (!((node -> messages).empty())) 
                {
                    (pool -> messages).push_back((node -> messages).front()); 
                    (node -> messages).pop_front();
                }
            }
        }

        {
            std::unique_lock<std::mutex> lock(pool -> queue_mutex);
            (pool -> nBusy)--;
            (pool -> finished).notify_one();
        }
    }
#endif
}

NodePool::NodePool(Eigen::MatrixXd* rslt_ptr,
                   std::vector<simulationResultSet>* rslt_c_ptr,
                   std::vector<int>* idx_ptr,
                       int threads,
                       int sd,
                       Eigen::VectorXi s,
                       Eigen::VectorXi e,
                       Eigen::VectorXi i,
                       Eigen::VectorXi r,
                       Eigen::VectorXd offs,
                       Eigen::MatrixXi y,
                       MatrixXb nm,
                       std::vector<Eigen::MatrixXd> dmv,
                       std::vector<std::vector<Eigen::MatrixXd> > tdmv,
                       std::vector<int> tdme,
                       Eigen::MatrixXd x,
                       Eigen::MatrixXd x_rs,
                       std::string mode,
                       Eigen::MatrixXd ei_prior,
                       Eigen::MatrixXd ir_prior,
                       double avgI,
                       Eigen::VectorXd sp_prior,
                       Eigen::VectorXd se_prec,
                       Eigen::VectorXd rs_prec,
                       Eigen::VectorXd se_mean,
                       Eigen::VectorXd rs_mean,
                       double ph,
                       int dmc,
                       bool cmltv,
                       int m)
{
    result_pointer = rslt_ptr;
    result_complete_pointer = rslt_c_ptr;
    index_pointer = idx_ptr;
    exit = false;
    nBusy = 0;
#ifdef SPATIALSEIR_SINGLETHREAD
    // Single threaded mode only needs single worker
    nodes.push_back(NodeWorker(this,
                                           sd + 1000*(1),s,e,i,
                     r,offs,y,nm,dmv,tdmv,tdme,x,x_rs,mode,ei_prior,ir_prior,avgI,
                     sp_prior,se_prec,rs_prec,se_mean,rs_mean,ph,dmc,cmltv, m
                    ));
#else
    for (int itr = 0; itr < threads; itr++)
    {
        nodes.push_back(std::thread(NodeWorker(this,
                                               sd + 1000*(itr+1),s,e,i,
                         r,offs,y,nm,dmv,tdmv,tdme,x,x_rs,mode,ei_prior,ir_prior,avgI,
                         sp_prior,se_prec,rs_prec,se_mean,rs_mean,ph,dmc,cmltv, m
                        )));
    }
#endif
}

void NodePool::setResultsDest(Eigen::MatrixXd* rslt_ptr,
                              std::vector<simulationResultSet>* rslt_c_ptr)
{
    result_pointer = rslt_ptr;
    result_complete_pointer = rslt_c_ptr;
}


void NodePool::awaitFinished()
{
#ifdef SPATIALSEIR_SINGLETHREAD
    nodes[0]();
#else
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        finished.wait(lock, [this](){ resolveMessages();
                return tasks.empty() && (nBusy == 0); });
    }
#endif
}

void NodePool::resolveMessages()
{
    {
        std::unique_lock<std::mutex> lock(result_mutex);
        while (!(messages.empty())) 
        {
            Rcpp::Rcout << messages.front() << "\n"; 
            messages.pop_front();
        }
    }
}

void NodePool::enqueue(std::string action_type, int param_idx, Eigen::VectorXd params)
{
    instruction inst;
    inst.param_idx = param_idx;
    inst.action_type = action_type;
    inst.params = params;
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        tasks.push_back(inst);
    }
    resolveMessages();
    condition.notify_one();
}

NodePool::~NodePool()
{
    exit = true;
    condition.notify_all();
    for (unsigned int i = 0; i < nodes.size(); i++)
    {
#ifdef SPATIALSEIR_SINGLETHREAD
        nodes.pop_back();
#else
        nodes[i].join();
#endif
    }
}


SEIR_sim_node::SEIR_sim_node(NodeWorker* worker,
                             int sd,
                             Eigen::VectorXi s,
                             Eigen::VectorXi e,
                             Eigen::VectorXi i,
                             Eigen::VectorXi r,
                             Eigen::VectorXd offs,
                             Eigen::MatrixXi y,
                             MatrixXb nm,
                             std::vector<Eigen::MatrixXd> dmv,
                             std::vector<std::vector<Eigen::MatrixXd> > tdmv,
                             std::vector<int> tdme,
                             Eigen::MatrixXd x,
                             Eigen::MatrixXd x_rs,
                             std::string mode,
                             Eigen::MatrixXd ei_prior,
                             Eigen::MatrixXd ir_prior,
                             double avgI,
                             Eigen::VectorXd sp_prior,
                             Eigen::VectorXd se_prec,
                             Eigen::VectorXd rs_prec,
                             Eigen::VectorXd se_mean,
                             Eigen::VectorXd rs_mean,
                             double ph,
                             int dmc,
                             bool cmltv,
                             int m_
                             ) : parent(worker),
                                 random_seed(sd),
                                 S0(s),
                                 E0(e),
                                 I0(i),
                                 R0(r),
                                 offset(offs),
                                 Y(y),
                                 na_mask(nm),
                                 DM_vec(dmv),
                                 TDM_vec(tdmv),
                                 TDM_empty(tdme),
                                 X(x),
                                 X_rs(x_rs),
                                 transitionMode(mode),
                                 E_to_I_prior(ei_prior),
                                 I_to_R_prior(ir_prior),
                                 inf_mean(avgI),
                                 spatial_prior(sp_prior),
                                 exposure_precision(se_prec),
                                 reinfection_precision(rs_prec),
                                 exposure_mean(se_mean),
                                 reinfection_mean(rs_mean),
                                 phi(ph),
                                 data_compartment(dmc),
                                 cumulative(cmltv),
                                 m(m_)
{
    try
    {
        std::minstd_rand0 lc_generator(sd);
        std::uint_least32_t seed_data[std::mt19937::state_size];
        std::generate_n(seed_data, std::mt19937::state_size, std::ref(lc_generator));
        std::seed_seq q(std::begin(seed_data), std::end(seed_data));
        generator = new mt19937{q};   
        int i;

        E_paths = std::vector<Eigen::MatrixXi>();
        I_paths = std::vector<Eigen::MatrixXi>();

        if (transitionMode == "weibull")
        {
            EI_transition_dist = std::unique_ptr<weibullTransitionDistribution>(
                    new weibullTransitionDistribution(E_to_I_prior.col(0)));
            IR_transition_dist = std::unique_ptr<weibullTransitionDistribution>(
                    new weibullTransitionDistribution(I_to_R_prior.col(0)));
            for (i = 0; i < m; i++){
                E_paths.push_back(Eigen::MatrixXi((int) E_to_I_prior(4,0), Y.cols()));
                I_paths.push_back(Eigen::MatrixXi((int) I_to_R_prior(4,0), Y.cols()));
            }
        }
        else if (transitionMode == "path_specific")
        {
            for (i = 0; i < m; i++){
                E_paths.push_back(Eigen::MatrixXi(E_to_I_prior.rows(),Y.cols()));
                I_paths.push_back(Eigen::MatrixXi(I_to_R_prior.rows(),Y.cols()));
            }
        }
        else
        {
            for (i = 0; i < m; i++)
            {
                E_paths.push_back(Eigen::MatrixXi(1,Y.cols()));
                I_paths.push_back(Eigen::MatrixXi(1,Y.cols()));
            }
        }
        if (phi > 0)
        {   
            // We take the floor of the resulting continuous normal, so shift by 0.5
            overdispersion_distribution = std::normal_distribution<double>(0.5, 
                    1.0/phi);
        }
    }
    catch (int e)
    {
        // TODO: handle these errors
        //aout(this) << "Error in constructor: " << e << "\n"; 
    }
    has_reinfection = (reinfection_precision(0) > 0); 
    has_spatial = (Y.cols() > 1);
    has_ts_spatial = (TDM_vec[0].size() > 0);
    const int nRho = (has_spatial && has_ts_spatial ? DM_vec.size() + TDM_vec[0].size() :
                     (has_spatial ? DM_vec.size() : 0));
    const int nReinf = (has_reinfection ? X_rs.cols() : 0);
    const int nBeta = X.cols();
    const int nTrans = (transitionMode == "exponential" ? 2 : 
                       (transitionMode == "weibull" ? 4 : 0));
    total_size = nRho + nReinf + nBeta + nTrans;
}

simulationResultSet SEIR_sim_node::simulate(Eigen::VectorXd params, bool keepCompartments)
{
    // Params is a vector made of:
    // [Beta, Beta_RS, rho, gamma_ei, gamma_ir]    
    int time_idx, i, j, k;   
    unsigned int idx; 

    simulationResultSet compartmentResults;
 
    Eigen::VectorXd results = Eigen::VectorXd::Zero(m); 
    Eigen::VectorXd beta = params.segment(0, X.cols()); 
    Eigen::VectorXd beta_rs;

    if (has_reinfection) 
    {
        beta_rs = params.segment(X.cols(), X_rs.cols());
    }
    else 
    {
        beta_rs = Eigen::VectorXd(1);
    }

    int nRho = (has_spatial && has_ts_spatial ? DM_vec.size() + TDM_vec[0].size() :
                     (has_spatial ? DM_vec.size() : 0));
    Eigen::VectorXd rho;
    if (has_spatial && has_reinfection)
    {
        rho = params.segment(X.cols() + X_rs.cols(), nRho);
    }
    else if (has_spatial)
    {
        rho = params.segment(X.cols(), nRho);
    }
    else
    {
        rho = Eigen::VectorXd(1.0);
    }
   
    // Should really unify these two code paths...
    double gamma_ei = (transitionMode == "exponential" ? 
            params(params.size() - 2) : -1.0);
    double gamma_ir = (transitionMode == "exponential" ? 
            params(params.size() - 1) : -1.0);

    Eigen::VectorXd EI_params;
    Eigen::VectorXd IR_params;
    if (transitionMode == "weibull")
    {
        EI_params = params.segment(params.size() - 4, 2);
        IR_params = params.segment(params.size() - 2, 2); 
        EI_transition_dist -> setCurrentParams(EI_params);
        IR_transition_dist -> setCurrentParams(IR_params);
    } 
    else
    {
        EI_params = Eigen::VectorXd::Zero(1);
        IR_params = Eigen::VectorXd::Zero(1);
    }

    // Both Weibull and arbitrary path specific priors require
    // Empty paths at beginning of sim
    if (transitionMode != "exponential")
    {
        // Clear out paths
        for (k = 0; k < m; k++)
        {
            for (i = 0; i < E_paths[k].cols(); i++)
            {
                for(j = 0; j < E_paths[k].rows(); j++)
                {
                    E_paths[k](j,i) = 0;
                }
            }
            for (i = 0; i < I_paths[k].cols(); i++)
            {
                for(j = 0; j < I_paths[k].rows(); j++)
                {
                    I_paths[k](j,i) = 0;
                }
            }
        }
    }

    Eigen::VectorXi N(S0.size());
    N = (S0 + E0 + I0 + R0);

    Eigen::MatrixXd p_se_cache(S0.size(), m);

    Eigen::MatrixXi current_S(S0.size(), m);
    Eigen::MatrixXi current_E(S0.size(), m);
    Eigen::MatrixXi current_I(S0.size(), m);
    Eigen::MatrixXi current_R(S0.size(), m);

    std::vector<compartment_tap> I_lag;
    for (i = 0; i < m; i++)
    {
        I_lag.push_back(compartment_tap(TDM_vec[0].size(), S0.size()));
    }

    Eigen::MatrixXi previous_S(S0.size(), m);
    Eigen::MatrixXi previous_E(S0.size(), m);
    Eigen::MatrixXi previous_I(S0.size(), m);
    Eigen::MatrixXi previous_R(S0.size(), m);


    Eigen::MatrixXi previous_S_star(S0.size(), m);
    Eigen::MatrixXi previous_E_star(S0.size(), m);
    Eigen::MatrixXi previous_I_star(S0.size(), m);
    Eigen::MatrixXi previous_R_star(S0.size(), m); 

    Eigen::MatrixXi cumulative_compartment(S0.size(), m);

    Eigen::MatrixXi* comparison_compartment = (data_compartment == 0 ?
                                               &previous_I_star : 
                                              (data_compartment == 1 ? 
                                               &previous_R_star : 
                                              (data_compartment == 2 ? 
                                               &previous_I : &previous_I_star)));

    // Calculate probabilities
    // p_se calculation
    // Equivalent R expression: 
    // exp(matrix(X %*% beta, nrow = nrow(Y), ncol = ncol(Y)))
    Eigen::MatrixXd eta = (X*beta).unaryExpr([](double elem){return(
                std::exp(elem));
            });

    Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor> p_se_components(eta.data(), 
                Y.rows(), Y.cols());

    time_idx = 0;
    for (i = 0; i < m; i++)
    {
        previous_S.col(i) = S0;
        previous_E.col(i) = E0;
        previous_I.col(i) = I0;
        previous_R.col(i) = R0;
        previous_S_star.col(i) = Eigen::VectorXi::Zero(m);
        previous_E_star.col(i) = Eigen::VectorXi::Zero(m);
        previous_I_star.col(i) = Eigen::VectorXi::Zero(m);
        previous_R_star.col(i) = Eigen::VectorXi::Zero(m);
    }
    if (has_ts_spatial)
    {
        for (i = 0; i < m; i++)
        {
            I_lag[i].push(previous_I.col(i));
        }
    }

    if (transitionMode != "exponential")
    {
        for (i = 0; i < m; i++)
        {
            I_paths[i].row(0) = I0;
            E_paths[i].row(0) = E0;
        }
    }

    p_se_cache = ((previous_I.cast<double>().array().colwise())
        /N.cast<double>().array()).array().colwise()*
        p_se_components.row(0).transpose().array();
    printDVector(p_se_components.row(0), "p_se_components.row(0)");
    printDVector(p_se_cache, "p_se_cache");

    Eigen::MatrixXd p_se = 1*p_se_cache; 

    if (has_spatial)
    {
        for (idx = 0; idx < DM_vec.size(); idx++)
        {
        
            printDMatrix(DM_vec[idx], "DM_vec[i]");
            printDVector(DM_vec[idx]*(p_se_cache), "DM_vec[idx]*(p_se_cache)");
            p_se += rho[idx]*(DM_vec[idx] * (p_se_cache));
            printDVector(p_se, "p_se");

        }
    }
    /*
    if (has_ts_spatial)
    {
        if (!TDM_empty[0])
        {
            p_se += rho[DM_vec.size()]*(TDM_vec[0][0] * p_se_cache);
        }
    }
    */

    p_se = (((-1.0*p_se.array()) * (offset(0)))).unaryExpr([](double e){
            return(1-std::exp(e));
            }); 
    printDVector(p_se, "p_se");


    // Not used if transitionMode != "exponential"
    Eigen::VectorXd p_ei = (-1.0*gamma_ei*offset)
                            .unaryExpr([](double e){return(1-std::exp(e));});
    
    Eigen::VectorXd p_ir = (-1.0*gamma_ir*offset)
                            .unaryExpr([](double e){return(1-std::exp(e));}); 

    Eigen::VectorXd p_rs;
    if (has_reinfection)
    {
        p_rs = ((((X_rs*beta_rs).unaryExpr([](double e){return(std::exp(e));})).array() 
                    * offset(0)).matrix()).unaryExpr([](double e){return(1-std::exp(-e));});
    }
    else
    {
        p_rs = Eigen::VectorXd::Zero(Y.rows());
    }

    // Initialize debug info if applicable
    if (keepCompartments)
    {
        compartmentResults.result = Eigen::VectorXd(m);
        for (i = 0; i < m; i++)
        {
            compartmentResults.result(i) = 0.0;
        }
        compartmentResults.S = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.E = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.I = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.R = Eigen::MatrixXi(Y.rows(), Y.cols());

        compartmentResults.S_star = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.E_star = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.I_star = Eigen::MatrixXi(Y.rows(), Y.cols());
        compartmentResults.R_star = Eigen::MatrixXi(Y.rows(), Y.cols());
        
        compartmentResults.X = Eigen::MatrixXd(X.rows(), X.cols());  
        compartmentResults.X = X; 
        compartmentResults.beta = Eigen::MatrixXd(1, beta.size());
        compartmentResults.beta = beta.transpose(); 
        // Todo: add reinfection component
        compartmentResults.rEA = Eigen::MatrixXd(Y.rows(), Y.cols());
        compartmentResults.r0t = Eigen::MatrixXd(Y.rows(), Y.cols());
        compartmentResults.effR0 = Eigen::MatrixXd(Y.rows(), Y.cols());
        compartmentResults.p_se = Eigen::MatrixXd(Y.rows(), Y.cols());
        compartmentResults.p_se.row(0) = p_se.col(0);
        if (transitionMode == "exponential")
        {
            compartmentResults.p_ei = p_ei.transpose(); 
            compartmentResults.p_ir = p_ir.transpose(); 
        }
        else
        {
            compartmentResults.p_ei = EI_params;
            compartmentResults.p_ir = IR_params; 
        }
        if (has_spatial)
        {
            compartmentResults.rho = rho.transpose();
        }
        else
        {
            compartmentResults.rho = Eigen::MatrixXd(1,1);
            compartmentResults.rho(0, 0) = 0.0; 
        }
    }

    // Initialize Random Draws
    time_idx = 0;     

    int tmpDraw, w;
    for (w = 0; w < m; w++)
    {
        for (i = 0; i < Y.cols(); i++)
        {
            previous_S_star(i, w) = std::binomial_distribution<int>(
                    previous_R(i, w), p_rs(0))(*generator);
            previous_E_star(i, w) = std::binomial_distribution<int>(
                    previous_S(i, w), p_se(i, w))(*generator);

            if (transitionMode == "exponential")
            {
                previous_I_star(i, w) = std::binomial_distribution<int>(
                        previous_E(i, w), p_ei(0))(*generator);
                previous_R_star(i, w) = std::binomial_distribution<int>(
                        previous_I(i, w), p_ir(0))(*generator);
            }
            else if (transitionMode == "path_specific")
            {
                // Updating E_paths and I_paths is repetative - factor out into a function?
                // That might be costly, unless it's inlined...
                previous_I_star(i, w) = E_paths[w](E_paths[w].rows() - 1, i); // could take outside loop
                E_paths[w](E_paths[w].rows() -1, i) = 0;
                for (j = 0; j < offset(0); j++)
                {
                    // TODO: stop early when possible
                    // idea: cache previous max?
                    for (k = E_paths[w].rows() - 2; 
                            k >= 0; k--)
                    {
                        if (E_paths[w](k,i) > 0)
                        {
                            tmpDraw = std::binomial_distribution<int>(
                                    E_paths[w](k,i), E_to_I_prior(k,5))(*generator);
                            previous_I_star(i, w) += tmpDraw;
                            E_paths[w](k,i) -= tmpDraw;
                            E_paths[w](k+1, i) = E_paths[w](k,i);
                            E_paths[w](k,i) = 0; // not needed?
                        }
                    }
                }
                previous_R_star(i, w) = I_paths[w](I_paths[w].rows() - 1, i); // could take outside loop
                I_paths[w](I_paths[w].rows() -1, i) = 0;
                for (j = 0; j < offset(0); j++)
                {
                    // TODO: stop early when possible
                    // idea: cache previous max?
                    for (k = I_paths[w].rows() - 2; 
                            k >= 0; k--)
                    {
                        if (I_paths[w](k,i) > 0)
                        {
                            tmpDraw = std::binomial_distribution<int>(
                                    I_paths[w](k,i), I_to_R_prior(k,5))(*generator);
                            previous_R_star(i,w) += tmpDraw;
                            I_paths[w](k, i) -= tmpDraw;
                            I_paths[w](k+1, i) = I_paths[w](k,i);
                            I_paths[w](k, i) = 0; // not needed?
                        }
                    }
                }
            }
            else
            {
                // Updating E_paths and I_paths is repetative - factor out into a function?
                // That might be costly, unless it's inlined...
                previous_I_star(i, w) = E_paths[w](E_paths[w].rows() - 1, i); // could take outside loop
                E_paths[w](E_paths[w].rows() - 1, i) = 0;
                for (j = 0; j < offset(0); j++)
                {
                    // TODO: stop early when possible
                    // idea: cache previous max?
                    for (k = E_paths[w].rows() - 2; 
                            k >= 0; k--)
                    {
                        if (E_paths[w](k,i) > 0)
                        { 
                            tmpDraw = std::binomial_distribution<int>(
                                        E_paths[w](k,i),
                                        EI_transition_dist -> getTransitionProb(k, k+1)
                                        )(*generator);
                            previous_I_star(i,w) += tmpDraw;
                            E_paths[w](k,i) -= tmpDraw;
                            E_paths[w](k+1, i) = E_paths[w](k,i);
                            E_paths[w](k,i) = 0; // not needed?
                        }
                    }
                }
                previous_R_star(i,w) = I_paths[w](I_paths[w].rows() - 1, i); // could take outside loop
                I_paths[w](I_paths[w].rows() -1, i) = 0;
                for (j = 0; j < offset(0); j++)
                {
                    // TODO: stop early when possible
                    // idea: cache previous max?
                    for (k = I_paths[w].rows() - 2; 
                            k >= 0; k--)
                    {
                        if (I_paths[w](k,i) > 0)
                        {
                            tmpDraw = std::binomial_distribution<int>(
                                        I_paths[w](k,i),
                                        IR_transition_dist -> getTransitionProb(k, k+1)
                                        )(*generator);
                            previous_R_star(i,w) += tmpDraw;
                            I_paths[w](k,i) -= tmpDraw;
                            I_paths[w](k+1, i) = I_paths[w](k,i);
                            I_paths[w](k,i) = 0; // not needed?
                        }
                    }
                }
            }

            if (cumulative)
            {
                cumulative_compartment(i) = (*comparison_compartment)(i);
            }

            results(w) += (na_mask(0,i) ? 0 : 
                    pow(((*comparison_compartment)(i,w) + 
                            (phi > 0 ? 
                             std::floor(overdispersion_distribution(*generator)) : 0)
                            - Y(0, i)), 2.0)); 
        }// End i loop
    }// End w loop

    if (keepCompartments)
    {
        for (i = 0; i < S0.size(); i++)
        {
            compartmentResults.S(time_idx, i) = S0(i);
            compartmentResults.E(time_idx, i) = E0(i);
            compartmentResults.I(time_idx, i) = I0(i);
            compartmentResults.R(time_idx, i) = R0(i);

            compartmentResults.S_star(time_idx, i) = previous_S_star(i,0);
            compartmentResults.E_star(time_idx, i) = previous_E_star(i,0);
            compartmentResults.I_star(time_idx, i) = previous_I_star(i,0);
            compartmentResults.R_star(time_idx, i) = previous_R_star(i,0);
        }
    }

    current_S = previous_S + previous_S_star - previous_E_star;
    current_E = previous_E + previous_E_star - previous_I_star;
    current_I = previous_I + previous_I_star - previous_R_star;
    current_R = previous_R + previous_R_star - previous_S_star;

    if (transitionMode != "exponential")
    {
        for (w = 0; w < m; w++)
        {
            E_paths[w].row(0) = previous_E_star.col(w);
            I_paths[w].row(0) = previous_I_star.col(w);
        }
    }

    previous_S = current_S;
    previous_E = current_E;
    previous_I = current_I;
    previous_R = current_R;

    if (has_ts_spatial)
    {
        for (w = 0; w < m; w++)
        {
            I_lag[w].push(previous_I.col(w));
        }
    }

    // Simulation: iterative case
    int lag;
    for (w = 0; w < m; w++)
    {
        for (time_idx = 1; time_idx < Y.rows(); time_idx++)
        {

            /* 
            p_se_cache = ((previous_I.cast<double>().array().colwise())
                /N.cast<double>().array()).array().colwise()*
                p_se_components.row(0).transpose().array();


               */
            printDVector(
                (previous_I.cast<double>().array().col(w)).array(),
                " previous_I.cast<double>().array().col(w)).array()");


            p_se_cache = (previous_I.cast<double>().array().col(w)).array()
                /N.cast<double>().array()*p_se_components.row(time_idx).transpose().array();
            printDVector(p_se_cache, "p_se_cache");

            p_se = 1*p_se_cache; 
            if (has_spatial)
            {
                for (idx = 0; idx < DM_vec.size(); idx++)
                {
                    printDMatrix(DM_vec[idx], "DM_vec[idx]");
                    printDMatrix(rho[idx]*(DM_vec[idx] * p_se_cache),
                            "rho[idx]*(DM_vec[idx] * p_se_cache)");

                    p_se += rho[idx]*(DM_vec[idx] * p_se_cache);
                    printDVector(p_se, "p_se");
                }
            }

            if (has_ts_spatial && !TDM_empty[time_idx])
            {
                for (lag = 0; time_idx - lag >= 0 && lag < (int) TDM_vec[0].size(); lag++)
                {
                    p_se_cache = (I_lag[w].get(lag).cast<double>().array())
                        /N.cast<double>().array()*p_se_components.row(time_idx - lag).transpose().array();
                    p_se += rho[DM_vec.size() + lag]*(TDM_vec[time_idx-lag][lag] * p_se_cache);
                }
            }


            p_se = ((-1.0*p_se.array() * offset(time_idx)).matrix()
                    ).unaryExpr([](double e){return(1-std::exp(e));});
     
            printDMatrix(p_se, "p_se");
            for (i = 0; i < Y.cols(); i++)
            {
                previous_S_star(i,w) = std::binomial_distribution<int>(previous_R(i,w), p_rs(time_idx))(*generator);
                previous_E_star(i,w) = std::binomial_distribution<int>(previous_S(i,w), p_se(i,0))(*generator);

                if (transitionMode == "exponential")
                {
                    previous_I_star(i,w) = std::binomial_distribution<int>(previous_E(i,w), p_ei(time_idx))(*generator);
                    previous_R_star(i,w) = std::binomial_distribution<int>(previous_I(i,w), p_ir(time_idx))(*generator);
                }
                else if (transitionMode == "path_specific")
                {
                    // Updating E_paths and I_paths is repetative - factor out into a function?
                    // That might be costly, unless it's inlined...
                    previous_I_star(i,w) = E_paths[w](E_paths[w].rows() - 1, i); // could take outside loop
                    E_paths[w](E_paths[w].rows() -1, i) = 0;
                    for (j = 0; j < offset(time_idx); j++)
                    {
                        // TODO: stop early when possible
                        // idea: cache previous max?
                        for (k = E_paths[w].rows() - 2; 
                                k >= 0; k--)
                        {
                            if (E_paths[w](k,i) > 0)
                            {
                                tmpDraw = std::binomial_distribution<int>(E_paths[w](k,i), E_to_I_prior(k,5))(*generator);
                                previous_I_star(i,w) += tmpDraw;
                                E_paths[w](k,i) -= tmpDraw;
                                E_paths[w](k+1, i) = E_paths[w](k,i);
                                E_paths[w](k,i) = 0; // not needed?
                            }
                        }
                    }
                    previous_R_star(i,w) = I_paths[w](I_paths[w].rows() - 1, i); // could take outside loop
                    I_paths[w](I_paths[w].rows() -1, i) = 0;
                    for (j = 0; j < offset(time_idx); j++)
                    {
                        // TODO: stop early when possible
                        // idea: cache previous max?
                        for (k = I_paths[w].rows() - 2; 
                                k >= 0; k--)
                        {
                            if (I_paths[w](k,i) > 0)
                            {
                                tmpDraw = std::binomial_distribution<int>(I_paths[w](k,i), I_to_R_prior(k,5))(*generator);
                                previous_R_star(i,w) += tmpDraw;
                                I_paths[w](k,i) -= tmpDraw;
                                I_paths[w](k+1, i) = I_paths[w](k,i);
                                I_paths[w](k,i) = 0; // not needed?
                            }
                        }
                    }
                }
                else
                {
                    // Updating E_paths and I_paths is repetative - factor out into a function?
                    // That might be costly, unless it's inlined...
                    previous_I_star(i,w) = E_paths[w](E_paths[w].rows() - 1, i); // could take outside loop
                    E_paths[w](E_paths[w].rows() - 1, i) = 0;
                    for (j = 0; j < offset(time_idx); j++)
                    {
                        // TODO: stop early when possible
                        // idea: cache previous max?
                        for (k = E_paths[w].rows() - 2; 
                                k >= 0; k--)
                        {
                            if (E_paths[w](k,i) > 0)
                            {
                                tmpDraw = std::binomial_distribution<int>(
                                        E_paths[w](k,i), 
                                        EI_transition_dist -> getTransitionProb(k, k+1))(*generator); 
                                previous_I_star(i,w) += tmpDraw;
                                E_paths[w](k,i) -= tmpDraw;
                                E_paths[w](k+1, i) = E_paths[w](k,i);
                                E_paths[w](k,i) = 0; // not needed?
                            }
                        }
                    }
                    previous_R_star(i,w) = I_paths[w](I_paths[w].rows() - 1, i); // could take outside loop
                    I_paths[w](I_paths[w].rows() -1, i) = 0;
                    for (j = 0; j < offset(time_idx); j++)
                    {
                        // TODO: stop early when possible
                        // idea: cache previous max?
                        for (k = I_paths[w].rows() - 2; 
                                k >= 0; k--)
                        {
                            if (I_paths[w](k,i) > 0)
                            {
                                tmpDraw = std::binomial_distribution<int>(
                                        I_paths[w](k,i), 
                                        IR_transition_dist -> getTransitionProb(k,k+1)
                                        )(*generator); 
                                previous_R_star(i,w) += tmpDraw;
                                I_paths[w](k,i) -= tmpDraw;
                                I_paths[w](k+1, i) = I_paths[w](k,i);
                                I_paths[w](k,i) = 0; // not needed?
                            }
                        }
                    }

                }
                if (cumulative)
                {
                    cumulative_compartment(i,w) += (*comparison_compartment)(i,w); 

                    results(w) += (na_mask(time_idx, i) ? 0 : 
                        pow(((cumulative_compartment)(i,w) 
                                + (phi > 0 ? 
                                    std::floor(overdispersion_distribution(*generator)) 
                                    : 0)
                                - Y(time_idx, i)), 2.0)); 
                }
                else
                {
                    results(w) += (na_mask(time_idx, i) ? 0 : 
                            pow(((*comparison_compartment)(i,w) + (phi > 0 ? 
                                        std::floor(overdispersion_distribution(*generator))
                                        : 0)- Y(time_idx, i)), 2.0)); 
                }
            }

            if (keepCompartments)
            {
                for (i = 0; i < S0.size(); i++)
                {
                    compartmentResults.S(time_idx, i) = previous_S(i,0);
                    compartmentResults.E(time_idx, i) = previous_E(i,0);
                    compartmentResults.I(time_idx, i) = previous_I(i,0);
                    compartmentResults.R(time_idx, i) = previous_R(i,0);

                    compartmentResults.S_star(time_idx, i) = previous_S_star(i,0);
                    compartmentResults.E_star(time_idx, i) = previous_E_star(i,0);
                    compartmentResults.I_star(time_idx, i) = previous_I_star(i,0);
                    compartmentResults.R_star(time_idx, i) = previous_R_star(i,0);

                    compartmentResults.p_se.row(time_idx) = p_se.col(0);
                }
            }

            current_S.col(w) = previous_S.col(w) + previous_S_star.col(w) - previous_E_star.col(w);
            current_E.col(w) = previous_E.col(w) + previous_E_star.col(w) - previous_I_star.col(w);
            current_I.col(w) = previous_I.col(w) + previous_I_star.col(w) - previous_R_star.col(w);
            current_R.col(w) = previous_R.col(w) + previous_R_star.col(w) - previous_S_star.col(w);

            if (transitionMode != "exponential")
            {
                E_paths[w].row(0) = previous_E_star.col(w);
                I_paths[w].row(0) = previous_I_star.col(w);
            }

            previous_S.col(w) = current_S.col(w);
            previous_E.col(w) = current_E.col(w);
            previous_I.col(w) = current_I.col(w);
            previous_R.col(w) = current_R.col(w);
        }
    }

    if (keepCompartments)
    {
        calculateReproductiveNumbers(&compartmentResults);
    }
    compartmentResults.result = results.unaryExpr([](double elem){
            return(std::sqrt(elem));
            }); 
    return(compartmentResults);
}

void SEIR_sim_node::nodeMessage(std::string msg)
{
    messages.push_back(msg);
}

void SEIR_sim_node::calculateReproductiveNumbers(simulationResultSet* results)
{
    int i, l, k;
    Eigen::VectorXi N(S0.size());
    N = S0 + E0 + I0 + R0;

    int time_idx;
    Eigen::VectorXd beta = (*results).beta; 

    Eigen::MatrixXd eta = (((*results).X*beta).unaryExpr([](double elem){return(std::exp(elem));}));
    Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor> p_se_components(eta.data(), 
                (*results).S.rows(), (*results).S.cols());
    Eigen::MatrixXd p_se_cache(1, (*results).S.cols());

    int nLoc = (*results).S.cols();
    int nTpt = (*results).S.rows();
    std::vector<Eigen::MatrixXd> GVector;
    // Fill in R0(t) and eff R0(t)
    // I_to_R_prior(k, 4)
    if (transitionMode == "weibull")
    {
        IR_transition_dist -> setCurrentParams((*results).p_ir.col(0));
    }
    for (time_idx = 0; time_idx < nTpt; time_idx++)
    {
        (*results).r0t.row(time_idx) = p_se_components.row(time_idx);
        for (l = 0; l < p_se_components.cols(); l++)
        {
            if (transitionMode == "exponential")
            {
                (*results).r0t(time_idx, l) /= (
                        -std::log(1.0-(*results).p_ir(time_idx))/offset(time_idx));
                (*results).effR0(time_idx, l) = (*results).r0t(time_idx, l)
                                            *((*results).S(time_idx, l))/N(l);
            }
            else if (transitionMode == "path_specific")
            {
                (*results).r0t(time_idx, l) *= inf_mean;
                (*results).effR0(time_idx, l) = (*results).r0t(time_idx, l)
                                            *((*results).S(time_idx, l))/N(l);
               
            }
            else // Weibull
            {
                (*results).r0t(time_idx, l) *= 
                    IR_transition_dist -> getAvgMembership();
                (*results).effR0(time_idx, l) = (*results).r0t(time_idx, l)
                                            *((*results).S(time_idx, l))/N(l);
            }
        }
    }

    // Calculate next generation matrices 
    unsigned int size_idx;
    has_ts_spatial = (TDM_vec[0].size() > 0);
    Eigen::MatrixXd CombinedDM = Eigen::MatrixXd::Zero(nLoc, nLoc);
    Eigen::MatrixXd tmpDM = Eigen::MatrixXd::Zero(nLoc, nLoc);
    for (size_idx = 0; size_idx < DM_vec.size(); size_idx++) 
    {
        CombinedDM += (DM_vec[size_idx]*((*results).rho(size_idx)));
    }
    for (k = 0; k < nLoc; k++)
    {
        CombinedDM(k,k) = 1.0;
    }
    tmpDM = CombinedDM;

    for (time_idx = 0; time_idx < nTpt; time_idx++)
    {
        Eigen::MatrixXd G(nLoc, nLoc); 
        for (i = 0; i < nLoc; i++)
        {
            for (l = 0; l < nLoc; l++)
            {
                if (has_ts_spatial)
                {
                    tmpDM = CombinedDM;
                    for (size_idx = 0; 
                         size_idx < TDM_vec[time_idx].size(); size_idx++) 
                    {
                        tmpDM += (TDM_vec[time_idx][size_idx]*
                                ((*results).rho(size_idx + DM_vec.size())));
                    }
                }
                G(i,l) = ((*results).I(time_idx, l) == 0 ? 
                            0.0 :
                            (*results).S(time_idx, i)/
                            (1.0*(*results).I(time_idx, l)) *
                            (1.0 - std::exp(-offset(time_idx) * tmpDM(i,l)
                                            * ((*results).I(time_idx, l) 
                                            * (p_se_components(time_idx, l)
                                            /  N(l)))
                                            )));
            }
        }
        Eigen::VectorXd colsums = G.colwise().sum();
        (*results).rEA.row(time_idx) = colsums;
    }
          
    double _1mpIR_cum = 1.0;
    int infTime = 0;
    Eigen::MatrixXd finalEARVal = (*results).rEA.row(nTpt - 1);
    if (transitionMode == "exponential")
    {
        for (time_idx = 0; time_idx < nTpt; time_idx++)
        {
            _1mpIR_cum = (1-(*results).p_ir(time_idx, 0));
            for (l = time_idx+1; l < nTpt; l++)
            {
                (*results).rEA.row(time_idx) += 
                    _1mpIR_cum*(*results).rEA.row(l);
                _1mpIR_cum *= (1 - (*results).p_ir(l, 0)); 
            }
            while (_1mpIR_cum > 1e-12)
            {
                (*results).rEA.row(time_idx) += _1mpIR_cum*finalEARVal; 
                _1mpIR_cum*=(1-(*results).p_ir(nTpt - 1,0));
            }
        }
    }
    else
    {
        for (time_idx = 0; time_idx < nTpt; time_idx++)
        {
            _1mpIR_cum = 1.0;
            infTime = 0;
            for (i = 0; i < offset(time_idx); i++)
            {
                _1mpIR_cum *= (1 - I_to_R_prior(infTime, 5)); 
                infTime ++;
            }
            if (transitionMode == "path_specific")
            {
                for (l = time_idx+1; (l < nTpt && infTime < I_paths[0].rows()); l++)
                {
                    (*results).rEA.row(time_idx) += 
                        _1mpIR_cum*(*results).rEA.row(l);
                    for (i = 0; i < offset(time_idx); i++)
                    {
                        _1mpIR_cum *= (1 - I_to_R_prior(infTime, 5)); 
                        infTime ++;
                    }
                }
                while (_1mpIR_cum > 1e-12 && infTime < I_paths[0].rows())
                {
                    (*results).rEA.row(time_idx) += _1mpIR_cum*finalEARVal; 
                    _1mpIR_cum *= (1 - I_to_R_prior(infTime, 5)); 
                    infTime ++;
                }
            }
            else
            {
                for (l = time_idx+1; (l < nTpt && infTime < I_paths[0].rows()); l++)
                {
                    (*results).rEA.row(time_idx) += 
                        _1mpIR_cum*(*results).rEA.row(l);
                    for (i = 0; i < offset(time_idx); i++)
                    {
                        _1mpIR_cum *= (1 - 
                            IR_transition_dist -> getTransitionProb(
                                infTime, infTime+1)); 
                        infTime ++;
                    }
                }
                while (_1mpIR_cum > 1e-12 && infTime < I_paths[0].rows())
                {
                    (*results).rEA.row(time_idx) += _1mpIR_cum*finalEARVal; 
                    _1mpIR_cum *= (1 - IR_transition_dist -> getTransitionProb(
                                infTime, infTime+1));
                    infTime ++;
                }
            }
        }
    }
}

SEIR_sim_node::~SEIR_sim_node()
{
    delete generator;
}

