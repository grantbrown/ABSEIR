#ifndef ACTOR_SEIRSIM_HEADER
#define ACTOR_SEIRSIM_HEADER
#include <map>
#include <vector>
#include <random>
#include <sstream>
#include <iostream>
#include <Eigen/Core>
#include <ABSEIR_constants.hpp>
#include <dataModel.hpp>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <memory>

using namespace std;

/*
using sim_atom = "sim";
using sim_result_atom = "sim_rslt";
using sample_atom = "sample";
using wakeup_atom = "wakeup";
using exit_atom = "exit";
*/

struct simulationResultSet;
class transitionDistribution;
class NodePool;

struct instruction{
   int param_idx; 
   std::string action_type;
   Eigen::VectorXd params;
};

class SEIR_sim_node {
    public:
        SEIR_sim_node(int random_seed,
                      Eigen::VectorXi S0,
                      Eigen::VectorXi E0,
                      Eigen::VectorXi I0,
                      Eigen::VectorXi R0,
                      Eigen::VectorXd offset,
                      Eigen::MatrixXi Y,
                      MatrixXb na_mask,
                      std::vector<Eigen::MatrixXd> DM_vec,
                      std::vector<std::vector<Eigen::MatrixXd> > tdm_vec,
                      Eigen::MatrixXd X, 
                      Eigen::MatrixXd X_rs,
                      std::string transitionMode,
                      Eigen::MatrixXd ei_prior,
                      Eigen::MatrixXd ir_prior,
                      double avgI,
                      Eigen::VectorXd sp_prior,
                      Eigen::VectorXd se_prec,
                      Eigen::VectorXd rs_prec,
                      Eigen::VectorXd se_mean,
                      Eigen::VectorXd rs_mean,
                      double phi,
                      int data_compartment,
                      bool cumulative);
        ~SEIR_sim_node();
        simulationResultSet simulate(Eigen::VectorXd param_vals, bool keepCompartments);

    private: 
        void calculateReproductiveNumbers(simulationResultSet* input);
        int sim_width;
        unsigned int random_seed;
        Eigen::VectorXi S0;
        Eigen::VectorXi E0;
        Eigen::VectorXi I0;
        Eigen::VectorXi R0;
        Eigen::VectorXd offset;
        Eigen::MatrixXi Y;
        Eigen::MatrixXi E_paths;
        Eigen::MatrixXi I_paths;
        MatrixXb na_mask;
        std::vector<Eigen::MatrixXd> DM_vec;
        std::vector<std::vector<Eigen::MatrixXd> > TDM_vec;
        Eigen::MatrixXd X;
        Eigen::MatrixXd X_rs;
        std::string transitionMode;
        Eigen::MatrixXd E_to_I_prior;
        Eigen::MatrixXd I_to_R_prior;
        double inf_mean;
        Eigen::VectorXd spatial_prior;
        Eigen::VectorXd exposure_precision;
        Eigen::VectorXd reinfection_precision;
        Eigen::VectorXd exposure_mean;
        Eigen::VectorXd reinfection_mean;
        /** General E to I transition Distribution*/
        std::unique_ptr<transitionDistribution> EI_transition_dist;
        /** General I to R transition Distribution*/
        std::unique_ptr<transitionDistribution> IR_transition_dist;
        double phi;
        int seed;
        double value;
        bool has_spatial;
        bool has_reinfection;
        int total_size;
        int data_compartment;
        bool cumulative;
        mt19937* generator;
        std::normal_distribution<double> overdispersion_distribution;
};



class NodeWorker{
    public:
        NodeWorker(NodePool* pl, 
                   int random_seed,
                   Eigen::VectorXi S0,
                   Eigen::VectorXi E0,
                   Eigen::VectorXi I0,
                   Eigen::VectorXi R0,
                   Eigen::VectorXd offset,
                   Eigen::MatrixXi Y,
                   MatrixXb na_mask,
                   std::vector<Eigen::MatrixXd> DM_vec,
                   std::vector<std::vector<Eigen::MatrixXd> > TDM_vec,
                   Eigen::MatrixXd X, 
                   Eigen::MatrixXd X_rs,
                   std::string transitionMode,
                   Eigen::MatrixXd ei_prior,
                   Eigen::MatrixXd ir_prior,
                   double avgI,
                   Eigen::VectorXd sp_prior,
                   Eigen::VectorXd se_prec,
                   Eigen::VectorXd rs_prec,
                   Eigen::VectorXd se_mean,
                   Eigen::VectorXd rs_mean,
                   double phi,
                   int data_compartment,
                   bool cumulative);
        void operator()();
    private:
        NodePool* pool;
        std::unique_ptr<SEIR_sim_node> node;
};

class NodePool{
    public:
        NodePool(std::vector<double>* result_pointer,
                 std::vector<simulationResultSet>* result_complete_pointer,
                 std::vector<int>* index_pointer,
                 int,
                 int random_seed,
                 Eigen::VectorXi S0,
                 Eigen::VectorXi E0,
                 Eigen::VectorXi I0,
                 Eigen::VectorXi R0,
                 Eigen::VectorXd offset,
                 Eigen::MatrixXi Y,
                 MatrixXb na_mask,
                 std::vector<Eigen::MatrixXd> DM_vec,
                 std::vector<std::vector<Eigen::MatrixXd> > TDM_vec,
                 Eigen::MatrixXd X, 
                 Eigen::MatrixXd X_rs,
                 std::string transitionMode,
                 Eigen::MatrixXd ei_prior,
                 Eigen::MatrixXd ir_prior,
                 double avgI,
                 Eigen::VectorXd sp_prior,
                 Eigen::VectorXd se_prec,
                 Eigen::VectorXd rs_prec,
                 Eigen::VectorXd se_mean,
                 Eigen::VectorXd rs_mean,
                 double phi,
                 int data_compartment,
                 bool cumulative
              );
        void awaitFinished();
        void enqueue(std::string action_type, int param_idx, Eigen::VectorXd params);
        std::vector<double>* result_pointer;
        std::vector<simulationResultSet>* result_complete_pointer;
        std::vector<int>* index_pointer;
        ~NodePool();

    private:
        friend class NodeWorker;
        std::vector<std::thread> nodes;
        std::deque<instruction>  tasks;
        std::atomic_int nBusy;

        std::mutex queue_mutex;
        std::mutex result_mutex;
        std::condition_variable condition;
        std::condition_variable finished;
        bool exit; 
};




#endif
