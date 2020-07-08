#ifndef ACTOR_SEIRSIM_HEADER
#define ACTOR_SEIRSIM_HEADER

#include <map>
#include <vector>
#include <random>
#include <sstream>
#include <iostream>
#include <Eigen/Core>
#include <ABSEIR_constants.hpp>
#include "./dataModel.hpp"
#include "./distanceModel.hpp"
#include "./exposureModel.hpp"
#include "./initialValueContainer.hpp"
#include "./reinfectionModel.hpp"
#include "./samplingControl.hpp"
#include "./SEIRSimNodes.hpp"
#include "./transitionPriors.hpp"
#include "./transitionDistribution.hpp"
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
class NodeWorker;

struct instruction{
   int param_idx; 
   std::string action_type;
   Eigen::VectorXd params;
};

class SEIR_sim_node {
    public:
        SEIR_sim_node(NodeWorker* worker,
                      size_t seed_offs,
                      dataModel* dataModelPtr,
                      exposureModel* exposureModelPtr, 
                      reinfectionModel* reinfectionModelPtr,
                      distanceModel* distanceModelPtr,
                      transitionPriors* transitionPriorsPtr, 
                      initialValueContainer* initialValueContainerPtr,
                      samplingControl* samplingControlPtr
                );
        ~SEIR_sim_node();
        std::deque<std::string> messages;
        simulationResultSet simulate(Eigen::VectorXd param_vals, bool keepCompartments);

    private: 
        NodeWorker* parent;
        unsigned int random_seed;
        /* GB Notes, 2020-07-08
         *
         * Since we're moving to having a thread specific copy of model 
         * components, we should get rid of all the duplicate stuff stored 
         * in this class. Leaving it for now during the transition. 
         *
         */

        
        dataModel* dataModelInstance;
        exposureModel* exposureModelInstance; 
        reinfectionModel* reinfectionModelInstance;
        distanceModel* distanceModelInstance;
        transitionPriors* transitionPriorsInstance; 
        initialValueContainer* initialValueContainerInstance;
        samplingControl* samplingControlInstance;

        Eigen::VectorXi S0;
        Eigen::VectorXi E0;
        Eigen::VectorXi I0;
        Eigen::VectorXi R0;
        Eigen::VectorXd offset;
        Eigen::MatrixXi Y;
        Eigen::VectorXd loc_weights;
        MatrixXb na_mask;
        int dataModelType;
        std::vector<Eigen::MatrixXd> DM_vec;
        std::vector<std::vector<Eigen::MatrixXd> > TDM_vec;
        std::vector<int> TDM_empty;

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
        double phi;
        int data_compartment;
        bool cumulative;
        int m;
		double lpow;

        std::vector<Eigen::MatrixXi> E_paths;
        std::vector<Eigen::MatrixXi> I_paths;

        /** General E to I transition Distribution*/
        std::unique_ptr<transitionDistribution> EI_transition_dist;
        /** General I to R transition Distribution*/
        std::unique_ptr<transitionDistribution> IR_transition_dist;
        void nodeMessage(std::string);

        int seed;
        double value;
        bool has_spatial;
        bool has_ts_spatial;
        bool has_reinfection;
        bool has_report_fraction;
        int total_size;

        int sim_width;
        mt19937* generator;
        std::normal_distribution<double> overdispersion_distribution;
};



class NodeWorker{
    public:
        NodeWorker(NodePool* pl, 
                   size_t seed_offs,
                   dataModel* dataModel,
                   exposureModel* exposureModel, 
                   reinfectionModel* reinfectionModel,
                   distanceModel* distanceModel,
                   transitionPriors* transitionPriors, 
                   initialValueContainer* initialValueContainer,
                   samplingControl* samplingControl
                   );

        void operator()();
        void addMessage(std::string);

    private:
        friend class SEIR_sim_node;
        NodePool* pool;
        std::unique_ptr<SEIR_sim_node> node;
};

class NodePool{
    public:
        NodePool(Eigen::MatrixXd* result_pointer,
                 std::vector<simulationResultSet>* result_complete_pointer,
                 std::vector<int>* index_pointer,
                 dataModel* dataModel,
                 exposureModel* exposureModel, 
                 reinfectionModel* reinfectionModel,
                 distanceModel* distanceModel,
                 transitionPriors* transitionPriors, 
                 initialValueContainer* initialValueContainer,
                 samplingControl* samplingControl
                 );

        void setResultsDest(Eigen::MatrixXd* result_pointer,
                            std::vector<simulationResultSet>* result_complete_pointer);
        void awaitFinished();
        void resolveMessages();
        void enqueue(std::string action_type, int param_idx, Eigen::VectorXd params);
        Eigen::MatrixXd* result_pointer;
        std::deque<std::string> messages;
        std::vector<simulationResultSet>* result_complete_pointer;
        std::vector<int>* index_pointer;
        ~NodePool();

    private:
        friend class NodeWorker;
        

#ifdef SPATIALSEIR_SINGLETHREAD
        std::vector<NodeWorker> nodes;
#else
        std::vector<std::thread> nodes;
#endif
        std::deque<instruction>  tasks;
        std::atomic_int nBusy;

        std::mutex queue_mutex;
        std::mutex result_mutex;
        std::condition_variable condition;
        std::condition_variable finished;
        bool exit; 
};




#endif
