#ifndef ACTOR_SEIRSIM_HEADER
#define ACTOR_SEIRSIM_HEADER
#include <map>
#include <vector>
#include <random>
#include <sstream>
#include <Rcpp.h>
#include <iostream>
#include <Eigen/Core>
#include "caf/all.hpp"

using namespace Rcpp;
using namespace std;
using namespace caf;

using sim_atom = atom_constant<atom("sim")>;
using wakeup_atom = atom_constant<atom("wakeup")>;
using exit_atom = atom_constant<atom("exit")>;

class SEIR_sim_node : public event_based_actor {
    public:
        SEIR_sim_node(int sim_width,
                      int random_seed,
                      std::vector<int> S0,
                      std::vector<int> E0,
                      std::vector<int> I0,
                      std::vector<int> R0,
                      Eigen::MatrixXi I_star,
                      std::vector<Eigen::MatrixXd> DM_vec,
                      Eigen::MatrixXd X, 
                      Eigen::MatrixXd X_rs,
                      std::vector<double> ei_prior,
                      std::vector<double> ir_prior,
                      std::vector<double> se_prec,
                      std::vector<double> rs_prec,
                      std::vector<double> se_mean,
                      std::vector<double> rs_mean,
                      double phi,
                      actor parent);
        ~SEIR_sim_node();
    protected:
        double simulate(Eigen::VectorXd param_vals);
        behavior make_behavior() override;
    private: 
        unsigned int random_seed;
        std::vector<int> S0;
        std::vector<int> E0;
        std::vector<int> I0;
        std::vector<int> R0;
        Eigen::MatrixXi I_star;
        std::vector<Eigen::MatrixXd> DM_vec;
        Eigen::MatrixXd X;
        Eigen::MatrixXd X_rs;
        std::vector<double> E_to_I_prior;
        std::vector<double> I_to_R_prior;
        std::vector<double> exposure_precision;
        std::vector<double> reinfection_precision;
        std::vector<double> exposure_mean;
        std::vector<double> reinfection_mean;
        double phi;
        int sim_width;
        int seed;
        double value;
        actor parent;
        scoped_actor* self;
        mt19937* generator;
        behavior    alive;
};




#endif
