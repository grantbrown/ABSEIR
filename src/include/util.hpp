#ifndef SPATIALSEIR_UTILDEF
#define SPATIALSEIR_UTILDEF
#include <Rcpp.h>
#include <Eigen/Core>

using namespace Rcpp;

class compartment_tap{
    public:
        compartment_tap(int nrow, int ncol);
        virtual void push(Eigen::VectorXi current_comp);
        Eigen::VectorXi get(int lag);

    private:
        int idx;
        int nLags;
        std::vector<int> beenSet;
        Eigen::MatrixXi compartment;
};


#endif
