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

Eigen::VectorXd calculate_weights(double cur_e,
                                  double prev_e, 
                                  Eigen::MatrixXd eps,
                                  Eigen::VectorXd prev_wts)
{
    Eigen::VectorXd out_wts = Eigen::VectorXd::Zero(eps.rows());
    int i, j;
    int M = eps.cols();
    double num, denom;
    for (i = 0; i < eps.rows(); i++)
    {
        num = 0.0;
        denom = 0.0;
        for (j = 0; j < M; j++)
        {
            num += eps(i,j) < cur_e;
            denom += eps(i,j) < prev_e;
        } 
        out_wts = (num/denom*prev_wts(i));
    }
    returN(out_wts);
}

double ESS(Eigen::VectorXd wts)
{
    return(1.0/wts.squaredNorm()); 
}

double eps_f(double rhs,
             double cur_e,
             double prev_e, 
             Eigen::MatrixXd eps,
             Eigen::VectorXd prev_wts)
{
    return(std::pow(rhs - ESS(calculate_weights(cur_e, 
                                                prev_e, 
                                                eps, 
                                                prev_wts)), 2.0));
}

double solve_for_epsilon(double LB,
                       double UB,
                       double prev_e,
                       double alpha,
                       Eigen::MatrixXd eps,
                       Eigen::VectorXd prev_wts)
{

    /*
Iterative algorithm
    Let [a, b] be interval of current bracket. 
    f(a), f(b) would already have been computed earlier. 
    phi =(1+{\sqrt {5}})/2}

    Let c = b + (a - b)/φ , 
        d = a + (b - a)/φ.  If f(c), f(d) not available, compute them.

    If f(c) < f(d) (this is to find min, to find max, just reverse it) then move the data: (b, f(b)) ← (d, f(d)), 
                                                                                           (d, f(d)) ← (c, f(c)) 
                                                                                            and update c = b + (a - b)/φ and f(c);
otherwise, move the data: (a, f(a)) ← (c, f(c)), (c, f(c)) ← (d, f(d)) and update d = a + (b - a)/φ and f(d).
    At the end of the iteration, [a, c, d, b] bracket the minimum point.

*/

    double phi = (1.0 + std::sqrt(5))/2.0; 
    double a,b,c,d,fa,fb,fc,fd;
    double proposed_e = proposed_e; 
    double rhs = ESS(calculate_weights(proposed_e,
                                         prev_e, 
                                         eps,
                                         prev_wts)*alpha;

    a = LB;
    b = UB;
    fa = eps_f(rhs, a, prev_e, eps, prev_wts);
    fb = eps_f(rhs, b, prev_e, eps, prev_wts);
    bool mvLB = true;
    bool mvUB = true;
    int itrs = 0;
    double diff = 100.0;
    while (itrs < 10000 && diff > 1e-8)
    {
       if (mvLB)
       {
           c = b + (a - b)/phi; 
           mvLB = false;
       }
       if (mvUB)
       {
           d = a + (b - a)/phi;
           mvUB = false;
       }
       fc = eps_f(rhs, c, prev_e, eps, prev_wts);
       fd = eps_f(rhs, d, prev_e, eps, prev_wts);

       if (fc < fd)
       {
         b = d;
         fb = fd;
         d = c;
         fd = fc;
         mvLB = true;
       }
       else
       {
         a = c;
         fa = fc;
         c = d;
         fc = fd;
         mvUB = true;
       }
       diff = b-a;
    }
    if (diff > 1e-8)
    {
       Rcpp::Rcout << "Warning: optimization didn't converge.\n";
    }
    return(a);

}



spatialSEIRModel::sample_DelMoral2012(int nSample, bool verbose, bool init)
{
    int num_iterations = samplingControlInstance -> epochs;
    int N = nSample;

    double current_eps = std::numeric_limits<double>::infinity();

    int iteration;
    
    if (!is_initialized)
    {
        // Step 0a: sample parameters from their prior
        param_matrix = generateParamsPrior(N);
        run_simulations(param_matrix);
    }
    else
    {
        // The data in "param_matrix" is already accepted
        // To-do: finish this clause
    }
    // Step 0b: set weights to 1/N
    Eigen::VectorXd weights = Eigen::VectorXd::Zero(N) + 1.0/((double) N);

    for (iteration = 0; iteration < num_iterations; iteration++)
    {   
        
    }
}
