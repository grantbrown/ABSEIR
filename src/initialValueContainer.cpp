#include <Rcpp.h>
#include <initialValueContainer.hpp>


using namespace Rcpp;

initialValueContainer::initialValueContainer()
{
    // Do nothing
}

int initialValueContainer::getModelComponentType()
{
    return(LSS_INIT_CONTAINER_TYPE);
}

void initialValueContainer::summary()
{
    Rcpp::Rcout << "Initial Compartment Values - Summary\n" << 
                   "------------------------------------\n";
    int i;
    int sz = S0.size();
    Rcpp::Rcout << "    S0, E0, I0, R0 \n";
    for (i = 0; i < sz; i++)
    {
        Rcpp::Rcout << "    " <<
                       S0(i) << ", " <<
                       E0(i) << ", " <<
                       I0(i) << ", " <<
                       R0(i) << "\n";
    }
    Rcpp::Rcout << "\n";
}

void initialValueContainer::setInitialValues(SEXP S0_, SEXP E0_, SEXP I0_, SEXP R0_)
{

    Rcpp::IntegerVector S0_vec(S0_);
    Rcpp::IntegerVector E0_vec(E0_);
    Rcpp::IntegerVector I0_vec(I0_);
    Rcpp::IntegerVector R0_vec(R0_);
    if (S0_vec.length() != E0_vec.length() || 
        E0_vec.length() != I0_vec.length() ||
        I0_vec.length() != R0_vec.length())
    {
        Rcpp::stop("Init compartment lengths do not match\n");
    }
    int i;
    S0 = Eigen::VectorXi(E0_vec.length());
    E0 = Eigen::VectorXi(E0_vec.length());
    I0 = Eigen::VectorXi(E0_vec.length());
    R0 = Eigen::VectorXi(E0_vec.length());

    for (i = 0; i < S0_vec.length(); i++){
        S0(i) = S0_vec(i);
        E0(i) = E0_vec(i);
        I0(i) = I0_vec(i);
        R0(i) = R0_vec(i);
    }
}


initialValueContainer::~initialValueContainer()
{
    if (prot !=0 ){
        Rcpp::stop("can't delete initialValueContainer, still being used.\n");
    }
}

RCPP_MODULE(mod_initialValueContainer)
{
    using namespace Rcpp;
    class_<initialValueContainer>( "initialValueContainer" )
    .constructor()
    .method("setInitialValues", &initialValueContainer::setInitialValues);
}


