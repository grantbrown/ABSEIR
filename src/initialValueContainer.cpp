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
        ::Rf_error("Init compartment lengths do not match\n");
    }
    int i;
    for (i = 0; i < S0_vec.length(); i++){
        S0.push_back(S0_vec[i]);
        E0.push_back(E0_vec[i]);
        I0.push_back(I0_vec[i]);
        R0.push_back(R0_vec[i]);
    }
}


initialValueContainer::~initialValueContainer()
{
    if (prot !=0 ){
        ::Rf_error("can't delete initialValueContainer, still being used.\n");
    }
}

RCPP_MODULE(mod_initialValueContainer)
{
    using namespace Rcpp;
    class_<initialValueContainer>( "initialValueContainer" )
    .constructor()
    .method("setInitialValues", &initialValueContainer::setInitialValues);
}


