#include <Rcpp.h>
#include <string>
#include <Rmath.h>
#include <transitionPriors.hpp>


using namespace Rcpp;

transitionPriors::transitionPriors(SEXP _mode)
{
    Rcpp::CharacterVector inMode(_mode);
    mode = inMode(0);
    // Store dummy transition info
    inf_mean = -1.0;
    if (mode == "exponential")
    {
        setUniformExpPriors();
    }
    else
    {
        gamma_ei_params = Eigen::MatrixXd(1,5);
        gamma_ir_params = Eigen::MatrixXd(1,5);
        for (int i = 0; i < 5; i++)
        {
            gamma_ei_params(0,i) = 1.0;
            gamma_ir_params(0,i) = 1.0;
        }
    }
}

int transitionPriors::getModelComponentType()
{
    return(LSS_TRANSITION_MODEL_TYPE);
}

void transitionPriors::setUniformExpPriors()
{
    gamma_ei_params = Eigen::MatrixXd(2, 1);
    gamma_ir_params = Eigen::MatrixXd(2, 1);

    gamma_ei_params(0,0) = 1.0; 
    gamma_ei_params(1,0) = 1.0;
    gamma_ir_params(0,0) = 1.0; 
    gamma_ir_params(1,0) = 1.0;
}

void transitionPriors::setPathSpecificPriors(SEXP _Zmat1, SEXP _Zmat2, 
        SEXP _avgI)
{
    Rcpp::NumericMatrix Zmat1(_Zmat1); 
    Rcpp::NumericMatrix Zmat2(_Zmat1); 
    Rcpp::NumericVector avgI(_avgI);
    inf_mean = avgI(0);

    gamma_ei_params = Eigen::MatrixXd(Zmat1.nrow(), Zmat1.ncol());   
    gamma_ir_params = Eigen::MatrixXd(Zmat2.nrow(), Zmat2.ncol());   
    int i,j;
    for (i = 0; i < Zmat1.ncol(); i++)
    {
        for (j = 0; j < Zmat1.nrow(); j++)
        {
            gamma_ei_params(j,i) = Zmat1(j,i);
        }
    }
    for (i = 0; i < Zmat2.ncol(); i++)
    {
        for (j = 0; j < Zmat2.nrow(); j++)
        {
            gamma_ir_params(j,i) = Zmat2(j,i);
        }
    }
    max_latent = gamma_ei_params.rows();
    max_infectious = gamma_ir_params.rows();
}

void transitionPriors::setPriorsFromProbabilities(SEXP p_ei, SEXP p_ir, 
                                                  SEXP p_ei_ess, SEXP p_ir_ess)
{
    gamma_ei_params = Eigen::MatrixXd(2, 1);
    gamma_ir_params = Eigen::MatrixXd(2, 1);

    double pEI, pIR;
    int pEIess, pIRess;
    double gamma_ei, gamma_ir;
    Rcpp::NumericVector p_ei_vec(p_ei);
    Rcpp::NumericVector p_ir_vec(p_ir);
    Rcpp::IntegerVector p_ei_ess_vec(p_ei_ess);
    Rcpp::IntegerVector p_ir_ess_vec(p_ir_ess);

    pEI = p_ei_vec(0); pIR = p_ir_vec(0);
    pEIess = p_ei_ess_vec(0); pIRess = p_ir_ess_vec(0);


    gamma_ei = -std::log(1-pEI);
    gamma_ir = -std::log(1-pIR);

    gamma_ei_params(0,0) = pEIess;
    gamma_ei_params(1,0) = pEIess/(gamma_ei);

    gamma_ir_params(0,0) = pIRess;
    gamma_ir_params(1,0) = pIRess/(gamma_ir); 
}

void transitionPriors::summary()
{
    if (gamma_ei_params.size() == 2){
        Rcpp::Rcout << "gamma_ei parameters: " << gamma_ei_params(0) << ", " << 1/gamma_ei_params(1) << "\n";
        Rcpp::Rcout << "gamma_ir parameters: " << gamma_ir_params(0) << ", " << 1/gamma_ir_params(1) << "\n";
    }
}

transitionPriors::~transitionPriors()
{
    if (prot !=0 ){
        Rcpp::stop("can't delete transitionPriors, still being used.\n");
    }
}

RCPP_MODULE(mod_transitionPriors)
{
    using namespace Rcpp;
    class_<transitionPriors>( "transitionPriors" )
    .constructor<SEXP>()
    .method("setUniformExpPriors", &transitionPriors::setUniformExpPriors)
    .method("setPathSpecificPriors", &transitionPriors::setPathSpecificPriors)
    .method("setPriorsFromProbabilities", &transitionPriors::setPriorsFromProbabilities)
    .method("summary", &transitionPriors::summary);
}


