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
        E_to_I_params = Eigen::MatrixXd(1,5);
        I_to_R_params = Eigen::MatrixXd(1,5);
        for (int i = 0; i < 5; i++)
        {
            E_to_I_params(0,i) = 1.0;
            I_to_R_params(0,i) = 1.0;
        }
    }
}

transitionPriors::transitionPriors(transitionPriors* tocopy)
{
    Eigen::MatrixXd E_to_I_paramsc = (tocopy -> E_to_I_params);
    E_to_I_params = E_to_I_paramsc;

    Eigen::MatrixXd I_to_R_paramsc = (tocopy -> I_to_R_params);
    I_to_R_params = I_to_R_paramsc;

    mode = (tocopy -> mode);
    max_latent = (tocopy -> max_latent);
    max_infectious = (tocopy -> max_infectious);
    inf_mean = (tocopy -> inf_mean);
}

int transitionPriors::getModelComponentType()
{
    return(LSS_TRANSITION_MODEL_TYPE);
}

void transitionPriors::setUniformExpPriors()
{
    E_to_I_params = Eigen::MatrixXd(2, 1);
    I_to_R_params = Eigen::MatrixXd(2, 1);

    E_to_I_params(0,0) = 1.0; 
    E_to_I_params(1,0) = 1.0;
    I_to_R_params(0,0) = 1.0; 
    I_to_R_params(1,0) = 1.0;
}

void transitionPriors::setPriorsForWeibull(SEXP _E_to_I_params,
                                           SEXP _I_to_R_params,
                                           SEXP _maxEI,
                                           SEXP _maxIR)
{
    Rcpp::NumericVector EIP(_E_to_I_params);
    Rcpp::NumericVector IRP(_I_to_R_params);
    Rcpp::NumericVector maxEI(_maxEI);
    Rcpp::NumericVector maxIR(_maxIR);
    E_to_I_params = Eigen::MatrixXd(5, 1);   
    I_to_R_params = Eigen::MatrixXd(5, 1);   
    E_to_I_params(0,0) = EIP(0);
    E_to_I_params(1,0) = EIP(1);
    E_to_I_params(2,0) = EIP(2);
    E_to_I_params(3,0) = EIP(3);
    E_to_I_params(4,0) = maxEI(0);

    I_to_R_params(0,0) = IRP(0);
    I_to_R_params(1,0) = IRP(1);
    I_to_R_params(2,0) = IRP(2);
    I_to_R_params(3,0) = IRP(3);
    I_to_R_params(4,0) = maxIR(0);
    inf_mean = IRP(1)*std::tgamma(1+1.0/IRP(0));
}


void transitionPriors::setPathSpecificPriors(SEXP _Zmat1, SEXP _Zmat2, 
        SEXP _avgI)
{
    Rcpp::NumericMatrix Zmat1(_Zmat1); 
    Rcpp::NumericMatrix Zmat2(_Zmat1); 
    Rcpp::NumericVector avgI(_avgI);
    inf_mean = avgI(0);

    E_to_I_params = Eigen::MatrixXd(Zmat1.nrow(), Zmat1.ncol());   
    I_to_R_params = Eigen::MatrixXd(Zmat2.nrow(), Zmat2.ncol());   
    int i,j;
    for (i = 0; i < Zmat1.ncol(); i++)
    {
        for (j = 0; j < Zmat1.nrow(); j++)
        {
            E_to_I_params(j,i) = Zmat1(j,i);
        }
    }
    for (i = 0; i < Zmat2.ncol(); i++)
    {
        for (j = 0; j < Zmat2.nrow(); j++)
        {
            I_to_R_params(j,i) = Zmat2(j,i);
        }
    }
    max_latent = E_to_I_params.rows();
    max_infectious = I_to_R_params.rows();
}

void transitionPriors::setPriorsFromProbabilities(SEXP p_ei, SEXP p_ir, 
                                                  SEXP p_ei_ess, SEXP p_ir_ess)
{
    E_to_I_params = Eigen::MatrixXd(2, 1);
    I_to_R_params = Eigen::MatrixXd(2, 1);

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

    E_to_I_params(0,0) = pEIess;
    E_to_I_params(1,0) = pEIess/(gamma_ei);

    I_to_R_params(0,0) = pIRess;
    I_to_R_params(1,0) = pIRess/(gamma_ir); 
}

void transitionPriors::summary()
{
    Rcpp::Rcout << "Transition Priors Summary\n";
    Rcpp::Rcout << "-------------------------\n";
    Rcpp::Rcout << "    transition mode: " << mode << "\n";
    Rcpp::Rcout << "    E to I transition prior params: (";
    for (int i = 0; i < E_to_I_params.rows() - 1; i++)
    {
        Rcpp::Rcout << E_to_I_params(i,0) << ", ";
    }
    Rcpp::Rcout << E_to_I_params(E_to_I_params.rows() - 1,0) << ")\n ";

    Rcpp::Rcout << "    I to R transition prior params: (";
    for (int i = 0; i < I_to_R_params.rows() - 1; i++)
    {
        Rcpp::Rcout << I_to_R_params(i,0) << ", ";
    }
    Rcpp::Rcout << I_to_R_params(I_to_R_params.rows() - 1,0) << ")\n ";
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
    .method("setPriorsForWeibull", &transitionPriors::setPriorsForWeibull)
    .method("summary", &transitionPriors::summary);
}


