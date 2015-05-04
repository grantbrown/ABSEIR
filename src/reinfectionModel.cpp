#include <Rcpp.h>
#include <reinfectionModel.hpp>


using namespace Rcpp;


reinfectionModel::reinfectionModel(SEXP reinfectMode)
{
     Rcpp::IntegerVector modeVec(reinfectMode);
     reinfectionMode = modeVec[0];
     betaPriorPrecision.push_back(-1); 
}

int reinfectionModel::getModelComponentType()
{
    return(LSS_REINFECTION_MODEL_TYPE);
}

void reinfectionModel::buildReinfectionModel(SEXP _X, SEXP _priorMean, SEXP _prec)
{
    Rcpp::NumericMatrix inX(_X);
    Rcpp::NumericVector inPrecision(_prec);
    Rcpp::NumericVector priorMeans(_priorMean);


    if (priorMeans.length() != inPrecision.length() || inPrecision.length() != inX.ncol())
    {
        ::Rf_error("Number of parameters, prior means, or precisions does not equal the number of supplied covariates.\n");
    }
    int i,j;
    for (i = 0; i < inX.ncol(); i++)
    {
        betaPriorPrecision.push_back(inPrecision[i]);
        betaPriorMean.push_back(priorMeans[i]);
    }
    X_rs = Eigen::MatrixXd(inX.nrow(), inX.ncol());
    for (i = 0; i < inX.nrow(); i++)
    {
        for (j = 0; j < inX.ncol(); j++)
        {
           X_rs(i,j) = inX[i,j]; 
        }
    }
}

reinfectionModel::~reinfectionModel()
{
    if (prot !=0 ){
        ::Rf_error("can't delete reinfectionModel, still being used.\n");
    }
}

RCPP_MODULE(mod_reinfectionModel)
{
    using namespace Rcpp;
    class_<reinfectionModel>( "reinfectionModel" )
    .constructor<SEXP>()
    .method("buildReinfectionModel", &reinfectionModel::buildReinfectionModel);
}


