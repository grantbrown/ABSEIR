#include <Rcpp.h>
#include <exposureModel.hpp>

using namespace Rcpp;

exposureModel::exposureModel(SEXP _X, SEXP _ntpt, SEXP _nloc, SEXP _priorMean, SEXP _prec)
{
    int i,j;
    Rcpp::NumericMatrix inX(_X);
    Rcpp::NumericVector inTpt(_ntpt);
    Rcpp::NumericVector inLoc(_nloc);
    Rcpp::NumericVector inPrecision(_prec);
    Rcpp::NumericVector priorMeans(_priorMean);

    nTpt = inTpt[0];
    nLoc = inLoc[0];

    X = Eigen::MatrixXd(inX.nrow(), inX.ncol());

    if (inX.nrow() != ((nTpt) *(nLoc)))
    {
        Rcpp::Rcout << "Error: Covariate matrix has invalid dimensions." 
                    << " X should be (n*t) by p\n";
        ::Rf_error("error building exposure model");
    }
    if (inPrecision.length() != priorMeans.length() || priorMeans.length() != X.cols())
    {
        Rcpp::Rcout << "Initial parameters have different length then number of columns of X\n";
        Rcpp::Rcout << "ncol(X): " << X.cols() << "\n";
        Rcpp::Rcout << "length(betaPriorMean): " << priorMeans.length() << "\n";
        Rcpp::Rcout << "length(betaPriorPrecision): " << inPrecision.length() << "\n";
        ::Rf_error("error building exposure model");
    }
    for (i = 0; i < nTpt; i++)
    {
        offset.push_back(1.0);
    }
    for (i = 0; i < X.rows(); i++)
    {
        for (j = 0; j < X.cols(); j++)
        {
            X(i,j) = inX[i,j];
        }
    }
    for (i = 0; i < X.cols(); i++)
    {
        betaPriorMean.push_back(priorMeans[i]);
        betaPriorPrecision.push_back(inPrecision[i]);
    }
}

int exposureModel::getModelComponentType()
{
    return(LSS_EXPOSURE_MODEL_TYPE);
}


void exposureModel::setOffset(NumericVector offsets)
{
    if (offsets.length() != (nTpt))
    {
        Rcpp::Rcout << "Error: offsets must have length equal to the number of time points.\n";
        ::Rf_error("Invalid offsets.");
    }
    int i;
    for (i = 0; i < offsets.length(); i++)
    {
        offset.push_back(offsets[i]);
    }
}

Rcpp::NumericVector exposureModel::getOffset()
{
    return(Rcpp::NumericVector(offset.begin(), offset.end()));
}

exposureModel::~exposureModel()
{
    if (prot !=0 ){
        ::Rf_error("can't delete exposureModel, still being used.\n");
    }
}

RCPP_MODULE(mod_exposureModel)
{
    using namespace Rcpp;
    class_<exposureModel>( "exposureModel" )
    .constructor<SEXP,SEXP,SEXP,SEXP,SEXP>()
    .property("offsets", &exposureModel::getOffset, &exposureModel::setOffset);
}


