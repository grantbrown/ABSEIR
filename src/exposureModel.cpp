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

    nTpt = inTpt(0);
    nLoc = inLoc(0);

    X = Eigen::MatrixXd(inX.nrow(), inX.ncol());
    offset = Eigen::VectorXd(nTpt);
    betaPriorMean = Eigen::VectorXd(priorMeans.size());
    betaPriorPrecision = Eigen::VectorXd(inPrecision.size());


    if (inX.nrow() != ((nTpt) *(nLoc)))
    {
        Rcpp::Rcout << "Error: Covariate matrix has invalid dimensions." 
                    << " X should be (n*t) by p\n";
        Rcpp::stop("error building exposure model");
    }
    if (inPrecision.length() != priorMeans.length() || priorMeans.length() != X.cols())
    {
        Rcpp::Rcout << "Initial parameters have different length then number of columns of X\n";
        Rcpp::Rcout << "ncol(X): " << X.cols() << "\n";
        Rcpp::Rcout << "length(betaPriorMean): " << priorMeans.length() << "\n";
        Rcpp::Rcout << "length(betaPriorPrecision): " << inPrecision.length() << "\n";
        Rcpp::stop("error building exposure model");
    }
    for (i = 0; i < nTpt; i++)
    {
        offset(i) = 1.0;
    }
    for (i = 0; i < X.rows(); i++)
    {
        for (j = 0; j < X.cols(); j++)
        {
            X(i,j) = inX(i,j);
        }
    }
    for (i = 0; i < X.cols(); i++)
    {
        betaPriorMean(i) = priorMeans(i);
        betaPriorPrecision(i) = inPrecision(i);
    }
}

exposureModel::exposureModel(exposureModel* tocopy)
{
    nTpt = (tocopy -> nTpt);
    nLoc = (tocopy -> nLoc);

    Eigen::MatrixXd Xc = (tocopy -> X);
    X = Xc;

    Eigen::VectorXd betaPriorPrecisionc = (tocopy -> betaPriorPrecision);
    betaPriorPrecision = betaPriorPrecisionc;

    Eigen::VectorXd betaPriorMeanc = (tocopy -> betaPriorMean);
    betaPriorMean = betaPriorMeanc;

    Eigen::VectorXd offsetc = (tocopy -> offset);
    offset = offsetc; 
}

int exposureModel::getModelComponentType()
{
    return(LSS_EXPOSURE_MODEL_TYPE);
}

void exposureModel::summary()
{
    Rcpp::Rcout << "Exposure Model Summary\n" << 
                   "----------------------" << 
                   "    dim(X): (" << X.rows() << ", " <<
                   X.cols() << ")\n" << 
                   "    beta prior means: (";
    int sz = betaPriorMean.size();
    int i;
    for (i = 0; i < sz; i++)
    {
        if (i != 0) Rcpp::Rcout << "                       ";
        Rcpp::Rcout << betaPriorMean(i);
        if (i+1 < sz)
        {
            Rcpp::Rcout << ",\n";
        }
        else
        {
            Rcpp::Rcout << ")\n";
        }
    }
    Rcpp::Rcout << "    beta prior precision: (";
    for (i = 0; i < sz; i++)
    {
        if (i != 0) Rcpp::Rcout << "                           ";
        Rcpp::Rcout << betaPriorPrecision(i);
        if (i+1 < sz)
        {
            Rcpp::Rcout << ",\n";
        }
        else
        {
            Rcpp::Rcout << ")\n";
        }
    }
    Rcpp::Rcout << "\n";

                   
}

void exposureModel::setOffset(NumericVector offsets)
{
    if (offsets.length() != (nTpt))
    {
        Rcpp::Rcout << "Error: offsets must have length equal to the number of time points.\n";
        Rcpp::stop("Invalid offsets.");
    }
    int i;
    for (i = 0; i < offsets.length(); i++)
    {
        offset(i) = offsets(i);
    }
}

Rcpp::NumericVector exposureModel::getOffset()
{
    Rcpp::NumericVector out(offset.size());
    int i;
    for (i = 0; i < offset.size(); i++)
    {
        out(i) = offset(i);
    }
    return(out);
}

exposureModel::~exposureModel()
{
    if (prot !=0 ){
        Rcpp::stop("can't delete exposureModel, still being used.\n");
    }
}

RCPP_MODULE(mod_exposureModel)
{
    using namespace Rcpp;
    class_<exposureModel>( "exposureModel" )
    .constructor<SEXP,SEXP,SEXP,SEXP,SEXP>()
    .property("offsets", &exposureModel::getOffset, &exposureModel::setOffset);
}


