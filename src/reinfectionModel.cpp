#include <Rcpp.h>
#include <reinfectionModel.hpp>


using namespace Rcpp;


reinfectionModel::reinfectionModel(SEXP reinfectMode)
{
     Rcpp::IntegerVector modeVec(reinfectMode);
     reinfectionMode = modeVec(0);
     betaPriorPrecision = Eigen::VectorXd(2);
     betaPriorMean = Eigen::VectorXd(2);
     X_rs = Eigen::MatrixXd(1,1); X_rs(0,0) = -1.0;
     betaPriorMean(0) = -1;
     betaPriorMean(1) = -1;
     betaPriorPrecision(0) = -1;
     betaPriorPrecision(1) = -1;
}

reinfectionModel::reinfectionModel(reinfectionModel* tocopy){
    reinfectionMode = (tocopy -> reinfectionMode);

    Eigen::MatrixXd X_rsc = (tocopy -> X_rs); 
    X_rs = X_rsc;

    Eigen::VectorXd betaPriorPrecisionc = (tocopy -> betaPriorPrecision);
    betaPriorPrecision = betaPriorPrecisionc;

    Eigen::VectorXd betaPriorMeanc = (tocopy -> betaPriorMean);
    betaPriorMean = betaPriorMeanc;
}

int reinfectionModel::getModelComponentType()
{
    return(LSS_REINFECTION_MODEL_TYPE);
}

void reinfectionModel::summary()
{
    Rcpp::Rcout << "Reinfection Model Summary\n" << 
                   "-------------------------\n";
    Rcpp::Rcout << "    reinfection mode: ";
    if (reinfectionMode == 1) Rcpp::Rcout << "SEIRS" << "\n";
    if (reinfectionMode == 2) Rcpp::Rcout << "SEIRS (fixed)" << "\n";
    if (reinfectionMode == 3) Rcpp::Rcout << "SEIR" << "\n";

    if (reinfectionMode != 3)
    {
        Rcpp::Rcout << "    dim(X_rs): (" << X_rs.rows() << ", " << X_rs.cols() << ")\n"; 

        Rcpp::Rcout << "    beta prior means: (";
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
    }
    Rcpp::Rcout << "\n";
}

void reinfectionModel::buildReinfectionModel(SEXP _X, SEXP _priorMean, SEXP _prec)
{
    Rcpp::NumericMatrix inX(_X);
    Rcpp::NumericVector inPrecision(_prec);
    Rcpp::NumericVector priorMeans(_priorMean);

    if (reinfectionMode != 3)
    {
        if (priorMeans.length() != inPrecision.length() || inPrecision.length() != inX.ncol())
        {
            Rcpp::stop("Number of parameters, prior means, or precisions does not equal the number of supplied covariates.\n");
        }
        int i,j;
        betaPriorPrecision = Eigen::VectorXd(inX.ncol());
        betaPriorMean = Eigen::VectorXd(inX.ncol());

        for (i = 0; i < inX.ncol(); i++)
        {
            betaPriorPrecision(i) = inPrecision(i);
            betaPriorMean(i) = priorMeans(i);
        }
        X_rs = Eigen::MatrixXd(inX.nrow(), inX.ncol());
        for (i = 0; i < inX.nrow(); i++)
        {
            for (j = 0; j < inX.ncol(); j++)
            {
               X_rs(i,j) = inX(i,j); 
            }
        }
    }
}

reinfectionModel::~reinfectionModel()
{
    if (prot !=0 ){
        Rcpp::stop("can't delete reinfectionModel, still being used.\n");
    }
}

RCPP_MODULE(mod_reinfectionModel)
{
    using namespace Rcpp;
    class_<reinfectionModel>( "reinfectionModel" )
    .constructor<SEXP>()
    .method("buildReinfectionModel", &reinfectionModel::buildReinfectionModel);
}


