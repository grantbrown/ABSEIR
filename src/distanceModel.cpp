#include <Rcpp.h>
#include <distanceModel.hpp>

using namespace Rcpp;


distanceModel::distanceModel()
{
    numLocations=-1;
    spatial_prior = Eigen::VectorXd(2);
    dm_list = std::vector<Eigen::MatrixXd>();
    tdm_list = std::vector<std::vector<Eigen::MatrixXd> >();
    tdm_empty = std::vector<int>();
    currentTDistIdx = 0;
}

int distanceModel::getModelComponentType()
{
    return(LSS_DISTANCE_MODEL_TYPE);
}

void distanceModel::setPriorParameters(double _priorAlpha, double _priorBeta)
{
    spatial_prior(0) = _priorAlpha;
    spatial_prior(1) = _priorBeta;
}

void distanceModel::setupTemporalDistanceMatrices(int nTpt)
{
    for (int i = 0; i < nTpt; i++)
    {
        tdm_list.push_back(std::vector<Eigen::MatrixXd>());
        tdm_empty.push_back(1);
    }
}

void distanceModel::addTDistanceMatrix(int tpt, NumericMatrix distMat)
{
    Eigen::MatrixXd new_mat(distMat.nrow(), distMat.ncol());
    int i,j;
    tpt--;
    bool empty = true;
    if (tpt < 0 || tpt >= (int) tdm_list.size())
    {
        Rcpp::stop("Invalid Time Index");
    }
    for (i = 0; i < distMat.nrow(); i++)
    {
        for (j = 0; j < distMat.ncol(); j++)
        {
            if (distMat(i,j) != 0){
                empty = false;
            }
            new_mat(i,j) = distMat(i,j);
        }
    }
    tdm_list[tpt].push_back(new_mat);
    if (!empty)
    {
        tdm_empty[tpt + tdm_list[tpt].size()] = 0;
    }
}

void distanceModel::addDistanceMatrix(NumericMatrix distMat)
{
    if (distMat.nrow() != distMat.ncol())
    {
        Rcpp::stop("Distance matrix must be square.\n");
    }
    else if (numLocations != -1 && distMat.nrow() != (numLocations))
    {
        Rcpp::stop("Dimension does not match previously added distance matrix\n");
    }
    Eigen::MatrixXd new_mat(distMat.nrow(), distMat.ncol());
    int i,j;
    for (i = 0; i < distMat.nrow(); i++)
    {
        for (j = 0; j < distMat.ncol(); j++)
        {
            new_mat(i,j) = distMat(i,j);
        }
    }

    dm_list.push_back(new_mat);
    numLocations = distMat.nrow();

}

void distanceModel::summary()
{
    Rcpp::Rcout << "Distance Model Summary\n" <<
                   "----------------------\n";
    if (numLocations > 1)
    {
        Rcpp::Rcout << "Number of locations: " << numLocations << "\n";
        Rcpp::Rcout << "Number of distance structures: " << (dm_list.size()) << "\n";
        Rcpp::Rcout << "Number of time varying distance structures: " 
            << (tdm_list.size()) << "\n";
        if (tdm_list.size() > 0)
        {
            Rcpp::Rcout << "Number of time varying distance lags: " 
                << tdm_list[0].size();
        }
    }
    else
    {
        Rcpp::Rcout << "    no distance structure.\n";
    }
    Rcpp::Rcout << "\n";
}

distanceModel::~distanceModel()
{
    if (prot !=0 ){
        Rcpp::stop("can't delete distanceModel, still being used.\n");
    }
}

int distanceModel::getNumDistanceMatrices()
{
    return(dm_list.size()); 
}

RCPP_MODULE(mod_distanceModel)
{
    using namespace Rcpp;
    class_<distanceModel>( "distanceModel" )
    .constructor()
    .method("addDistanceMatrix", &distanceModel::addDistanceMatrix)
    .method("addTDistanceMatrix", &distanceModel::addTDistanceMatrix)
    .method("setupTemporalDistanceMatrices", 
            &distanceModel::setupTemporalDistanceMatrices)
    .method("summary", &distanceModel::summary)
    .method("setPriorParameters", &distanceModel::setPriorParameters)
    .property("numMatrices", &distanceModel::getNumDistanceMatrices, "Number of distict distance matrices.");
}


