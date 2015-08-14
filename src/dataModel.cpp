#include <Rcpp.h>
#include <dataModel.hpp>


using namespace Rcpp;

dataModel::dataModel(SEXP _Y, SEXP type, SEXP compartment, SEXP _phi, 
        SEXP _na_mask)
{
    Rcpp::NumericMatrix input(_Y);
    Rcpp::IntegerMatrix input_na_mask(_na_mask);

    nLoc = input.ncol();
    nTpt = input.nrow();
    Rcpp::NumericVector in_phi(_phi);

    phi = in_phi(0);

    Rcpp::StringVector inputType(type);
    Rcpp::StringVector inputCompartment(compartment);

    if (inputType(0) == "identity")
    {
        dataModelType = 0;
    }
    else if (inputType(0) == "overdispersion")
    {
        dataModelType = 1;
    }
    else
    {
        Rcpp::Rcout << "Unrecognized data model type: " << type << "\n";
        Rcpp::Rcout << "Falling back to default: identity\n";
        dataModelType = 0;
    }

    if (inputCompartment(0) == "I_star")
    {
        dataModelCompartment = 0;
    }
    else if (inputCompartment(0) == "R_star")
    {
        dataModelCompartment = 1;
    }
    else
    {
        Rcpp::Rcout << "Unrecognized data model compartment: " << compartment << "\n";
        Rcpp::Rcout << "Falling back to default: I_star\n";
        dataModelCompartment = 0;
    }
    Y = Eigen::MatrixXi(nTpt, nLoc);
    na_mask = MatrixXb(nTpt, nLoc);
    int i,j;
    for (i = 0; i < (nTpt); i++)
    {
        for (j = 0; j < (nLoc); j++)
        {
            Y(i,j) = input(i, j); 
            na_mask(i, j) = (input_na_mask(i, j) == 1); 
        }
    }
}

int dataModel::getModelComponentType()
{
    return(LSS_DATA_MODEL_TYPE);
}

void dataModel::summary()
{
    Rcpp::Rcout << "Number of locations: " << nLoc << "\n";
    Rcpp::Rcout << "Number of time points: " << nTpt << "\n";
    Rcpp::Rcout << "Data Model Compartment: ";
    if (dataModelCompartment == 0)
    {
        Rcpp::Rcout << "I_star\n";
    }
    else if (dataModelCompartment == 1)
    {
        Rcpp::Rcout << "R_star\n";
    }
    else
    {
        Rcpp::Rcout << "Unknown\n";
    }
    Rcpp::Rcout << "Data Model Type: ";
    if (dataModelType == 0)
    {
        Rcpp::Rcout << "identity\n";
    }
    else if (dataModelType == 1)
    {
        Rcpp::Rcout << "overdispersion\n";
    }
    else
    {
        Rcpp::Rcout << "Unknown\n";
    }

}
dataModel::~dataModel()
{
    if (prot !=0 ){
        Rcpp::stop("can't delete dataModel, still being used.\n");
    }
}

RCPP_MODULE(mod_dataModel)
{
    using namespace Rcpp;
    class_<dataModel>( "dataModel" )
    .constructor<SEXP,SEXP,SEXP,SEXP,SEXP>()
    .method("summary", &dataModel::summary);
}


