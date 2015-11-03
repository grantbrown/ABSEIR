#include<transitionDistribution.hpp>
#include <cmath>

double gammapdf(double x, double shape, double scale, bool log)
{
    if (log)
    {
        return(
            -std::lgamma(shape) - shape*std::log(scale)
            + (shape - 1.0)*std::log(x) - x/scale);
    }
    return(
        1.0/(std::tgamma(shape)*std::pow(scale, shape))
        *std::pow(x,shape-1.0)
        *std::exp(-1.0*x/scale));
}

weibullTransitionDistribution::weibullTransitionDistribution(
        Eigen::VectorXd priorParams)
{
    shapePriorAlpha = priorParams(0);
    shapePriorBeta  = priorParams(1);
    scalePriorAlpha = priorParams(2);
    scalePriorBeta  = priorParams(3);
    currentShape = 1.0;
    currentScale = 1.0;
}

weibullTransitionDistribution::~weibullTransitionDistribution()
{
}

double weibullTransitionDistribution::evalParamPrior(Eigen::VectorXd params)
{
    return(gammapdf(params(0), shapePriorAlpha, shapePriorBeta, false)
          *gammapdf(params(1), scalePriorAlpha, scalePriorBeta, false));
}

double weibullTransitionDistribution::getTransitionProb(int startIdx,
                                                        int stopIdx)
{
    return(1.0 - std::exp(std::pow(-stopIdx/currentScale, currentShape) 
                        + std::pow(startIdx/currentScale, currentShape) ));
}

void weibullTransitionDistribution::setCurrentParams(
                                                Eigen::VectorXd currentParams)
{
    currentShape = currentParams(0);
    currentScale = currentParams(1);
}



