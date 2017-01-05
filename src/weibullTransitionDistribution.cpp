#include <transitionDistribution.hpp>
#include <cmath>

double gammapdf(double x, double shape, double rate, bool log)
{
    if (x <= 0)
    {
        return(log ? -std::numeric_limits<double>::infinity() : 0);
    }
    double out = -std::lgamma(shape) + shape*std::log(rate)
            + (shape - 1.0)*std::log(x) - x*rate;
    return (log? out : std::exp(out));
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
    return(gammapdf(params(0), shapePriorAlpha, shapePriorBeta, true)
          +gammapdf(params(1), scalePriorAlpha, scalePriorBeta, true));
}

double weibullTransitionDistribution::getTransitionProb(int startIdx,
                                                        int stopIdx)
{
    return(1.0 - std::exp(std::pow(startIdx/currentScale, currentShape) 
                         - std::pow(stopIdx/currentScale, currentShape) ));
}

double weibullTransitionDistribution::getAvgMembership()
{
    return(currentScale*std::exp(std::lgamma(1 + 1.0/currentShape)));
}

void weibullTransitionDistribution::setCurrentParams(
                                                Eigen::VectorXd currentParams)
{
    currentShape = currentParams(0);
    currentScale = currentParams(1);
}



