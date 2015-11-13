#ifndef SPATIALSEIR_TRANSITION_DISTRIBUTIONS
#define SPATIALSEIR_TRANSITION_DISTRIBUTIONS

#include<Eigen/Core>

class transitionDistribution 
{
    public:
        virtual ~transitionDistribution(){};
        virtual double evalParamPrior(Eigen::VectorXd params){};
        virtual void setCurrentParams(Eigen::VectorXd currentParams){};
        virtual double getTransitionProb(int startIdx, 
                                         int stopIdx){}; 
        virtual double getAvgMembership(){};
};

class weibullTransitionDistribution : public transitionDistribution
{
    public:
        weibullTransitionDistribution(Eigen::VectorXd priorParams);
        double evalParamPrior(Eigen::VectorXd params);
        void setCurrentParams(Eigen::VectorXd currentParams);
        double getTransitionProb(int startIdx, 
                                 int stopIdx);
        double getAvgMembership();
        ~weibullTransitionDistribution();
    private:
        double shapePriorAlpha;
        double shapePriorBeta;
        double scalePriorAlpha;
        double scalePriorBeta;
        double currentShape;
        double currentScale;
};

#endif
