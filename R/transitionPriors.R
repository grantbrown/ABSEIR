# transitionPriors module helper functions
TransitionPriors = function(p_ei, p_ir, p_ei_ess, p_ir_ess)
{
    structure(list(mode="prob",
                   p_ei=p_ei,
                   p_ir=p_ir,
                   p_ei_ess=p_ei_ess,
                   p_ir_ess=p_ir_ess,
                   priorAlpha_gammaEI=NA,
                   priorBeta_gamma_EI=NA,
                   priorAlpha_gammaIR=NA,
                   priorBeta_gamma_IR=NA), class = "TransitionPriors")
}

