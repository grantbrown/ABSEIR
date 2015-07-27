# initialValueContainer module helper function
InitialValueContainer = function(S0, E0, I0, R0) 
{
    # get rid of any data frame nonsense
    structure(list(S0 = as.numeric(S0),
                   E0 = as.numeric(E0),
                   I0 = as.numeric(I0),
                   R0 = as.numeric(R0)), class="InitialValueContainer")
}



