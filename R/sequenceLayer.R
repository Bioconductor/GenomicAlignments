### =========================================================================
### sequenceLayer()
### -------------------------------------------------------------------------


sequenceLayer <- function(x, cigar, from="query", to="reference",
                          D.letter="-", N.letter=".",
                          I.letter="-", S.letter="+", H.letter="+")
{
    call_new_fun_in_cigarillo("sequenceLayer", "project_sequences",
                              x=x, cigars=cigar,
                              from=from, to=to,
                              I.letter=I.letter, D.letter=D.letter,
                              N.letter=N.letter,
                              S.letter=S.letter, H.letter=H.letter)
}

