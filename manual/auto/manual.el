(TeX-add-style-hook "manual"
 (lambda ()
    (LaTeX-add-labels
     "sec:credits"
     "sec:start"
     "synthmode"
     "modelfile"
     "instprof"
     "depcoef"
     "profilefile"
     "invmode"
     "errorbars"
     "sec:tips"
     "sec:restarting"
     "sec:compiling"
     "sec:MPI")
    (TeX-run-style-hooks
     "latex2"
     "bk10"
     "book")))

