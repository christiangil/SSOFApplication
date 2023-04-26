SSOFA/NEID
========
You may want to run
```julia
using Pkg
Pkg.activate("NEID")
Pkg.develop(;path=".")
```
from the SSOFA directory to make these scripts work as intended

This folder has scripts which analyze NEID data with SSOF
- init.jl: Reformats the data into an SSOF-acceptable form
- analysis.jl : Performs the model fitting
