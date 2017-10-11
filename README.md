# rdSensitivity

R function for testing the sensitivity of regression discontinuity effect estimates to variations in kernel, bandwidth, and polynomial order. Output is in a kable-type table

Dependencies: rdrobust, knitr

Function arguments: 

y= the outcome of interest, 
x= the forcing variable (variable used to decide who gets treatment/control), c= cutoff value for forcing variable; those above get treatment, those below get control; default is zero
fuzzy= vector indicating who actually recieved treatment, if not all those assigned to treatment (i.e., above cutoff) ended up receiving it and/or if not all those assigned to control (below cutoff) did not recieve treatment. Estimates when fuzzy is specified are marginal complier average treatment effects