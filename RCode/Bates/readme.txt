The code comprised of two parts: (1) scripts for generating figures and tables and (2) method source code.


1.  Figures and table generation.

 * pros_example.Rmd: code for the prostate cancer data analysis.
 * imputation_investigation.Rmd: code for evaluating different zero-replacement schemes for the prostate cancer data.
 * comprehensive_simulation.Rmd: code for the log-ratio lasso MSE and support recovery simulation experiments.
 * runtime_sims.R: code for the runtime simulation experiment.
 * post_selective_inference.R: code for the selective inference example.


2.  Method source code:

Files containing code to run the approximate forward stepwise procedure and choose the optimal number of steps by cross validation.
 * approximate_fs.r
 * cv_fs.r

Files containing the one and two stage log-ratio lasso procedures, with associated CV routines:
 * logRatioLasso.r
 * two_step.r
 * two_step_cv.r

File containing the selective inference code:
 * si.r

Miscellaneous helper functions for the prostate data analysis:
 * myFuns.r

