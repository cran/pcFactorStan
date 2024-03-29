# Generated by rstantools.  Do not edit by hand.

# names of stan models
stanmodels <- c("correlation", "correlation_ll", "factor", "factor1", "factor1_ll", "factor_ll", "unidim", "unidim_adapt", "unidim_ll")

# load each stan module
Rcpp::loadModule("stan_fit4correlation_mod", what = TRUE)
Rcpp::loadModule("stan_fit4correlation_ll_mod", what = TRUE)
Rcpp::loadModule("stan_fit4factor_mod", what = TRUE)
Rcpp::loadModule("stan_fit4factor1_mod", what = TRUE)
Rcpp::loadModule("stan_fit4factor1_ll_mod", what = TRUE)
Rcpp::loadModule("stan_fit4factor_ll_mod", what = TRUE)
Rcpp::loadModule("stan_fit4unidim_mod", what = TRUE)
Rcpp::loadModule("stan_fit4unidim_adapt_mod", what = TRUE)
Rcpp::loadModule("stan_fit4unidim_ll_mod", what = TRUE)

# instantiate each stanmodel object
stanmodels <- sapply(stanmodels, function(model_name) {
  # create C++ code for stan model
  stan_file <- if(dir.exists("stan")) "stan" else file.path("inst", "stan")
  stan_file <- file.path(stan_file, paste0(model_name, ".stan"))
  stanfit <- rstan::stanc_builder(stan_file,
                                  allow_undefined = TRUE,
                                  obfuscate_model_name = FALSE)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name,
                            model_cppcode = stanfit$cppcode)
  # create stanmodel object
  methods::new(Class = "stanmodel",
               model_name = stanfit$model_name,
               model_code = stanfit$model_code,
               model_cpp = stanfit$model_cpp,
               mk_cppmodule = function(x) get(paste0("rstantools_model_", model_name)))
})
