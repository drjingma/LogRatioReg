## Default repo
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)
})

# need to do library("igraph"), 
#   but for that, need a fortran compiler, so outside of R, may need to do:
#     sudo apt-get install libpoppler-cpp-dev
#     sudo apt-get install libapparmor-dev
# (note: this might be specific to ubuntu)

packages <- c(
  "mvtnorm",
  "balance", 
  "selbal", 
  "propr", 
  "Matrix", 
  "glmnet", 
  "compositions", 
  "stats", 
  "limSolve", 
  "microbenchmark",
  "ggplot2",
  "logratiolasso",
  "foreach",
  "future", 
  "parallel", 
  "doFuture",
  "rngtools",
  "doRNG"
)

sapply(packages, function (pkgname) {
  status <- find.package(pkgname, quiet=TRUE)
  if (identical(status, character(0))) {
    install.packages(pkgname)  
  }
})

status <- find.package("data.table", quiet=TRUE)
if (identical(status, character(0))) {
  
  packageurl <- "https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.12.2.tar.gz"
  
  install.packages(packageurl, repos=NULL, type="source")
}
