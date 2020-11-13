## Default repo
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)
})

packages <- c(
  "mvtnorm",
  "balance",
  "selbal",
  "propr",
  "microbenchmark",
  "ggplot2",
  "logratiolasso",
  "foreach",
  "doFuture",
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