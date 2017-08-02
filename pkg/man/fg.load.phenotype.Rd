\name{fg.load.phenotype}
\alias{fg.load.phenotype}
\title{
Loading phenotypic data
}
\description{
Loading longitudinal phenotypic curve and covariate files.
}
\usage{
fg.load.phenotype(file.phe.long, file.phe.cov=NULL, file.phe.time=NULL, 
    curve.type = NULL, 
    covariance.type = NULL, 
    file.plot.pdf = NULL, 
    intercept = TRUE,
    options = list(verbose=F))
}
\arguments{
  \item{file.phe.long}{string, indicating the file name of longitudinal curve for each individual.}
  \item{file.phe.cov}{string, indicating the file name of covariates.}
  \item{file.phe.time}{string, indicating the file name of observed time points for each individual. If this file is not specified, the column indexes are used as measured times. }
  \item{curve.type}{string, indicating the curve type, available options are in the details. Default is 'auto' which means the package selects the curve type using curve fitting.}
  \item{covariance.type}{string, indicating the covariance structure type, available options are in the details. Default is 'auto' which means the package selects the covariance type using MLE.}
  \item{intercept}{boolean, indicating whether intercept is included in the model.}
  \item{file.plot.pdf}{string, indicating a PDF file name to illustrate the performance of curve fitting.}
  \item{options}{list, including \code{max.optim.failure},  \code{min.optim.success},  \code{R2.loop}, and \code{verbose}}
}
\value{

This function returns a S3 object with the class label of \code{fgwas.phe.obj}:

\item{pheY}{Matrix, the longitidual curve data with the rowname indicating the individuals' id.}
\item{pheX}{Matrix, the covariate with the rowname indicating the individuals' id.}
\item{pheT}{Matrix, the observed time points with the rowname indicating the individuals' id.}
\item{ids}{Vector, the common individuals' id in all data file}
\item{obj.curve}{Curve object, S4 Object inherited from \code{fg.curve.base}  }
\item{obj.covar}{Covariance object, S4 Object inherited from \code{fg.covariance.base} }
\item{est.curve}{List, including the curve type and estimated parameters.}
\item{est.covar}{List, including the covariance type and estimated parameters. }
\item{summary.curve}{List, including the result of curve fitting.}
\item{summary.covar}{List, including the result of covariance fitting.}
\item{params}{List, the parameters of function calling, including the file names of longitudinal data and covariate, curve type and covariance type}
\item{options}{List, not used currently }

You can plot the curve by the function \code{\link[=plot.fgwas.phe.obj]{plot}} or print summary information by the function \code{\link[=print.fgwas.phe.obj]{print}}.

}

\details{


The phenotype file(\code{file.phe.long}), observed time file(\code{file.phe.time}) and covariate file(\code{file.phe.cov}) must be a CSV file.
The following sections illustrate the format of each data file.

1) The phenotype file.
The first column is individual ID and the rest columns are sample data for each measurement. It
looks like the following file. Please note missing data is coded as space or NA in all data files.

\preformatted{
ID. 1st 2nd 3rd 4th 5th 6th 7th
1 2.9033 4.118 6.1495 7.8161 9.8379 12.963 14.918
2 4.3306 5.3783 7.0647 9.3624 11.439 NA 15.701
3 2.3997 4.052 5.5431 7.6933 9.8471 NA 12.849
4 3.3044 4.154 5.8924 7.7133 9.2144 10.945 NA
...
}

2) The measurement time file. The first column is individual ID and the rest columns are observed times as the following format.

\preformatted{
ID. 1st 2nd 3rd 4th 5th 6th 7th
1 1 2 3 4 5 6 7
2 1 2 3 4 5 NA 7
3 1 2 3 4 5 NA 7
4 1 2 3 4 5 6 NA
...
}

3) The covariate file. The first column is individual ID and the rest columns are covariate values as the following format.

\preformatted{
ID. X1 X2
1 1.0 0 
2 1.1 1 
3 2.1 1 
4 3.1 0 
...
}


The function returns a S3 object which details can be checked by the command \code{print} or \code{str}. \cr

The following conents are exported by \code{\link[=print.fgwas.phe.obj]{print}} command. \cr

\preformatted{
== Phenotype Object in fGWAS ==
  Longitudal value : /tmp/Rtmp0lhMx2/file58b1776d8395.csv
   -- Individual count : 1678
  Longitudal time : /tmp/Rtmp0lhMx2/file58b1131b9b1.csv
   -- Time count : 8
  Covariate file : /tmp/Rtmp0lhMx2/file58b1246a670b.csv
   -- Covariate count : 6
   -- Intercept : YES
   -- Estimate values : 36917.16 -1.553526 -2.702489 7.852712 -5.91549 7.210022 3.805883
  Curve type : auto
   -- Estimate type : Legendre2
   -- Estimate values : -36888.65 1.081318 -1.852532
  Covariate type : auto
   -- Estimated type : TOEPH
   -- Estimate values : 0.9082402 0.8627414 0.8142826 0.7713539 0.7200614 0.6669603 ...
}



9 curves have been implemented in current version, including: \cr

\bold{1) "Logistic" }
        \deqn{g(t) = \frac{a}{1+b*e^{-r*t}} }{g(t)=a/(1+b*exp(-rt))}


\bold{2) "Bi-Logistic"}
        \deqn{g(t) = \frac{a1}{1+b1*e^{-r1*t}} + \frac{a2}{1+b2*e^{-r2*t}} }{ g(t)=a1/(1+b1*exp(-r1*t)) + a2/(1+b2*exp(-r2*t)) }

\bold{3) "Pharmacology"}
        \deqn{g(t) = \frac{E_{max}*t}{Ec_{50}+t} + E_{0} }{ g(t)=Emax*t/(Ec50+t) + E0 }
        
\bold{4) "Exponential"}
        \deqn{g(t)= a*e^{-r*t} }{ g(t)=a1*exp(-r*t) }

\bold{5) "Bi-Exponential"}
        \deqn{g(t) = a_{1}*e^{-r_{1}*t} + a_{2}*e^{-r_{2}*t} }{ g(t)=a1*exp(-r1*t) + a2*exp(-r2*t) }

\bold{6) "Power"}
        \deqn{g(t)= a*t^b }{ g(t)=a*t^b }

\bold{7) "Legendre2"}, Legendre Polynomial(2nd-order)
        \deqn{ g(t) = u_{0}+u_{1}*t + u_{2}*(3*t^2-1)/2 }{ y = u0 + u1*t + u2*1/2*(3*t^2-1) }

\bold{8) "Legendre3"}, Legendre Polynomial(3rd-order)
        \deqn{ g(t) = u_{0}+u_{1}*t + u_{2}*(2*t^2-1)/2 + u_{3}*(5*t^3-3t)/2 }{ y = u0 + u1*t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) }

\bold{9) "Legendre4"}, Legendre Polynomial(4th-order)
        \deqn{ g(t) = u_{0}+u_{1}*t + u_{2}*(2*t^2+1)/2 + u_{3}*(5*t^3-3t)/2 + ... }{ y = u0 + u1*t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) + u4*1/8*(35*ti^4-30*ti^2+3) }

}

\examples{

objx <- fg.simulate( "Logistic", "SAD1", 2000, 800, 1:6, phe.missing=0.05, snp.missing=0.05, sig.pos=301, plink.format=TRUE, file.prefix = "test1.fgwas" );

obj.phe <- fg.load.phenotype("test1.fgwas.pheY.csv", NULL, "test1.fgwas.pheT.csv", curve.type = "Logistic", "SAD1", file.plot.pdf = NULL, intercept = FALSE,  options = list(verbose=TRUE))

obj.phe;

plot( obj.phe, curve.fitting=FALSE, file.pdf = "test1.fgwas.phe.pdf");

#Remove all test files
unlink("test1.fgwas.*");

}
\keyword{ Data }
