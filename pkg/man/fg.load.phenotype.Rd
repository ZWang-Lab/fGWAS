\name{fg.load.phenotype}
\alias{fg.load.phenotype}
\title{
Loading phenotypic data
}
\description{
Loading longitudinal phenotypic curve and covariate files.
}
\usage{
fg.load.phenotype(file.phe.long, file.phe.cov, file.phe.time, 
    curve.type = NULL, 
    covariance.type = NULL, 
    file.plot.pdf = NULL, 
    options = list())
}
\arguments{
  \item{file.phe.long}{string, indicating the file name of longitudinal curve for each individual.}
  \item{file.phe.cov}{string, indicating the file name of covariates.}
  \item{file.phe.time}{string, indicating the file name of observed time points for each individual.}
  \item{curve.type}{string, indicating the curve type, available options are "auto", "Logistic", "DoubleLogistic" and "ABRK". The full list in the details. Default is 'auto' which means the package selects the curve type using curve fitting.}
  \item{covariance.type}{string, indicating the covariance structure type, available options are "auto", "SAD1" and "AR1". The full list in the details. Default is 'auto' which means the package selects the covariance type using MLE.}
  \item{file.plot.pdf}{string, indicating a PDF file name to illustrate the performance of curve fitting.}
  \item{options}{list, not used currently }
}
\value{

This function returns a S3 object with the class label of \code{fgwas.phe.obj} which structure is identical with the result from \code{\link{fg.simulate}}. 
The phenotype curve can be plotted by the function \code{\link[=plot.fgwas.phe.obj]{plot}}.

\item{pheY}{Matrix, the longitidual curve data with the rowname indicating the individuals' id.}
\item{pheX}{Matrix, the covariate with the rowname indicating the individuals' id.}
\item{pheT}{Matrix, the observed time points with the rowname indicating the individuals' id.}
\item{ids}{Vector, the common individuals' id in all data file}
\item{obj.curve}{Curve object, S4 Object inherited from \code{fg.curve.base}  }
\item{obj.covar}{Covariance object, S4 Object inherited from \code{fg.covariance.base} }
\item{est.curve}{List, including the curve type and estimated parameters.}
\item{est.covar}{List, including the covariance type and estimated parameters. }
\item{params}{List, the parameters of function calling, including the file names of longitudinal data and covariate, curve type and covariance type}
\item{options}{List}

}

\details{


The phenotype file(\code{file.phe.long}), observed time file(\code{file.phe.time}), covariate file(\code{file.phe.cov}) must be a CSV file.
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

3) The covariate file. The first column is individual ID and the rest columns are covariate values as the followingformat.

\preformatted{
ID. X1 X2
1 1.0 0 
2 1.1 1 
3 2.1 1 
4 3.1 0 
...
}


\preformatted{
== Phenotype Object in fGWAS ==
  Longitudal value : /tmp/RtmpBp6LEy/fileb7112f1c3957.csv 
   -- Individual count : 1678 
  Longitudal time : /tmp/RtmpBp6LEy/fileb71152bf4cee.csv 
   -- Time count : 8 
  Covariate file : /tmp/RtmpBp6LEy/fileb7117dadd1cb.csv 
   -- Covariate count : 6 
  Curve type : auto 
   -- Estimate type : Legendre2 
   -- Estimate values : 31.69093 1.081318 -1.852532 
  Covariate type : auto 
   -- Estimate type : TOEPH 
   -- Estimate values : 0.9082581 0.8627943 0.8143804 0.7714791 0.7202178 0.6671215 0.6295313 16.17749 18.59827 18.51997 19.13514 20.15374 20.10135 20.20764 22.55035 
}

The function returns a S3 object, the details can be checked by the command \code{print} or \code{str}. The following conents are exported by \code{str} command

\preformatted{
#List of 10
# $ options  : list()
# $ params   :List of 5
#  ..$ file.phe.long  : chr "/tmp/RtmpYrkzlf/fileda14176b9c0d.csv"
#  ..$ file.phe.cov   : chr "/tmp/RtmpYrkzlf/fileda147a3872df.csv"
#  ..$ file.phe.time  : chr "/tmp/RtmpYrkzlf/fileda148740849.csv"
#  ..$ curve.type     : chr "auto"
#  ..$ covariance.type: chr "auto"
# $ pheY     : num [1:1678, 1:8] 20.5 28.5 19.5 26.1 28.4 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:1678] "8228" "2294" "7395" "309" ...
#  .. ..$ : chr [1:8] "Y_1" "Y_2" "Y_3" "Y_4" ...
# $ pheT     : int [1:1678, 1:8] 40 44 33 42 47 51 28 29 57 45 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:1678] "8228" "2294" "7395" "309" ...
#  .. ..$ : chr [1:8] "Z_1" "Z_2" "Z_3" "Z_4" ...
# $ pheX     : num [1:1678, 1:6] 2 1 2 2 2 1 2 2 1 1 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:1678] "8228" "2294" "7395" "309" ...
#  .. ..$ : chr [1:6] "X_1" "X_2" "X_3" "X_4" ...
# $ ids      : chr [1:1678] "8228" "2294" "7395" "309" ...
# $ obj.curve:Formal class 'fg.curve.log' [package "fGWAS"] with 2 slots
#  .. ..@ type       : chr "Logistic"
#  .. ..@ description: chr "logistic curve"
# $ obj.covar:Formal class 'fg.covariance.AR1' [package "fGWAS"] with 3 slots
#  .. ..@ par_num    : int(0)
#  .. ..@ type       : chr "AR1"
#  .. ..@ description: chr "AR1"
# $ est.curve:List of 7
#  ..$ type       : chr "Logistic"
#  ..$ param      : Named num [1:3] 36 39.8 29.5
#  .. ..- attr(*, "names")= chr [1:3] "" "" ""
#  ..$ param.lower: num [1:3] 34.1 35.9 26.6
#  ..$ param.upper: num [1:3] 37.8 43.4 32.4
#  ..$ parx       : Named num [1:7] -6.93 -1.59 -1.34 8.28 -5.06 ...
#  .. ..- attr(*, "names")= chr [1:7] "" "X_1" "X_2" "X_3" ...
#  ..$ parX.lower : num [1:7] -8.8828 -1.9554 -9.0214 0.0122 -12.7226 ...
#  ..$ parX.upper : num [1:7] -5.1 -1.29 13.85 15.07 1.86 ...
# $ est.covar:List of 2
#  ..$ type : chr "AR1"
#  ..$ param: num [1:2] 0.907 19.916
#
#
# - attr(*, "class")= chr "fgwas.phe.obj"
}


9 curves have been implemented in current version, including:

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
}
\keyword{ fGWAS }
\keyword{ Data }
