\name{print.fgwas.phe.obj}
\alias{print.fgwas.phe.obj}

\title{Printing brief information for the phenotype object}

\description{ 
Printing brief information for the phenotype object.
} 

\usage{ \method{print}{fgwas.phe.obj}( x, ..., useS4=FALSE ) }

\arguments{
\item{x}{a phenotype data object returned by \code{\link{fg.load.phenotype}} or \code{obj.phe} in the \code{\link{fg.simulate}}}
\item{...}{additional arguments affecting the summary produced.}
\item{useS4}{an argument used to match showDefault function. Fixed as FALSE.}
}

\value{
No return values, only output the brief information on the R console.
}

\details{

The phenotype object is described in \code{\link{fg.load.phenotype}}.

The following example demonstrates the output of summary command for a phenotype object.\cr

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

}


\examples{

r<-fg.simulate("Logistic", "AR1", 2000, 500, 1:7 );
r$obj.phe;


}

\keyword{Data}
\keyword{Print}
