\name{summary.fgwas.scan.obj}
\alias{summary.fgwas.scan.obj}

\title{Summary of the data object}

\description{ 
  Summarize information for the data object.
} 

\usage{ \method{summary}{fgwas.scan.obj}( object, ... ) }

\arguments{
\item{object}{a data object return by \code{\link{fg.snpscan}} }
\item{...}{additional arguments affecting the summary produced.}
}

\value{
No return values, only output the summary information on the R console.
}

\details{

The data object is described in \code{\link{fg.snpscan}}.

The following example demonstrates summary command for a data object.

\preformatted{
Data set for FunMap model:
------------------------------------
          Date: 2010-03-19 03:49:36
          Model: Logistic Curve
          Cross: F2
    Pheno. file: simu.pheno.LC.F2
     Geno. file: simu.geno.LC.F2
     Maker file: simu.marker.LC.F2
    Sample size: 100       
   Sample times: 7         
   Marker count: 6         
          LC  a: 19.83678  
              b: 8.96370   
              r: 0.47202   
            rho: 0.75430   
         sigma2: 0.58849   
------------------------------------

}

}


\examples{
#check the codes in fg.snpscan()
}

\keyword{Data}
\keyword{Summary}
