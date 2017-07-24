\name{summary.fgwas.phe.obj}
\alias{summary.fgwas.phe.obj}

\title{Summarizing phenotype object.}

\description{ 
  Summarizing phenotype object.
} 

\usage{ \method{summary}{fgwas.phe.obj}( object, ... ) }

\arguments{
\item{object}{a phenotype object returned by \code{\link{fg.load.phenotype}} or \code{obj.phe} in the \code{\link{fg.simulate}}}
\item{...}{additional arguments affecting the summary produced.}
}

\value{
This function returns an extend table containing the covariate data , the longtudinal curve data  and the measured times. \cr

}

\details{
None
}

\examples{

r<-fg.simulate("Logistic", "AR1", 2000, 500, 1:7 );
r$obj.phe;
ext_phe <- summary(r$obj.phe);
head(ext_phe);

}

\keyword{Data}
\keyword{Summary}
