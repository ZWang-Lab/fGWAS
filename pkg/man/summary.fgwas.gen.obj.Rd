\name{summary.fgwas.gen.obj}
\alias{summary.fgwas.gen.obj}

\title{Summarizing genotype object}

\description{ 
  Print the summary information of the genotype object.
} 

\usage{ \method{summary}{fgwas.gen.obj}( object, ... ) }

\arguments{
\item{object}{a data object returned by \code{\link{fg.load.plink}} or \code{\link{fg.load.simple}} }
\item{...}{additional arguments affecting the summary produced.}
}

\value{
A numeric vector including snp count, individual count and used individual count.
}

\details{
The genotype object is described in \code{\link{fg.load.plink}} or \code{\link{fg.load.simple}}.
}
\examples{

r<-fg.simulate("Logistic", "AR1", 2000, 500, 1:7 );
r$obj.gen;
summary(r$obj.gen);

}

\keyword{Data}
\keyword{Summary}
