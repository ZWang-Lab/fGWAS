\name{plot.fgwas.phe.obj}
\alias{plot.fgwas.phe.obj}

\title{Ploting the longitudinal curves}

\description{ 
  Ploting the longitudinal curves.
} 

\usage{ \method{plot}{fgwas.phe.obj}(x, y, ..., curve.fitting = T, file.pdf = NULL) }

\arguments{
\item{x}{a phenotype object return by \code{\link{fg.load.phenotype}} or \code{\link{fg.simulate}} }
\item{y}{required parameter in the generic plot function, not used.}
\item{...}{additional arguments affecting the summary produced.}
\item{curve.fitting}{...}
\item{file.pdf}{string indicating the output pdf file name.}
}

\value{
No return values, only figure is exported to PDF file. 
}

\details{
The data object is described in \code{\link{fg.load.phenotype}}.
}

\examples{
r<-fg.simulate( "Logistic", "SAD1", 2000, 1000, c(2,4,6,8,10) );
plot(r$obj.phe, file.pdf="file.test.pdf");
unlink("file.test.pdf");

}

\keyword{Data}
\keyword{Plot}
