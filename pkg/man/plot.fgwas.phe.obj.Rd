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
\item{curve.fitting}{boolean value indicating whether the fitted curve is plot.}
\item{file.pdf}{string indicating the output pdf file name.}
}

\value{
No return values, only figure is exported to PDF file. 
}

\details{
The phenotype object is described in \code{\link{fg.load.phenotype}}.
}

\examples{

# data simulation 
r<-fg.simulate( "Logistic", "SAD1", 2000, 1000, c(2,4,6,8,10) );
# plot phenotype traits
plot(r$obj.phe, file.pdf="temp.fg.test.pdf");
# remove the PDF file
unlink("temp.fg.test.pdf");

}

\keyword{Data}
\keyword{Plot}
