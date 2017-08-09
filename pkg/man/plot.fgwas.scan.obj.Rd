\name{plot.fgwas.scan.obj}
\alias{plot.fgwas.scan.obj}

\title{Drawing Manhattan plot}

\description{ 
  Drawing p-values in Manhattan plot.
} 

\usage{ \method{plot}{fgwas.scan.obj}(x, y, ..., file.pdf = NULL) }

\arguments{
\item{x}{a result object return by \code{\link{fg.snpscan}} }
\item{y}{required parameter in the generic plot function, but unused in this function.}
\item{...}{additional arguments affecting the summary produced.}
\item{file.pdf}{string indicating output PDF file name.}
}

\value{
No return values, only figure is exported into PDF file.
}

\details{
N/A
}
\examples{
# data simulation
r<-fg.simulate("Logistic", "AR1", 2000, 500, 1:7 );
# SNP scaning using 'fast' method
obj.scan <- fg.snpscan(r$obj.gen, r$obj.phe, method="fast");
# plot the Manhattan figure
plot(obj.scan, file.pdf="temp.fgwas.plot.pdf");
# remove the PDF file
unlink("temp.fgwas.plot.pdf");
}

\keyword{Scanning}
\keyword{Plot}
