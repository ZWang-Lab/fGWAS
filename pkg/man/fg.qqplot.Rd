\name{fg.qqplot}
\alias{fg.qqplot}
\title{
Drawing QQ plot for the fGWAS results.
}
\description{
Drawing QQ plot for the fGWAS results.
}
\usage{
fg.qqplot(object, png.file, title = "", width = 480)
}
\arguments{
  \item{object}{A S3 object of scaning result returned by \code{\link{fg.snpscan}}. }
  \item{png.file}{string indicating PNG file name.}
  \item{title}{string indicating the title in the figure.}
  \item{width}{integer indicating the width in pixel unit.}
}
\details{
}
\value{
No return values. PNG file is exported if this function calling is successful.
}

\seealso{
\code{\link{fg.snpscan}}.
}
\examples{

# simulate phenotype object and genotype object
r <- fg.simulate("Logistic", "AR1", 2000, 500, 1:7, sig.pos=250 );
# SNP scaning
obj.scan <- fg.snpscan(r$obj.gen, r$obj.phe, method="fgwas");

fg.qqplot( obj.scan, "fgwas-pkg-test.png", title="Simulation");

}

\keyword{Plot}
