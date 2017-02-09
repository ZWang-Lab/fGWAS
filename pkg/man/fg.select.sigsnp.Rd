\name{fg.select.sigsnp}
\alias{fg.select.sigsnp}
\title{
Selecting significant SNPs
}
\description{
Selecting significant SNPs
}
\usage{
fg.select.sigsnp( fgwas.scan.obj, sig.level=0.05, pv.adjust="bonferroni",  options=list() )
}
\arguments{
  \item{fgwas.scan.obj}{A S3 object of scaning result returned by \code{\link{fg.snpscan}}. }
  \item{sig.level}{Numeric value indicating the p-value of criteria. }
  \item{pv.adjust}{string indicating whether p-value is adjusted value for multiple Comparisons, the optional values are holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" defined in \code{\link{p.adjust}}.}
  \item{options}{None}
}
\details{
}
\value{
The function returns a sub-matrix including the significant SNPs.

}
\examples{
}
\keyword{ fGWAS }
\keyword{ Scaning }
