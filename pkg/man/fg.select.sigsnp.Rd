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
  \item{options}{Not used}
}
\details{

}
\value{
The function returns a sub-matrix including the significant SNPs.
}
\examples{

r<-fg.simulate("Logistic", "AR1", 2000, 500, 1:7 );
obj.scan <- fg.snpscan(r$obj.gen, r$obj.phe, covariance.type="AR1", method="fast");
tb<-fg.select.sigsnp(obj.scan, sig.level=1e-10, "fdr")
head(tb);
show(tb$INDEX);
show(tb$NAME);
}
\keyword{ Scanning }
