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

# simulate the phenotype object and genotype object.
r<-fg.simulate("Logistic", "AR1", 2000, 500, 1:7 );
# call SNP scaning using 'fast' method 
obj.scan <- fg.snpscan(r$obj.gen, r$obj.phe, covariance.type="AR1", method="fast");
# select the significant SNPs from the result object
tb<-fg.select.sigsnp(obj.scan, sig.level=1e-10, "fdr")

# show 5 lines of significant SNP table
head(tb);
# show SNP indexes in significant SNP table 
show(tb$INDEX);
# show SNP names in significant SNP table 
show(tb$NAME);
}
\keyword{ Scanning }
