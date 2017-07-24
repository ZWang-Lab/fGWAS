\name{plot.fgwas.curve}
\alias{plot.fgwas.curve}
\title{
Plot the curves associated with different genotype 
}
\description{
Plot the curves based on the parameter estimation in GWAS analysis for each genotype.
}
\usage{
plot.fgwas.curve( object, snp.sub, file.pdf = NULL)
}
\arguments{
  \item{object}{Result object returned by the \code{\link{fg.snpscan}}.}
  \item{snp.sub}{Vector of SNP index or SNP name.}
  \item{file.pdf}{String indicating the PDF file name.}
}
\details{
This function will apply to the result obtained from 'fgwas', 'fgwas-optin' or 'fast' methond.
}
\value{
No return values.
}
\examples{

r <- fg.simulate("Logistic", "AR1", 2000, 500, 1:7 );
obj.scan <- fg.snpscan(r$obj.gen, r$obj.phe, method="fast");
tb.sig <- fg.select.sigsnp(obj.scan)
plot.fgwas.curve( obj.scan, tb.sig$INDEX, file.pdf="test.plot.pdf");
unlink("test.plot.pdf")

}
\keyword{ Plot }
