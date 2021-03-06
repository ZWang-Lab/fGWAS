\name{plot.fgwas.curve}
\alias{plot.fgwas.curve}
\title{
Plotting the genetic curve effects associated with SNP. 
}
\description{
Plot the genetic curve effects based on the parameter estimation in the fGWAS model for each genotype.
}
\usage{
plot.fgwas.curve( object, snp.sub, file.pdf=NULL, draw.rawdata=TRUE, draw.meanvector=TRUE, ... )
}
\arguments{
  \item{object}{Result object returned by the \code{\link{fg.snpscan}}.}
  \item{snp.sub}{Vector of SNP index or SNP name.}
  \item{file.pdf}{String indicating the PDF file name.}
  \item{draw.rawdata}{Logical value indicating whether the raw phenotype curves are drawn.}
  \item{draw.meanvector}{Logical value indicating whether the mean vectors for each genotype are drawn.}
  \item{...}{additional arguments affecting the plot, including \code{xlab}, \code{ylab}, \code{xlim}, \code{ylim}, \code{title}. }
}
\details{
This function will apply to the result obtained from 'fgwas', 'fgwas-optim' or 'fast' methond.
}
\value{
No return values.
}
\examples{

# simulate phenotype object and genotype object
r <- fg.simulate("Logistic", "AR1", 2000, 500, 1:7, sig.pos=250 );
# SNP scaning
obj.scan <- fg.snpscan(r$obj.gen, r$obj.phe, method="fast", snp.sub=c(245:255) );
# select significant SNPs by default.
tb.sig <- fg.select.sigsnp(obj.scan)
# plot the genteic curves for each significant SNP.
plot.fgwas.curve( obj.scan, tb.sig$INDEX, file.pdf="test.plot.pdf");
# remove the PDF file
unlink("test.plot.pdf")

}
\keyword{ Plot }
