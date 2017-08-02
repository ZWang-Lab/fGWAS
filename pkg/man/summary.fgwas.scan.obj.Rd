\name{summary.fgwas.scan.obj}
\alias{summary.fgwas.scan.obj}

\title{Summarizing the result object}

\description{ 
 Summarizing the result object.
} 

\usage{ \method{summary}{fgwas.scan.obj}( object, ... ) }

\arguments{
\item{object}{a result object return by \code{\link{fg.snpscan}} }
\item{...}{additional arguments affecting the summary produced.}
}

\value{
This function returns a data frame including the results of scanned SNPs and excluding covariate and curve object. \cr
}

\details{

Below lists all columns for the analysis with curve 'Legendre2'  and covariance 'AR1':\cr

\preformatted{

#obj.fgwas is result object of 'fgwas' method.
> colnames(summary(obj.fgwas))
 [1] "INDEX"      "NAME"       "CHR"        "POS"        "MAF"       
 [6] "NMISS"      "SNP0"       "SNP1"       "SNP2"       "GENO"      
[11] "LR"         "pv"         "L0"         "h0.X0"      "h0.X1"     
[16] "h0.X2"      "h0.X3"      "h0.X4"      "h0.X5"      "h0.X6"     
[21] "h0.u0"      "h0.u1"      "h0.u2"      "h0X.rho"    "h0X.sigma2"
[26] "L1"         "h1.X0"      "h1.X1"      "h1.X2"      "h1.X3"     
[31] "h1.X4"      "h1.X5"      "h1.X6"      "h1.G0.u0"   "h1.G0.u1"  
[36] "h1.G0.u2"   "h1.G1.u0"   "h1.G1.u1"   "h1.G1.u2"   "h1.G2.u0"  
[41] "h1.G2.u1"   "h1.G2.u2"   "h1X.rho"    "h1X.sigma2"
}

}
\examples{
r<-fg.simulate("Logistic", "AR1", 2000, 500, 1:7, sig.pos=190 );
obj.scan <- fg.snpscan(r$obj.gen, r$obj.phe, covariance.type="AR1", method="fast");
tb <- summary(obj.scan);
head(tb);


}

\keyword{Scanning}
\keyword{Summary}
