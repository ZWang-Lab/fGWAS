\name{print.fgwas.scan.obj}
\alias{print.fgwas.scan.obj}

\title{Printing the significant SNPs or top SNPs.}

\description{ 
Printing the significant SNPs or top SNPs.
} 

\usage{ \method{print}{fgwas.scan.obj}( x,  ..., useS4 = FALSE ) }

\arguments{
\item{x}{a result object return by \code{\link{fg.snpscan}} }
\item{...}{additional arguments affecting the summary produced.}
\item{useS4}{an argument used to match showDefault function. Fixed as FALSE.}
}

\value{
No return values, only output the significant SNPs or top SNPs on the R console.
}

\details{

The result object is described in \code{\link{fg.snpscan}}.

The following example demonstrates the output of print command for a result object.

\preformatted{
== Result from 'fGWAS' method ==
 SNP =  431670 
            INDEX       NAME CHR       POS         MAF NMISS           pv
rs11081728 386746 rs11081728  18  29511663 0.169253731     3 5.413541e-15
rs11051162 301826 rs11051162  12  31074982 0.034090909     6 3.997973e-14
rs7876945  423711  rs7876945  23  13134844 0.292572464    22 5.341771e-14
rs12652390 146007 rs12652390   5 126322102 0.014379868     9 2.143879e-13
rs6029168  407076  rs6029168  20  39158207 0.324880668     2 2.182798e-13
rs7876326  424440  rs7876326  23  23626268 0.470481928    18 3.074301e-13
...

}

}


\examples{
r <- fg.simulate("Logistic", "AR1", 2000, 500, 1:7 );
obj.scan <- fg.snpscan(r$obj.gen, r$obj.phe, method="fast");
obj.scan;

}

\keyword{Scanning}
\keyword{Print}
