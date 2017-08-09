\name{print.fgwas.gen.obj}
\alias{print.fgwas.gen.obj}

\title{Printing brief information for the genotype object}

\description{ 
  Printing brief information for the genotype object.
} 

\usage{ \method{print}{fgwas.gen.obj}( x, ..., useS4 = FALSE ) }

\arguments{
\item{x}{a genotype object returned by \code{\link{fg.load.plink}} or \code{\link{fg.load.simple}} }
\item{...}{additional arguments affecting the summary produced.}
\item{useS4}{an argument used to match showDefault function. Fixed as FALSE.}
}

\value{
No return values, only output the brief information on the R console.
}

\details{

The genotype object is described in \code{\link{fg.load.plink}} or \code{\link{fg.load.simple}}.

The following example demonstrates the output of \code{print} command for a genotype object.\cr

\preformatted{

== Genotype Object in fGWAS ==
  Plink bed : /demo.bed 
  Plink bim : /demo.bim 
  Plink fam : /demo.fam 
  Data file :   
  SNP count : 431670 
  Total individuals : 1678 
  Reference matrix object of class "fg.dm.plink" 
  Data type:   
  Description:   
  Chromosome:  all 
  SNP Count:  431670 
  Individual Count:  1678 
  Individual Used:  1678 
  Plink Command:  plink 
}

}


\examples{
# data simulation
r<-fg.simulate("Logistic", "AR1", 2000, 500, 1:7 );
# call this function to print the genotype information
r$obj.gen;
# summarize the genotype object
summary(r$obj.gen);

}

\keyword{Data}
\keyword{Print}
