\name{fg.load.simple}
\alias{fg.load.simple}
\title{
Loading simple SNP file 
}
\description{
Loading the SNP data from very simple text file
}
\usage{
fg.load.simple( file.simple.snp, options=list())
}
\arguments{
  \item{file.simple.snp}{string, text file name of SNP data. This text file encodes genotype type as 0,1,2  and SNP information in the 1st column to 3rd column ( SNP name, chromosome and BP). The other column names are the individuals ID.}
  \item{options}{List}
}
\value{

A S3 object (\code{fgwas.gen.obj}) is returned by this function, including:\cr

\item{options}{List,}
\item{reader}{A reference class \code{}, indicating the plink reader object.}
\item{params}{List, including the plink file names and plink command assigned to this function.}
\item{n.snp}{Numeric, indicating the total SNP number.}
\item{n.ind.total}{Numeric, indicating the total individual number.}
\item{n.ind.used}{Nuermic, indicating the used individual number except the missing data. }

}
\details{
}

\examples{
}
\keyword{ fGWAS }
\keyword{ Loading }
