\name{fg.load.simple}
\alias{fg.load.simple}
\title{
Loading simple SNP file 
}
\description{
Loading the SNP data from simple text file
}
\usage{
fg.load.simple( file.simple.snp, options=list())
}
\arguments{
  \item{file.simple.snp}{string, CSV file name of SNP data, the format is described in the details. }
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

This CSV file contains SNP information and individual SNPs. 
The SNP information, including SNP name, chromosome and BP, are located from the 1st column to 3rd column. 
The individual SNPs follows the SNP information at each row, encoding genotype type as 0,1,2 or NA (for missing data).

\preformatted{	
CHR,POS,Sub1,Sub2,Sub3,Sub4,Sub5,...
0,1, 2, 0, 1, 1, 1, ...
0,2, 2, 1, 1, 0, 0, ...
...
}

}

\examples{
}
\keyword{ Data }
