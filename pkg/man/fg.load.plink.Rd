\name{fg.load.plink}
\alias{fg.load.plink}
\title{
Loading PLINK data set.
}
\description{
Loading genotype data from PLINK data set.
}
\usage{
fg.load.plink(file.plink.bed, file.plink.bim, file.plink.fam, 
    plink.command = NULL, 
    chr = NULL, 
    options = list())
}
\arguments{
  \item{file.plink.bed}{string, the name of PLINK bed file, containing the packed binary SNP genotype data}
  \item{file.plink.bim}{string, the name of PLINK bim file, containing the SNP descriptions}
  \item{file.plink.fam}{string, the name of PLINK fam file, containing subject(and, possibly, family) identifiers}
  \item{plink.command}{string, indicating PLINK command path, for the large PLINK data, the package loads partial plink data extracted by the plink command rather than the whole data set. }
  \item{chr}{vector of string, indicating the chromosome number involved to do hypothesis test.}
  \item{options}{list, including \code{force.split} and \code{verbose}}
}

\value{

This function returns a S3 object with the class label of \code{fgwas.gen.obj}, including:

\item{reader}{a reference class \code{fg.dm.plink}, the reader object for plink data, the structure is demostrated in the details. }
\item{n.snp}{integer value, indicating the total SNP number.}
\item{n.ind.total}{integer value, indicating the total individual number.}
\item{n.ind.used}{integer value, indicating the used individual number except the missing data. }
\item{params}{list, including the plink file names and plink command assigned to this function.}
\item{options}{list, two options: \code{force.split} and \code{verbose}, default values are TRUE.}

You can print summary information by the function \code{\link[=print.fgwas.gen.obj]{print}}.

}

\details{

The following conents are exported by \code{\link[=print.fgwas.gen.obj]{print}} command. \cr

\preformatted{
== Genotype Object in fGWAS ==
  Plink bed : /home/zw355/proj/gwas2/bmi-c1c2-qc2.bed 
  Plink bim : /home/zw355/proj/gwas2/bmi-c1c2-qc2.bim 
  Plink fam : /home/zw355/proj/gwas2/bmi-c1c2-qc2.fam 
  Data file :   
  SNP count : 431670 
  Total individuals : 1678 
  Reference matrix object of class "fg.dm.plink" 
  Data type:   
  Description:   
  SNP Count:  431670 
  Individual Count:  1678 
  Individual Used:  1678 
  Plink Command:  plink 
  chromosome:  all 

}

}
\examples{
}
\keyword{ Data }

