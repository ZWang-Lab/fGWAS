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

A S4 object (\code{fg.dm.plink}) is returned by this function, including:\cr

\item{reader}{a reference class \code{fg.dm.plink}, the reader object for plink data, the structure is demostrated in the \link{details} }
\item{n.snp}{integer value, indicating the total SNP number.}
\item{n.ind.total}{integer value, indicating the total individual number.}
\item{n.ind.used}{integer value, indicating the used individual number except the missing data. }
\item{params}{list, including the plink file names and plink command assigned to this function.}
\item{options}{list, two options: \code{force.split} and \code{verbose}, default values are TRUE.}
}

\details{


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


The following section shows an example of plink data object loaded by the function \code{fg.load.plink}.

\preformatted{
# list of 5
# $ reader     :Reference class 'fg.dm.plink' [package "fGWAS"] with 13 fields
#  ..$ type          : chr(0)
#  ..$ description   : chr(0)
#  ..$ n.snp         : int 274
#  ..$ n.ind.total   : int 1678
#  ..$ n.ind.used    : int 1678
#  ..$ ids.used      : chr [1:1678] "8228" "2294" "7395" "309" ...
#  ..$ file.plink.bed: chr "/home/zw355/proj/gwas2/bmi-c1c2-qc2.bed"
#  ..$ file.plink.bim: chr "/home/zw355/proj/gwas2/bmi-c1c2-qc2.bim"
#  ..$ file.plink.fam: chr "/home/zw355/proj/gwas2/bmi-c1c2-qc2.fam"
#  ..$ plink.command : chr "plink"
#  ..$ chromosome    : num 0
#  ..$ snp.blocksize : num 1000
#  ..$ snpdata       :List of 2
#  .. ..$ plink.bim:'data.frame':        274 obs. of  6 variables:
#  .. .. ..$ V1: int [1:274] 0 0 0 0 0 0 0 0 0 0 ...
#  .. .. ..$ V2: chr [1:274] "ss66037954" "ss66040774" "ss66043486" "ss66047046" ...
#  .. .. ..$ V3: int [1:274] 0 0 0 0 0 0 0 0 0 0 ...
#  .. .. ..$ V4: int [1:274] 0 0 0 0 0 0 0 0 0 0 ...
#  .. .. ..$ V5: chr [1:274] "A" "T" "C" "G" ...
#  .. .. ..$ V6: chr [1:274] "G" "C" "G" "C" ...
#  .. ..$ plink.fam:'data.frame':        1678 obs. of  6 variables:
#  .. .. ..$ V1: chr [1:1678] "1" "2" "8" "17" ...
#  .. .. ..$ V2: int [1:1678] 8228 2294 7395 309 18679 392 7079 19621 22156 22041 ...
#  .. .. ..$ V3: int [1:1678] 0 0 0 0 0 0 0 0 0 0 ...
#  .. .. ..$ V4: int [1:1678] 0 0 0 0 0 0 0 0 0 0 ...
#  .. .. ..$ V5: int [1:1678] 2 1 2 2 2 1 2 2 1 1 ...
#  .. .. ..$ V6: num [1:1678] 20.8 31.4 19.9 25.8 30.3 ...
#  ..and 26 methods, of which 12 are  possibly relevant:
#  ..  get_individuals, get_individuals#fg.dm.genodata, get_snpinfo,
#  ..  get_snpinfo#fg.dm.genodata, get_snpmat, get_snpmat#fg.dm.genodata,
#  ..  select_individuals, select_individuals#fg.dm.genodata, show#envRefClass,
#  ..  show#fg.dm.genodata, shrink, shrink#fg.dm.genodata
# $ n.snp      : int 274
# $ n.ind.total: int 1678
# $ n.ind.used:  int 1678
# $ options    :List of 2
#  ..$ force.split: logi TRUE
#  ..$ verbose    : logi TRUE
# $ params     :List of 4
#  ..$ file.plink.bed: chr "/home/zw355/proj/gwas2/bmi-c1c2-qc2.bed"
#  ..$ file.plink.bim: chr "/home/zw355/proj/gwas2/bmi-c1c2-qc2.bim"
#  ..$ file.plink.fam: chr "/home/zw355/proj/gwas2/bmi-c1c2-qc2.fam"
#  ..$ plink.command : NULL
# - attr(*, "class")= chr "fgwas.gen.obj"
}


}
\examples{
}
\keyword{ fGWAS }
\keyword{ Loading }
