\name{fg.load.simple}
\alias{fg.load.simple}
\title{
Loading simple SNP file 
}
\description{
Loading the SNP data from a text file in simple format.
}
\usage{
fg.load.simple( file.simple.snp, options=list(verbose=F))
}
\arguments{
  \item{file.simple.snp}{string, file name of SNP data table, the format is described in the details. }
  \item{options}{list, only including the item \code{verbose}.}
}
\value{

A S3 object (\code{fgwas.gen.obj}) is returned by this function, including:\cr

\item{options}{List,}
\item{reader}{A reference class \code{}, indicating the reader object for simple data format.}
\item{params}{List, including the file names and other optional parameters.}
\item{n.snp}{Numeric, indicating the total SNP number.}
\item{n.ind.total}{Numeric, indicating the total individual number.}
\item{n.ind.used}{Nuermic, indicating the used individual number except the missing data. }

}
\details{

This table file contains SNP information and individual SNPs. 
The SNP information, including SNP name, chromosome, position, reference allele and alternate allele are located from the 1st column to 5rd column. 
The individual SNPs follow the SNP information at each row, encoding genotype type as 0,1,2 or NA (for missing data). 
For example:

\preformatted{	
SNP  CHR POS RefAllele AltAllele Sub1 Sub2 Sub3 Sub4 Sub5 ...
SNP1 1 1  A T  2  0  1  1 1 ...
SNP2 1 2  A G  2  1  1  0 0 ...
...
}

}

\examples{


# make a simple SNP data table
snp <- matrix( round( runif(20000, 0, 2) ), nrow=100)
colnames(snp) <- paste("sub", 1:200, sep="");
snp.name <- paste("snp", 1:100, sep="");
# make a data frame including snp information and snp data table
snp.mat <- data.frame(snp.name, chr=1, pos=1:100, Allele1="A", Allele2="B", snp);
# write the data frame into a tab-seperated table file
write.table(snp.mat, file="temp.fgwas.snp.mat.csv", 
    row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t");
# load the snp data from a simple snp data file.    
obj.gen <- fg.load.simple("temp.fgwas.snp.mat.csv");

# simulate the phenotype and genotype data and generate a simple snp data file
objx <- fg.simulate( "Logistic", "SAD1", 2000, 800, 1:6, phe.missing=0.05, snp.missing=0.05, 
    sig.pos=301, plink.format=FALSE, file.prefix = "test1.fgwas" );
# load the snp data generated from the simulation
obj.gen <- fg.load.simple("test1.fgwas.geno.tab");
# call SNP scaning 
r1 <- fg.snpscan(obj.gen, objx$obj.phe, snp.sub=c(290:310), options=list(ncores=1));

# simulate the phenotype and genotype data and generate PLINK data file
objx <- fg.simulate( "Logistic", "SAD1", 2000, 800, 1:6, phe.missing=0.05, snp.missing=0.05, 
    sig.pos=301, plink.format=TRUE, file.prefix = "test1.fgwas" );

# load PLINK data file using the function in the snpStats package
obj.plink <- read.plink( "test1.fgwas.geno.bed", "test1.fgwas.geno.bim", "test1.fgwas.geno.fam");
# convert the PLINK data object into SNP data table
x0 <- as.data.frame( obj.plink$map[,c(2,1,4,5,6)])
x1 <-  matrix( as.numeric( (t(obj.plink$genotypes ))), nrow=800)
x1[which(x1==0)]<-NA
x1 <- x1 - 1;
x1 <-  matrix( x1, nrow=800);
colnames(x1) <- paste("N", 1:2000, sep="_");
snp.mat <- data.frame(x0, x1) ;
# write snp data table into a tab-seperated table file.
write.table(snp.mat, file="temp.fgwas.snp.mat.csv", 
    row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t");
    
# load the snp data converted from the PLINK files
obj.gen <- fg.load.simple("temp.fgwas.snp.mat.csv");
# call SNP scaning in a short range(290:310)
r2 <- fg.snpscan(obj.gen, objx$obj.phe, snp.sub=c(290:310), options=list(ncores=1));



}
\keyword{ Data }
