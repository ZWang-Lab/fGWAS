\name{fg.load.simple}
\alias{fg.load.simple}
\title{
Loading simple SNP file 
}
\description{
Loading the SNP data from simple text file
}
\usage{
fg.load.simple( file.simple.snp, options=list(verbose=F))
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
SNP, CHR,POS,RefAllel, AltAllel, Sub1,Sub2,Sub3,Sub4,Sub5,...
SNP1,1,1, A, T, 2, 0, 1, 1, 1, ...
SNP2,1,2, A, G, 2, 1, 1, 0, 0, ...
...
}

}

\examples{


snp <- matrix( round( runif(20000, 0, 2) ), nrow=100)
colnames(snp) <- paste("sub", 1:200, sep="");
snp.name <- paste("snp", 1:100, sep="");
snp.mat <- data.frame(snp.name, chr=1, pos=1:100, Allele1="A", Allele2="B", snp);
write.table(snp.mat, file="temp.fgwas.snp.mat.csv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t");
obj.gen <- fg.load.simple("temp.fgwas.snp.mat.csv");


objx <- fg.simulate( "Logistic", "SAD1", 2000, 800, 1:6, phe.missing=0.05, snp.missing=0.05, sig.pos=301, plink.format=FALSE, file.prefix = "test1.fgwas" );
obj.gen <- fg.load.simple("test1.fgwas.geno.tab");
r1 <- fg.snpscan(obj.gen, objx$obj.phe, snp.sub=c(290:310), options=list(ncores=1));


objx <- fg.simulate( "Logistic", "SAD1", 2000, 800, 1:6, phe.missing=0.05, snp.missing=0.05, sig.pos=301, plink.format=TRUE, file.prefix = "test1.fgwas" );
obj.plink <- read.plink( "test1.fgwas.geno.bed", "test1.fgwas.geno.bim", "test1.fgwas.geno.fam"  );
x0 <- as.data.frame( obj.plink$map[,c(2,1,4,5,6)])
x1 <-  matrix( as.numeric( (t(obj.plink$genotypes ))), nrow=800)
x1[which(x1==0)]<-NA
x1 <- x1 - 1;
x1 <-  matrix( x1, nrow=800);
colnames(x1) <- paste("N", 1:2000, sep="_");
snp.mat <- data.frame(x0, x1) ;
write.table(snp.mat, file="temp.fgwas.snp.mat.csv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t");
obj.gen <- fg.load.simple("temp.fgwas.snp.mat.csv");
r2 <- fg.snpscan(obj.gen, objx$obj.phe, snp.sub=c(290:310), options=list(ncores=1));



}
\keyword{ Data }
