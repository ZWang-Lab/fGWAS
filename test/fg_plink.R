
try( detach("package:fGWAS", unload=T) )
library(fGWAS)

#load genotype data
file.plink.bed = "/home/zw355/proj/gwas2/bmi-c1c2-qc2.bed"
file.plink.bim = "/home/zw355/proj/gwas2/bmi-c1c2-qc2.bim"
file.plink.fam = "/home/zw355/proj/gwas2/bmi-c1c2-qc2.fam"

obj.gen <- fg.load.plink( file.plink.bed, file.plink.bim, file.plink.fam, plink.command=NULL, chr=NULL, options=list())
obj.gen;

#load phenotype data
tb <- read.csv("/home/zw355/proj/gwas2/phe-cov-time-64.csv");
file.phe <- tempfile( fileext = ".csv")
file.phe.cov <- tempfile( fileext = ".csv")
file.phe.time <- tempfile( fileext = ".csv")
tb.y <- tb[,c("ID", "Y_1", "Y_2", "Y_3", "Y_4", "Y_5", "Y_6", "Y_7", "Y_8" )];
tb.y[,-1] <- log(tb.y[,-1])
write.csv(tb.y, file=file.phe, quote=F, row.names=F);
write.csv(tb[,c("ID", "Z_1", "Z_2", "Z_3", "Z_4", "Z_5", "Z_6", "Z_7", "Z_8" )], file=file.phe.time, quote=F, row.names=F);
write.csv(tb[,c("ID", "X_1", "X_2", "X_3", "X_4", "X_5", "X_6" )], file=file.phe.cov, quote=F, row.names=F);

obj.phe <- fg.load.phenotype( file.phe, file.phe.cov, file.phe.time );
obj.phe;
save(obj.phe, obj.gen, file="bmi.obj.phe.rdata");

obj.phe <- fg.data.estimate( obj.phe );
obj.phe;

library(fGWAS)
load("bmi.obj.phe.rdata");
obj.fast.rs770707 <- fg.snpscan( obj.gen, obj.phe, method="fast", snp.sub=c("rs770707"), curve.type="Legendre2", covariance.type="AR1", options=list(order=3, ncores=1))
obj.fast <- fg.snpscan( obj.gen, obj.phe, method="fast", snp.sub=c(1000:3000), covariance.type="AR1", options=list(order=3, ncores=20))
obj.fgwas <- fg.snpscan( obj.gen, obj.phe, snp.sub=c(1000:3000), covariance.type="AR1", options=list(order=3, ncores=20))

library(fGWAS)
load("bmi.obj.phe.rdata");
obj.phe$intercept=T
obj.fgwas <- fg.snpscan( obj.gen, obj.phe, method="fgwas", snp.sub=c(1:100), curve.type="Legendre2",covariance.type="AR1", options=list(order=3, ncores=10))
obj.optim <- fg.snpscan( obj.gen, obj.phe, method="optim-fgwas", snp.sub=c(1:100), curve.type="Legendre2",covariance.type="AR1", options=list(order=3, ncores=10))

tb.snps <- fg.select.sigsnp ( obj.optim, sig.level=0.05, pv.adjust="none"  )
plot.fgwas.curve(obj.optim, tb.snps$INDEX, file.pdf="test.pdf");

library(fGWAS);
load("bmi.obj.phe.rdata");

obj.fgwas2 <- fg.snpscan( obj.gen, obj.phe, method="optim-fgwas", snp.sub=c(1:30, 100:130, 5001:5030), covariance.type="AR1", options=list(order=3, ncores=3))
obj.fgwas <- fg.snpscan( obj.gen, obj.phe, method="fgwas", snp.sub=c(1:30, 100:130, 5001:5030), covariance.type="AR1", options=list(order=3, ncores=3))
obj.fast <- fg.snpscan( obj.gen, obj.phe, method="fast", snp.sub=c(1:30, 100:130, 5001:5030) covariance.type="AR1", options=list(order=3, ncores=3))
obj.mixed <- fg.snpscan( obj.gen, obj.phe, method="mixed", snp.sub=c(1:30, 100:130, 5001:5030), covariance.type="AR1", options=list(order=3, ncores=3))
obj.gls <- fg.snpscan( obj.gen, obj.phe, method="gls", snp.sub=c(1:30, 100:130, 5001:5030), covariance.type="AR1", options=list(order=3, ncores=10))

