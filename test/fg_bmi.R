
try( detach("package:fGWAS", unload=T) )
library(fGWAS)

NCORES<- 10;

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
write.csv(tb[,c("ID", "Y_1", "Y_2", "Y_3", "Y_4", "Y_5", "Y_6", "Y_7", "Y_8" )], file=file.phe, quote=F, row.names=F);
write.csv(tb[,c("ID", "Z_1", "Z_2", "Z_3", "Z_4", "Z_5", "Z_6", "Z_7", "Z_8" )], file=file.phe.time, quote=F, row.names=F);
write.csv(tb[,c("ID", "X_1", "X_2", "X_3", "X_4", "X_5", "X_6" )], file=file.phe.cov, quote=F, row.names=F);

obj.phe <- fg.load.phenotype( file.phe, file.phe.cov, file.phe.time );
obj.phe;

obj.phe <- fg.data.estimate( obj.phe );
obj.phe;

obj.fgwas.100 <- fg.snpscan( obj.gen, obj.phe, method="fgwas", snp.sub=c(1:100), curve.type="Legendre2",covariance.type="AR1", options=list(order=3, ncores=NCORES))
obj.optim.100 <- fg.snpscan( obj.gen, obj.phe, snp.sub=c(1:100), curve.type="Legendre2",covariance.type="AR1", options=list(order=3, ncores=NCORES))

obj.fgwas.100;
obj.optim.100;

tb.optim <- fg.select.sigsnp ( obj.optim.100, sig.level=0.05, pv.adjust="none"  )
plot.fgwas.curve(obj.optim.100, tb.optim$INDEX[1:20], file.pdf="bmi.sig.optim.100.pdf");

obj.fgwas <- fg.snpscan( obj.gen, obj.phe, method="fgwas",  curve.type="Legendre2",covariance.type="AR1", options=list(order=3, ncores=NCORES))
obj.optim <- fg.snpscan( obj.gen, obj.phe, curve.type="Legendre2",covariance.type="AR1", options=list(order=3, ncores=NCORES))
save(obj.fgwas, obj.optim, file="bmi.ret.rdata");

tb.optim <- fg.select.sigsnp ( obj.optim, sig.level=0.05, pv.adjust="none"  )
plot.fgwas.curve(obj.optim, tb.optim$INDEX[1:20], file.pdf="bmi.sig.optim.pdf");

tb.fgwas <- fg.select.sigsnp ( obj.fgwas, sig.level=0.05, pv.adjust="none"  )
plot.fgwas.curve(obj.optim, tb.fgwas$INDEX[1:20], file.pdf="bmi.sig.fgwas.pdf");




