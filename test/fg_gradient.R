library(fGWAS)
load("obj1.rdata");

system.time(log1 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="auto", covariance.type="auto", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))
show(log1);

system.time(log1 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="AR1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))
show(log1);

system.time(log2 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="SAD1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))
show(log2);

library(fGWAS)
load("obj1.rdata");
system.time(log3 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="ARMA(1,1)", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log4 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="ARH1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log5 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="CS", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log6 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="CSH", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log7 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="VS", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log8 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="SI", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log9 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="FA1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log10 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="FAH1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log11 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="TOEP", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log12 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="TOEPH", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log13 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="HF", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))


library(fGWAS)
load("obj1.rdata");
system.time(log13 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Bi-Logistic", covariance.type="AR1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))


library(fGWAS)
load("obj1.rdata");
system.time(log14 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="ABRK", covariance.type="AR1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))


library(fGWAS)
load("obj1.rdata");
system.time(log15 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Legendre2", covariance.type="AR1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))


library(fGWAS)
load("obj1.rdata");
system.time(log16 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Legendre3", covariance.type="AR1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log17 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Legendre4", covariance.type="AR1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))


library(fGWAS)
load("obj1.rdata");
system.time(log18 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Power", covariance.type="AR1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log19 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="ChapmanRichard", covariance.type="AR1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))


library(fGWAS)
load("obj1.rdata");
system.time(log20 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Exponential", covariance.type="AR1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log21 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Bi-Exponential", covariance.type="AR1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))

library(fGWAS)
load("obj1.rdata");
system.time(log22 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Pharmacology", covariance.type="AR1", method="fgwas", snp.sub=c(499:502), options=list(ncores=1, use.gradient=T, max.loop=1)))




