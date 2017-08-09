
try( detach("package:fGWAS", unload=T) )
library(fGWAS)

#Full test of GLS method
obj1 <- fg.simulate( "Logistic", "AR1", 2000, 1000, 1:8, phe.missing=0.05, snp.missing=0.05, sig.pos=501, plink.format=FALSE, file.prefix = NULL );
print(obj1$obj.gen);
print(obj1$obj.phe);
plot(obj1$obj.phe, curve.fitting=F);


# GLS method
obj1.gls <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, method="gls", options=list(ncores=20) )
obj1.gls;
plot(obj1.gls);
obj1.gls.sig <- fg.select.sigsnp( obj1.gls, sig.level=0.05, pv.adjust="bonferroni",  options=list() )


#Full test of FAST method
obj1.fast <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, method="fast", curve.type="Exponential", options=list(ncores=1))
obj1.fast;
plot(obj1.fast);
ret <- fg.sig.snp(obj1.fast);


#Full test of fGWAS method
obj1.fgwas <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, method="fgwas", options=list(ncores=20))
save(obj1.gls, obj1.fast, obj1.fgwas, file="temp.fg.simu.obj1.rdata");
obj1.fgwas;
plot(obj1.fgwas);
ret <- fg.sig.snp(obj1.fgwas);


#Full test of fGWAS method
obj1.log <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="AR1", options=list(ncores=20))
obj1$obj.pheY <- matrix( rowMeans(as.matrix(obj1$obj.pheY), na.rm=T), ncol=1)
obj1.scan <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, method="gls", options=list(ncores=1))

#Two variates
obj2 <- fg.simulate( "Logistic", "AR1", 500, 1000, 1:8, phe.missing=0.03, snp.missing=0.03, par.X = c(1,5), sig.pos=500, plink.format=TRUE, file.prefix="temp.fg.simu" );
obj2.gls1 <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, method="gls", options=list(ncores=1))
obj2.fast1 <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, method="fast", options=list(ncores=1))
obj2.gls <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, method="gls", options=list(ncores=20))
obj2.fast <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, method="fast", options=list(ncores=20))
obj2.fast <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, method="fgwas", options=list(ncores=20))

obj2.log <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, curve.type="Logistic", covariance.type="AR1", snp.idx=1:10, options=list(ncores=1))
obj2$obj.pheY <- matrix( rowMeans(as.matrix(obj2$obj.pheY), na.rm=T), ncol=1)
obj2.scan <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, method="gls", options=list(ncores=1))


library(fGWAS)
obj2 <- fg.simulate( "Logistic", "AR1", 500, 1000, 1:8, phe.missing=0.03, snp.missing=0.03, par.X = c(1,5), sig.pos=500, plink.format=TRUE, file.prefix="temp.fg.simu" );
obj2.log <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, curve.type="Logistic", covariance.type="AR1", snp.idx=1:10, options=list(ncores=5))
