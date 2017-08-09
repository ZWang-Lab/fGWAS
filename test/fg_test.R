try( detach("package:fGWAS", unload=T) )
library(fGWAS)

NCORES <- 20;

#No covariate
obj1 <- fg.simulate( "Logistic", "AR1", 2000, 1000, 1:8, phe.missing=0.05, snp.missing=0.05, sig.pos=501, plink.format=FALSE, file.prefix = NULL );
obj1 <- fg.simulate( "Logistic", "AR1", 2000, 1000, 1:8, phe.missing=0.05, snp.missing=0.05, sig.pos=501, plink.format=FALSE, file.prefix = NULL, par.X=c(2, 3.2) );
obj1.gls   <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, method="gls", options=list(ncores=20))
obj1.mixed <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, method="mixed",options=list(ncores=20))
obj1.fast  <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, method="fast", curve.type="auto", covariance.type="AR1", options=list(ncores=NCORES))
obj1.fgwas <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, method="fgwas",curve.type="auto", covariance.type="auto", options=list(ncores=NCORES))
obj1.opt   <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, method="optim-fgwas", curve.type="Logistic", covariance.type="auto", options=list(ncores=NCORES))
obj1.opt2  <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, curve.type="Logistic", covariance.type="AR1", options=list(ncores=20))

obj1$obj.pheY <- matrix( rowMeans(as.matrix(obj1$obj.pheY), na.rm=T), ncol=1)
obj1.gls2 <- fg.snpscan( obj1$obj.gen, obj1$obj.phe, method="gls", options=list(ncores=1))

#Two variates
obj2 <- fg.simulate( "Logistic", "AR1", 500, 1000, 1:8, phe.missing=0.03, snp.missing=0.03, par.X = c(1,5), sig.pos=500, plink.format=TRUE, file.prefix="gwas_test" );
obj2.gls   <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, method="gls",  options=list(ncores=NCORES))
obj2.mixed <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, method="mixed",options=list(ncores=NCORES))
obj2.fast  <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, method="fast",  curve.type="auto", options=list(ncores=NCORES))
obj2.fgwas <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, method="fgwas", curve.type="auto", options=list(ncores=NCORES))
obj2.fast  <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, curve.type="auto", options=list(ncores=NCORES))
obj2.opt   <- fg.snpscan( obj2$obj.gen, obj2$obj.phe, method="optim-fgwas", curve.type="Logistic", covariance.type="AR1", snp.sub=c(400:600), options=list(ncores=NCORES))

obj2$obj.pheY <- matrix( rowMeans(as.matrix(obj2$obj.pheY), na.rm=T), ncol=1)
obj2.scan <- fg.snpscan( obj2$obj.gen, obj2$obj.phe,  method="gls", options=list(ncores=1))

