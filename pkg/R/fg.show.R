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

print_item <- function(item.name, item.value)
{
	if(!is.null(item.value))
		cat(" ", item.name, ":", item.value, "\n")
	else
		cat(" ", item.name, ":", " ", "\n");
}

print_fgwas_gen_obj <- function(object)
{
	cat("== Genotype Object in fGWAS ==\n");

	print_item( "Plink bed", object$params$file.plink.bed );
	print_item( "Plink bim", object$params$file.plink.bim );
	print_item( "Plink fam", object$params$file.plink.fam );
	print_item( "Data file", object$params$file.simple.snp );
	print_item( "SNP count", object$n.snp );
	print_item( "Total individuals", object$n.ind.total );
	if(!is.null(object$reader))
		object$reader$show();
}


summary_fgwas_gen_obj <- function(object)
{
	return(c(n.snp      = object$n.snp,
			n.ind.total = object$n.ind.total,
			n.ind.used  = object$n.ind.used));
}

#List of 10
# $ options  : list()
# $ params   :List of 5
#  ..$ file.phe.long  : chr "/tmp/RtmpYrkzlf/fileda14176b9c0d.csv"
#  ..$ file.phe.cov   : chr "/tmp/RtmpYrkzlf/fileda147a3872df.csv"
#  ..$ file.phe.time  : chr "/tmp/RtmpYrkzlf/fileda148740849.csv"
#  ..$ curve.type     : chr "auto"
#  ..$ covariance.type: chr "auto"
# $ pheY     : num [1:1678, 1:8] 20.5 28.5 19.5 26.1 28.4 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:1678] "8228" "2294" "7395" "309" ...
#  .. ..$ : chr [1:8] "Y_1" "Y_2" "Y_3" "Y_4" ...
# $ pheT     : int [1:1678, 1:8] 40 44 33 42 47 51 28 29 57 45 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:1678] "8228" "2294" "7395" "309" ...
#  .. ..$ : chr [1:8] "Z_1" "Z_2" "Z_3" "Z_4" ...
# $ pheX     : num [1:1678, 1:6] 2 1 2 2 2 1 2 2 1 1 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:1678] "8228" "2294" "7395" "309" ...
#  .. ..$ : chr [1:6] "X_1" "X_2" "X_3" "X_4" ...
# $ ids      : chr [1:1678] "8228" "2294" "7395" "309" ...
# $ obj.curve:Formal class 'fg.curve.log' [package "fGWAS"] with 2 slots
#  .. ..@ type       : chr "Logistic"
#  .. ..@ description: chr "logistic curve"
# $ obj.covar:Formal class 'fg.covariance.AR1' [package "fGWAS"] with 3 slots
#  .. ..@ par_num    : int(0)
#  .. ..@ type       : chr "AR1"
#  .. ..@ description: chr "AR1"
# $ est.curve:List of 7
#  ..$ type       : chr "Logistic"
#  ..$ intercept  : True or FLASE
#  ..$ param      : Named num [1:3] 36 39.8 29.5
#  .. ..- attr(*, "names")= chr [1:3] "" "" ""
#  ..$ param.lower: num [1:3] 34.1 35.9 26.6
#  ..$ param.upper: num [1:3] 37.8 43.4 32.4
#  ..$ parX       : Named num [1:7] -6.93 -1.59 -1.34 8.28 -5.06 ...
#  .. ..- attr(*, "names")= chr [1:7] "" "X_1" "X_2" "X_3" ...
#  ..$ parX.lower : num [1:7] -8.8828 -1.9554 -9.0214 0.0122 -12.7226 ...
#  ..$ parX.upper : num [1:7] -5.1 -1.29 13.85 15.07 1.86 ...
#  ..$ R2         : num 0.846
# $ est.covar:List of 2
#  ..$ type : chr "AR1"
#  ..$ param: num [1:2] 0.907 19.916
# $ intercept
#
# - attr(*, "class")= chr "fgwas.phe.obj"


print_fgwas_phe_obj <- function(object)
{
	cat("== Phenotype Object in fGWAS ==\n");

	print_item( "Longitudal value", object$params$file.phe.long );
	print_item( " -- Individual count", ifelse( is.null(object$pheY), 0, NROW(object$pheY) ) );
	print_item( "Longitudal time", object$params$file.phe.time );
	print_item( " -- Time count", ifelse( is.null(object$pheT), 0, NCOL(object$pheT) ) );
	print_item( "Covariate file", object$params$file.phe.cov );
	print_item( " -- Covariate count", ifelse( is.null(object$pheX), 0, NCOL(object$pheX) ) );
	print_item( " -- Intercept", ifelse( object$intercept, "YES", "NO") );
	print_item( " -- Estimate values", if(!is.null(object$est.curve)) object$est.curve$parX else "N/A" );

	print_item( "Curve type", object$params$curve.type );
	if(is.null(object$est.curve))
		print_item( " -- Estimation", "Unknown" )
	else
	{
		print_item( " -- Estimate type", object$est.curve$type );
		print_item( " -- Estimate values", object$est.curve$param );
		print_item( " -- R2 ", object$est.curve$R2 );
	}

	print_item( "Covariate type", object$params$covariance.type );
	if(is.null(object$est.covar))
		print_item( " -- Estimation", "Unknown" )
	else
	{
		print_item( " -- Estimated type", object$est.covar$type );
		print_item( " -- Estimate values", object$est.covar$param ) ;
	}
}

summary_fgwas_phe_obj <- function(object)
{
	df <- object$pheY;
	if( !is.null(object$pheX))
		df <- cbind(df, object$pheX[match(rownames(object$pheY), rownames(object$pheX)), ,drop=F]);

	if( !is.null(object$pheT))
		df <- cbind(df, object$pheT[match(rownames(object$pheY), rownames(object$pheT)), ,drop=F]);

	return(df);
}

# 'fgwas.scan.obj':
#    $ params  : list
#  ..$ options : list
#  ..$ system.time : c()
#  ..$ obj.gen : fgwas.gen.obj
#  ..$ obj.phe : fgwas.phe.obj
#  ..$ obj.perm: fgwas.perm.obj
#  ..$ ret.gls : data.frame
#      ..$ result rowname=SNPNAME
#         ..$ INDEX
#         ..$ NAME
#         ..$ CHR
#         ..$ POS
#         ..$ NMISS
#         ..$ MAF
#         ..$ LR
#         ..$ pv
#  ..$ ret.mixed : data.frame
#      ..$ lmer
#      ..$ beta
#      ..$ result rowname=SNPNAME
#         ..$ INDEX
#         ..$ NAME
#         ..$ CHR
#         ..$ POS
#         ..$ MAF
#         ..$ Allel1
#         ..$ Allel2
#         ..$ NMISS
#         ..$ pv
#         ..$ P_min
#         ..$ P_join
#         ..$ P0
#         ..$ P1
#         ..$ P2
#         ..$ P3
#  ..$ ret.fast : data.frame
#  ..$ ret.fgwas : data.frame
#      ..$ result  rowname=SNPNAME
#         ..$ INDEX
#         ..$ NAME
#         ..$ CHR
#         ..$ POS
#         ..$ MAF
#         ..$ NMISS
#         ..$ SNP0
#         ..$ SNP1
#         ..$ SNP2
#         ..$ GENO
#         ..$ LR
#         ..$ pv
#         ..$ h0_X1
#         ..$ h0_X2
#         ..$ h0_cur1
#         ..$ h0_cur2
#         ..$ h0_cur3
#         ..$ h0_cor1
#         ..$ h0_cor2
#         ..$ h1_X1
#         ..$ h1_X2
#         ..$ h1_curA1
#         ..$ h1_curA2
#         ..$ h1_curA3
#         ..$ ..
#         ..$ h1_cor1
#         ..$ h1_cor2
#         ..$ h0_R2
#         ..$ h1_R2


print_fgwas_scan_obj <- function(object)
{
	show_sigsnp <- function( df, method )
	{
		cat("== Result from '", method, "' method ==\n");
		cat(" SNP = ", NROW(df),"\n")
		show(head(df[ order(df[,"pv"],decreasing=F),
			c("INDEX", "NAME","CHR", "POS", "MAF","NMISS", "pv")], n=20));
	}

	ret.item="";
	if(!is.null(object$ret.gls)) {
		ret.item <- "ret.gls";
		show_sigsnp( object$ret.gls$result, "GLS" );
	}

	if(!is.null(object$ret.mixed)) {
		ret.item <- "ret.mixed";
		show_sigsnp( object$ret.mixed$result, "Mixed" );
	}

	if(!is.null(object$ret.fast)) {
		ret.item <- "ret.fast";
		show_sigsnp( object$ret.fast$result, "FAST" );
	}

	if(!is.null(object$ret.fgwas)) {
		ret.item <- "ret.fgwas";
		show_sigsnp( object$ret.fgwas$result, "fGWAS"  );
	}

	cat(paste("\nCheck all SNPs using this variable: 'your_object$", ret.item, "$result'.\n", sep=""));
}

summary_fgwas_scan_obj <- function(object)
{
	ret <- list();
	if(!is.null(object$ret.gls))
		ret <- object$ret.gls$result;

	if(!is.null(object$ret.mixed))
		ret <- object$ret.mixed$result;

	if(!is.null(object$ret.fast))
		ret <- object$ret.fast$result;

	if(!is.null(object$ret.fgwas))
		ret <- object$ret.fgwas$result;

	return(ret);
}

fg_select_sigsnp <- function( obj.scan, sig.level = 0.05, pv.adjust="bonferoni",  options=list())
{
	r <- NULL;
	object <-
	if(!is.null(obj.scan$ret.fgwas)) r <- obj.scan$ret.fgwas$result
	else
	if(!is.null(obj.scan$ret.fast)) r <- obj.scan$ret.fast$result
	else
	if(!is.null(obj.scan$ret.mixed)) r <- obj.scan$ret.mixed$result
	else
	if(!is.null(obj.scan$ret.gls)) r <- obj.scan$ret.gls$result
	else
		return(NULL);

	r <- r[ order(r[,"pv"], decreasing=F), ];
	r$adj.pv <- p.adjust( r$pv, method = pv.adjust );
	r <- r[which(r$adj.pv <= sig.level), ,drop=FALSE]

	return(r);
}

plot_fgwas_scan_obj <- function(x, y=NULL, file.pdf=NULL, ... )
{
	object <- x;

	if( missing(file.pdf) || is.null(file.pdf) )
		file.pdf <- tempfile(pattern="fgwas.plot", tmpdir=getwd(), fileext=".pdf");

	pdf(file.pdf, width=7, height=4)

	if(!is.null(object$ret.gls))
	{
		filter.man <- object$ret.gls$result[, c("CHR", "POS", "pv"), drop=F]
		fpt.plot_manhattan( filter.man , p.05=NA, p.01=0.01/NROW(object$ret.gls$result), map.title="GLS" )
	}

	if(!is.null(object$ret.mixed))
	{
		filter.man <-  object$ret.mixed$result[, c("CHR", "POS", "pv"), drop=F]
		fpt.plot_manhattan( filter.man , p.05=NA, p.01=0.01/NROW(object$ret.mixed$result), map.title="Mixed Model")
	}

	if(!is.null(object$ret.fast))
	{
		filter.man <- object$ret.fast$result[, c("CHR", "POS", "pv"), drop=F]
		fpt.plot_manhattan( filter.man , p.05=NA, p.01=0.01/NROW(object$ret.fast$result), map.title="Fast fGWAS Model")
	}

	if(!is.null(object$ret.fgwas))
	{
		filter.man <- object$ret.fgwas$result[, c("CHR", "POS", "pv"), drop=F]
		fpt.plot_manhattan( filter.man, p.05=NA, p.01=0.01/NROW(object$ret.fgwas$result), map.title="fGWAS Model")
	}

	dev.off();

	cat(" Output manhattan figure(s) into ", file.pdf, "\n")
	invisible(file.pdf);
}

plot_fgwas_phe_obj<-function( obj.phe, file.pdf=NULL, curve.fitting=FALSE, ...)
{
	if( missing(file.pdf) || is.null(file.pdf) )
		file.pdf <- tempfile(pattern="fgwas.plot", fileext=".pdf");

	pdf(file.pdf, width=6, height=5)
	if( curve.fitting && (is.null(obj.phe$est.curve) || is.null(obj.phe$est.covar)))
		obj.phe <- fg_dat_est( obj.phe, obj.phe$obj.curve@type, obj.phe$obj.covar@type );

	par_h0 <- c( obj.phe$est.curve$parX, obj.phe$est.curve$param, obj.phe$est.covar$param);
	fpt.plot_curve ( obj.phe, unlist(par_h0), NULL, NULL, extra=list());

	dev.off();

	cat(" Output curve figure(s) into ", file.pdf, "\n")
}


plot_fgwas_curve<-function( object, snp.sub, file.pdf )
{
	if( missing(file.pdf) || is.null(file.pdf) )
		file.pdf <- tempfile(pattern="fgwas.plot", fileext=".pdf");

	if(is.null(object$obj.phe$intercept))
		object$obj.phe$intercept <- FALSE;

	method <- "";
	if(!is.null(object$ret.fast)) { ret <- object$ret.fast$result; method<-"FAST"; }
	if(!is.null(object$ret.fgwas)) { ret <- object$ret.fgwas$result; method<-"fGWAS"; }

	pdf(file.pdf, width=6, height=5)

	if( is.character( snp.sub) )
		snp.index <- match(snp.sub, ret$NAME)
	else
		snp.index <- match(snp.sub, ret$INDEX);

	snp.index <- snp.index [ !is.na(snp.index) ];
	if(length(snp.index)<=0)
		stop("No SNP is selected for the curve plot.");

	snp.mat  <- object$obj.gen$reader$get_snpmat( ret$INDEX[snp.index] );
	snp.info <- object$obj.gen$reader$get_snpinfo( ret$INDEX[snp.index] );
	ret.set  <- ret[snp.index,, drop=F]

	par_x_num <- ifelse( is.null( object$obj.phe$pheX ), 0 , NCOL( object$obj.phe$pheX ) ) +
	             ifelse (object$obj.phe$intercept, 1, 0) ;
	par_curve_num <- get_param_info( object$obj.phe$obj.curve, object$pheT )$count;
	par_covar_num <- get_param_info( object$obj.phe$obj.covar, object$pheT )$count;

	for(i in 1:length(snp.index))
	{
		par_h0 <- ret.set[i, 14:(13+par_x_num+par_curve_num)]
		par_h1 <- ret.set[i, (14 + par_x_num + par_curve_num + par_covar_num + 1):(14 + par_x_num + par_curve_num + par_covar_num + par_x_num + 3*par_curve_num)]
		fpt.plot_curve ( object$obj.phe, unlist(par_h0), unlist(par_h1), snp.mat$snpmat[,i],
		        extra=list(INDEX=ret.set[i,1],
		        NAME = as.character(ret.set[i,2]),
		        CHR  = ret.set[i,3],
		        POS  = ret.set[i,4],
		        MAF  = ret.set[i,5],
		        NMISS= ret.set[i,6] ,
		        LR2  = ret.set[i,11],
		        METHOD=method))
	}

	dev.off();

	cat(" Output curve figure(s) into ", file.pdf, "\n")

}


profile_fgwas_curve<-function( object, snp.sub )
{
	method <- "";
	if(!is.null(object$ret.fast)) { ret <- object$ret.fast$result; method<-"FAST"; }
	if(!is.null(object$ret.fgwas)) { ret <- object$ret.fgwas$result; method<-"fGWAS"; }

	if( is.character(snp.sub) )
		snp.index <- match(snp.sub, ret$NAME)
	else
		snp.index <- match(snp.sub, ret$INDEX);

	snp.index <- snp.index [ !is.na(snp.index) ];
	if(length(snp.index)<=0)
		stop("No SNP is selected for the curve plot.");

	snp.mat  <- object$obj.gen$reader$get_snpmat( ret$INDEX[snp.index] );
	snp.info <- object$obj.gen$reader$get_snpinfo( ret$INDEX[snp.index] );
	ret.set  <- ret[snp.index,, drop=F]

	par_x_num <- ifelse( is.null( object$obj.phe$pheX ), 1 , 1+NCOL( object$obj.phe$pheX ) );
	par_curve_num <- get_param_info( object$obj.phe$obj.curve, object$pheT )$count;
	par_covar_num <- get_param_info( object$obj.phe$obj.covar, object$pheT )$count;

	options <- list();
	options$max.time=max(object$obj.phe$pheT, na.rm=T);
	options$min.time=min(object$obj.phe$pheT, na.rm=T);

	h0.pv <- c();
	h1.pv <- c();
	for(i in 1:length(snp.index))
	{
		par_h0 <- ret.set[i, 14:(13+par_x_num+par_curve_num+par_covar_num)]
		par_h1 <- ret.set[i, (14 + par_x_num + par_curve_num + par_covar_num + 1):(14 + par_x_num + par_curve_num + par_covar_num + par_x_num + 3*par_curve_num+par_covar_num)]

		r <- fpt.profile ( object$obj.phe, unlist(par_h0), unlist(par_h1), snp.mat$snpmat[,i], options);
		h0.pv <- rbind(h0.pv, r[1,]);
		h1.pv <- rbind(h1.pv, r[2,]);
	}

	return(list(h0=h0.pv, h1=h1.pv, snp=t(snp.mat$snpmat)));
}

fpt.profile <- function( obj.phe, par_h0, par_h1, snp.vec, options)
{
	get_est_vector<-function( pheY0, pheX0, pheT0, par_X, par_curve, par_covar )
	{
		X.cov <- matrix( pheX0 %*% par_X, ncol=1, byrow=T) %*% matrix(rep(1,NCOL(pheT0)), nrow=1, byrow=T)
		Y.delt <- pheY0 - get_curve( obj.phe$obj.curve, par_curve, pheT0, options=options )  - X.cov;
		mat.cov <- get_matrix( obj.phe$obj.covar, par_covar, pheT0 );
		pv <- try( dmvnorm_fast( Y.delt, rep(0,NROW(mat.cov)), mat.cov, NULL, log=T), .RR("try.silent") );

		return(pv);
	}

	pheY <- obj.phe$pheY;
	pheX <- obj.phe$pheX;
	pheT <- obj.phe$pheT;

	if ( obj.phe$intercept )
		if (is.null(pheX)) pheX <- matrix( 1, nrow=NROW(pheY), ncol=1 ) else pheX <- cbind(1, pheX);

	par_x_num <- ifelse( is.null( pheX ), 0 , NCOL( pheX ) );
	par_curve_num <- get_param_info(obj.phe$obj.curve, pheT )$count;

	par_x <- c();
	if(par_x_num>0)	par_x <- par_h0[1:par_x_num]
	h0.pv <- get_est_vector( pheY, pheX, pheT, par_x, par_h0[(par_x_num+1):(par_x_num+par_curve_num)], par_h0[-c(1:(par_x_num+par_curve_num))] );

	snp0.idx <- which(snp.vec==0);
	snp1.idx <- which(snp.vec==1);
	snp2.idx <- which(snp.vec==2);

	par_x <- c();
	if(par_x_num>0)	par_x <- par_h1[1:par_x_num]

	h1.pv <- rep(NA, length(h0.pv))
	if(length(snp0.idx)>0)
		h1.pv[snp0.idx] <- get_est_vector( pheY[snp0.idx,,drop=F], pheX[snp0.idx,,drop=F], pheT[snp0.idx,,drop=F],
		       par_x, par_h1[(par_x_num+1):(par_x_num+par_curve_num)], par_h1[-c(1:(par_x_num+par_curve_num*3))] );

	if(length(snp1.idx)>0)
		h1.pv[snp1.idx] <- get_est_vector( pheY[snp1.idx,,drop=F], pheX[snp1.idx,,drop=F], pheT[snp1.idx,,drop=F],
		       par_x, par_h1[(par_x_num+1+par_curve_num):(par_x_num+par_curve_num*2)], par_h1[-c(1:(par_x_num+par_curve_num*3))]);

	if(length(snp2.idx)>0)
		h1.pv[snp2.idx] <- get_est_vector( pheY[snp2.idx,,drop=F], pheX[snp2.idx,,drop=F], pheT[snp2.idx,,drop=F],
		      par_x, par_h1[(par_x_num+1+par_curve_num*2):(par_x_num+par_curve_num*3)], par_h1[-c(1:(par_x_num+par_curve_num*3))]);

	return(rbind(h0.pv, h1.pv));
}
