# Phenotype File( CSV Format, Header=T  Comma seprated )
#
# Format 1:
# shareid, PHE1, PHE2, ...
#
# Covariate File( CSV Format, Header=T  Comma seprated  )
# Format 2:
#
# shareid, COV1,....
#
# Format 3:
# shareid, TIME1, TIME2, ...

check_pheno_file<-function( file.phe.long, file.phe.cov, file.phe.time)
{
	err.info <- "Error!";
	cat("Checking phenotype file......\n");
	cat("* PHE.LONG =", file.phe.long , "\n");

	if(!file.exists(file.phe.long))
		return(list(error=T, err.info=paste("Longitudinal file is not found, file=",file.phe.long, sep="") ) );

	phe.long <- try( read.csv(file.phe.long, header=T, stringsAsFactors=F, row.names=1) );
	if (class(phe.long)=="try-error")
	{
		cat("! Can not open file(", file.phe.long, ")\n");
		return(list(error=T, err.info=err.info));
	}
	else{
		cat("* Individuals =", NROW(phe.long), "\n");
		cat("* Times =", NCOL(phe.long), "\n");
		cat("* Mean =",  colMeans(phe.long, na.rm=T), "\n");
		cat("* SD =",    colSds(phe.long, na.rm=T),"\n");

		cat("  First 5 Items:\n");
		show(head(phe.long, n=5));
	}

	cat("* PHE.TIME =", file.phe.time , "\n");
	phe.time <- NULL;
	if(!is.null(file.phe.time)  && !is.na(file.phe.time) )
	{
		if(!file.exists(file.phe.time))
			return(list(error=T, err.info=paste("Time file is not found, file=",file.phe.time, sep="") ) );

		phe.time <- try( read.csv(file.phe.time, header=T, stringsAsFactors=F, row.names=1) );
		if (class(phe.time)=="try-error")
		{
			cat("! Can not open file(", file.phe.time, ")\n");
			return(list(error=T, err.info=err.info));
		}
		else
		{
			cat("* Individuals =", NROW(phe.time), "\n");
			cat("* Times =", NCOL(phe.time), "\n");
			cat("* Mean =",  colMeans(phe.time, na.rm=T), "\n");
			cat("* SD =",    colSds(phe.time, na.rm=T),"\n");
			cat("  First 5 Items:\n");
			show(head(phe.time, n=5));
		}

		idx.inter <- intersect( rownames(phe.long), rownames(phe.time) );
		if( !( length(idx.inter)==NROW(phe.long) && length(idx.inter)==NROW(phe.time) ) )
		{
			cat("! PHE.LONG don't have consistent IDs with PHE.TIME.\n" );
			return(list(error=T, err.info=err.info));
		}
	}

	cat("Checking covariate file......\n");

	cat("* COV.FILE =", file.phe.cov , "\n");
	if(!is.null(file.phe.cov) && !is.na(file.phe.cov) )
	{
		if(!file.exists(file.phe.cov))
			return(list(error=T, err.info=paste("Covariate file is not found, file=",file.phe.cov, sep="") ) );

		phe.cov <- try( read.csv(file.phe.cov, header=T, stringsAsFactors=F, row.names=1) );
		if (class(phe.cov)=="try-error")
		{
			cat("! Can not open file(", file.phe.cov, ")\n");
			return(list(error=T, err.info=err.info));
		}
		else{
			cat("* Individuals =", NROW(phe.cov), "\n");
			cat("* Covariate =", NCOL(phe.cov), "\n");
			cat("* Mean =",  colMeans(phe.cov, na.rm=T), "\n");
			cat("* SD =",    colSds(phe.cov, na.rm=T), "\n");
			cat("  First 5 Items:\n");
			show(head(phe.cov, n=5));
		}

		y.ncov <- NCOL(phe.cov);
		all.na <- which( is.na( rowSums(phe.cov, na.rm=T)	) )
		if (length(all.na)>0 )
		{
			cat("!", length(all.na), "IDs dont have non-missing data.\n" );
			cat("! First 5 IDs:\n");
			show( head( phe.cov[all.na, ], n=5) );
		}
	}


	all.na <- which( is.na( rowSums(phe.long, na.rm=T)	) )
	if (length(all.na)>0 )
	{
		cat("!", length(all.na), "IDs dont have non-missing data.\n" );
		cat("! First 5 IDs:\n");
		show( head( phe.long[all.na, ], n=5) );
	}

	return(list(error=F))
}

show_options<-function(options)
{
##@ TODO

}

fg_load_phenotype <-function( file.phe.long, file.phe.cov=NULL, file.phe.time=NULL, curve.type="auto", covariance.type="auto", file.plot.pdf=NULL, intercept=T, options=list())
{
	cat( "[ Loading Phenotype Data]\n");
	cat( "Checking the parameters ......\n");

	if ( missing(file.phe.long) )
		stop("! file.phe.long must be assigned with the valid file name.");

	if(missing(curve.type) || is.null(curve.type) || is.na(curve.type) ) curve.type <- "auto";
	if(missing(covariance.type) || is.null(covariance.type) || is.na(covariance.type) ) covariance.type <- "auto";
	if(missing(file.phe.cov) || is.null(file.phe.cov) || is.na(file.phe.cov) )  file.phe.cov <- NULL;
	if(missing(file.phe.time) || is.null(file.phe.time) || is.na(file.phe.time) ) file.phe.time <- NULL;
	if(missing(intercept) || is.null(intercept) || is.na(intercept) ) intercept <- TRUE;

	cat("* Phenotypic Data File = ",  file.phe.long, "\n");
	cat("* Covariate Data File = ",  file.phe.cov, "\n");
	cat("* Time Data File = ",  file.phe.time, "\n");
	cat("* Curve Type = ",  curve.type, "\n")
	cat("* Covariance Type = ",   covariance.type, "\n")
	cat("* Intercept = ",   intercept, "\n")

	default_options <- list();

	if (missing(options))
		options <- default_options
	else
	{
		options0 <- default_options;
		options0[names(options)] <- options;
		options <- options0;
	}

	cat( "Checking the optional items......\n");
	show_options( options);

	ret <- list();
	ret$options <- options;
	ret$params <- list( 	file.phe.long  = file.phe.long,
				file.phe.cov   = file.phe.cov,
				file.phe.time  = file.phe.time,
				curve.type     = curve.type,
				covariance.type= covariance.type,
				intercept      = intercept);

	r.chk <- check_pheno_file( file.phe.long, file.phe.cov, file.phe.time );
	if(r.chk$error)
		stop( r.chk$err.info );

	ret$pheY <- read.csv(file.phe.long, header=T, stringsAsFactors=F, row.names=1 );

	if (!is.null(file.phe.time) && !is.na(file.phe.time))
		ret$pheT <- read.csv(file.phe.time, header=T, stringsAsFactors=F, row.names=1)
	else
	{
		ret$pheT <- matrix( rep(c(1:NCOL(ret$pheY)), NROW(ret$pheY)), nrow = NROW(ret$pheY))
		rownames(ret$pheT) <- rownames(ret$pheY);
	}

	if (!is.null(file.phe.cov) && !is.na(file.phe.cov))
		ret$pheX <- read.csv(file.phe.cov,  header=T, stringsAsFactors=F, row.names=1)
	else
		ret$pheX <- NULL;

	##check the id consistent
	ret$ids <- rownames(ret$pheY);

	if(!is.null(ret$pheX) && is.data.frame(ret$pheX)) ret$pheX <- as.matrix(ret$pheX)
	if(!is.null(ret$pheY) && is.data.frame(ret$pheY)) ret$pheY <- as.matrix(ret$pheY)
	if(!is.null(ret$pheT) && is.data.frame(ret$pheT)) ret$pheT <- as.matrix(ret$pheT)

	r.est  <- fg_dat_est( ret, curve.type,  covariance.type, file.plot.pdf );
	if( r.est$error )
		stop(r.est$err.info)
	else
	{
		ret$obj.curve      <- r.est$obj.curve;
		ret$obj.covar      <- r.est$obj.covar;
		ret$est.curve      <- r.est$est.curve;
		ret$est.covar      <- r.est$est.covar;
		ret$summary.curve  <- r.est$summary.curve;
		ret$summary.covariance  <- r.est$summary.covar;
	}

	class(ret) <- "fgwas.phe.obj";
	return(ret);
}


## for loading data and compareing id i pheX and pheY, pheT
keep_gls.fgwas <- function( snp.mat, pheY, pheZ, pheX=NULL, op.cpu=0)
{
	sample.ids <- intersect( rownames(phe.mat), colnames(snp.mat)[-c(1:2)] );
	if(length(sample.ids)==0)
	{
		cat("! No same population data in the phenotypic data and genotypic data.\n");
		return(list(error=T, err.info="No same population data in the phenotypic data and genotypic data."))
	}

	phe.mat <- phe.mat[ sample.ids, , drop=F]
	snp.mat <- cbind( snp.mat[,c(1,2)], snp.mat[ , sample.ids, drop=F]);

	y.p <- grep( Y.prefix, colnames(phe.mat) );
	z.p <- grep( Z.prefix, colnames(phe.mat) );
	x.p <- c();
	if(!is.null(covar.names))
		x.p <- match( covar.names, colnames(phe.mat) );

	Y.sub <- substring( colnames(phe.mat)[y.p], nchar(Y.prefix)+1);
	Z.names <- paste( Z.prefix, Y.sub, sep="");
	z.p <- match( Z.names, colnames(phe.mat) );
	if( length( which(is.na(z.p)) )>0 )
	{
		cat("! The Z columns are not matched with Y columns.\n");
		return(list(error=T, err.info="The Z columns are not matched with Y columns."))
	}

	len.x <- length(x.p);
	len.y <- length(y.p);
	len.z <- length(z.p);
	if ( len.y != len.z)
	{
		cat("! The Z columns are not matched with Y columns.\n");
		return(list(error=T, err.info="The Z columns are not matched with Y columns."))
	}

	if( length(x.p)>0 && length(which(is.na(x.p)))>0)
	{
		cat("! The covariate names are not matched with phenotypical data.\n");
		return(list(error=T, err.info="The covariate names are not matched with phenotypical data."))
	}

	phe.mat <- cbind(phe.mat, ID = c(1:NROW(phe.mat)) );
	id.p <- NCOL(phe.mat);

	phe.gls.mat <- array(0,dim=c(0, 1 + len.x + 2 ))
	for(i in 1:len.y)
	{
		phe.tmp <- phe.mat[, c( ID=id.p, Y=y.p[i], Z=z.p[i], x.p), drop=F ];

		phe.tmp.col <- colnames(phe.tmp);
		phe.tmp.col[ c(1:3) ] <- c("ID", "Y", "Z");
		colnames(phe.tmp) <- phe.tmp.col;

		phe.tmp.missing <- which( is.na(phe.tmp$Z) | is.na(phe.tmp$Y) );
		if(length(phe.tmp.missing)>0)
			phe.tmp <- phe.tmp[-phe.tmp.missing, , drop=F ];
		phe.gls.mat <- rbind( phe.gls.mat, phe.tmp );
	}

	reg.str0 <- ifelse( len.x>0, paste("Y ~ ", paste(covar.names,collapse= "+")),  "Y ~ 1" );
	reg.str1 <- ifelse( len.x>0, paste("Y ~ ", paste(covar.names,collapse= "+"), "+ as.factor(SNP)") ,  "Y ~ 1 + as.factor(SNP)" );

	cat("* H0 =", as.character(reg.str0), "\n" );
	cat("* H1 =", as.character(reg.str1), "\n" );

	r0 <- try( gls( as.formula(reg.str0), phe.gls.mat, correlation = corAR1(form = ~ Z | ID ), method="ML" ) );
	if(class(r0)=="try-error")
	{
		cat("! Failed to call gls() method.\n");
		return(list( error=T, err.info="Failed to call gls() method.") )
	}

	cpu.fun<-function( sect )
	{
		if( (sect-1)*n.percpu+1 > NROW(snp.mat) )
			return(NULL);

		range.fr <- (sect-1)*n.percpu+1;
		range.to <- sect*n.percpu;
		range.to <- ifelse( range.to > NROW(snp.mat),  NROW(snp.mat), range.to );

		r.gls <- array( NA, dim = c( (range.to-range.fr+1), 7 ) );
		r.gls[,1]<-c(range.fr:range.to);

		reg.f0 <- as.formula(reg.str0);
		reg.f1 <- as.formula(reg.str1);

		for(i in c(range.fr:range.to) )
		{
			snp <- unlist( snp.mat[i, c(3:NCOL(snp.mat)) ] );
			snp.gls <- snp[ phe.gls.mat$ID ];
			gls.mat <- cbind( phe.gls.mat, SNP=snp.gls );

			maf <- mean(snp.gls, na.rm=T)/2;
			if(maf>0.5) maf <- 1-maf;

			pv.max <- NA;
			if( maf <= 10^(-4)) pv.max <- 1.0;

			# remove NA row
			na.row <- unique( which(is.na(gls.mat))%% NROW(gls.mat) );
			if(length(na.row)>0) gls.mat <- gls.mat[-na.row,,drop=F];

			r0 <- try( do.call("gls", args = list(reg.f0, gls.mat, correlation = corAR1(form = ~ Z | ID ), method="ML") ) );
			r1 <- try( do.call("gls", args = list(reg.f1, gls.mat, correlation = corAR1(form = ~ Z | ID ), method="ML") ) );
			if(any(class(r1)=="try-error") || any(class(r0)=="try-error") )
			{
				r.gls[i-range.fr+1,] <- c( i, snp.mat[i,1], snp.mat[i,2], length(na.row), maf, NA, pv.max );
				next;
			}

			r <- try( anova( r0,r1 ) );
			if(any(class(r)=="try-error"))
			{
				r.gls[i-range.fr+1,] <- c( i, snp.mat[i,1], snp.mat[i,2], length(na.row), maf,NA, pv.max );
				next;
			}

			r.gls[(i-range.fr+1),] <- c( i, snp.mat[i,1], snp.mat[i,2], length(na.row), maf, r[2,8], r[2,9]);
		}

		return(r.gls);
	}


	r.gls <- c();
	if( op.cpu>1 && require("snowfall") )
	{
		cat("Starting parallel computing, snowfall/snow......\n");
		snowfall::sfInit(parallel = TRUE, cpus = op.cpu, type = "SOCK")

		n.percpu <- ceiling( NROW(snp.mat) / op.cpu );
		snowfall::sfExport("n.percpu", "phe.gls.mat", "snp.mat", "covar.names", "reg.str0", "reg.str1" );

		gls.cluster <- snowfall::sfClusterApplyLB( 1:op.cpu, cpu.fun);
		snowfall::sfStop();

		cat("Stopping parallel computing......\n");

		r.gls <- do.call( rbind, gls.cluster );
	}
	else
	{
		cat("Starting the fGWAS estimate for each SNP......\n");
		n.percpu <- NROW(snp.mat);
		r.gls <- cpu.fun(1);
	}

	colnames(r.gls) <- c("SNP.IDX", "CHR", "POS", "N.miss", "MAF", "L.Ratio", "pv");
	rownames(r.gls) <- rownames(snp.mat)

	return(list(error=F, r=r.gls));
}
