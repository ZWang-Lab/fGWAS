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
	ret$intercept <- intercept;
	ret$options <- options;
	ret$params <- list( 	file.phe.long  = file.phe.long,
				file.phe.cov   = file.phe.cov,
				file.phe.time  = file.phe.time,
				curve.type     = curve.type,
				covariance.type= covariance.type);

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

