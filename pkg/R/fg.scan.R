check_ids<-function( obj.gen, obj.phe )
{
	stopifnot( all(rownames(obj.phe$pheT)==rownames(obj.phe$pheX)) && all(rownames(obj.phe$pheX)==rownames(obj.phe$pheY)) )

	ids.phe <- rownames(obj.phe$pheY);
	ids.gen <- obj.gen$reader$get_individuals();

	m.phe <- match( as.character(ids.phe), as.character(ids.gen)  );
	if (length(which(is.na(m.phe))) > 0 )
	{
		cat("!", length(which(is.na(m.phe))), "IDs, can not be found in the genotype file.\n" );
		cat("! First 5 IDs:\n");
		show( head(ids.gen[is.na(m.phe), ], n=5 ) );
	}

	m.snp <- match( as.character(ids.gen) , as.character(ids.phe) );
	if (length(which(is.na(m.snp))) > 0 )
	{
		cat("!", length(which(is.na(m.snp))), "IDs, can not be found in the phenotype file.\n" );
		cat("! First 5 IDs:\n");
		show( head(ids.phe[is.na(m.snp), ], n=5 ) );
	}

	ids.com <- intersect( ids.phe, ids.gen);
	return(ids.com);
}

fg_snpscan <-function( obj.gen, obj.phe, method, curve.type=NULL, covariance.type=NULL, snp.sub = NULL, permutation=NULL, options=list())
{
	if( class(obj.gen)!="fgwas.gen.obj" )
		stop("The first paramater should be genotype object.");

	if( class(obj.phe)!="fgwas.phe.obj" )
		stop("The second paramater should be phenotype object.");

	default_options <- list( max.optim.failure=20, min.optim.success=2, intercept=F, degree=3, ncores=1, verbose=FALSE, use.gradient=F, piecewise=1000, use.snowfall=TRUE);
	default_options[names(options)] <- options;
	options <- default_options;
	options$opt.method <- toupper(method);

	if(options$verbose)
		cat( "[ SNP Scanning ] \n");

	ids.com <- check_ids(obj.gen, obj.phe);
	obj.phe <- select_individuals( obj.phe, ids.com );
	obj.gen$reader$select_individuals( ids.com );

	if(toupper(method )!="GLS" && toupper(method )!="MIXED")
	{
		if(is.null(curve.type) && inherits( obj.phe$obj.curve, "fg.curve.base") )
			curve.type <- obj.phe$obj.curve@type;
		if(is.null(covariance.type) && inherits( obj.phe$obj.covar, "fg.covariance.base") )
			covariance.type <- obj.phe$obj.covar@type;
		if( toupper(curve.type )!= toupper(obj.phe$obj.curve@type) || toupper(covariance.type) != toupper(obj.phe$obj.covar@type ) )
			obj.phe <- fg_dat_est( obj.phe, curve.type, covariance.type, options=list(verbose=options$verbose) );
	}

	if(options$verbose)
		cat( "Genetic Effect Analysis by '", method,"' method......\n", sep="");

	ret <- list( error=FALSE );
	if(is.null(snp.sub)) snp.sub <- 1:obj.gen$n.snp;

	if(is.character(snp.sub))
	{
		snp.sub <- obj.gen$reader$get_snpindex( snp.sub );
		snp.sub <- snp.sub  [ which(!is.na(snp.sub)) ];
	}

	if (length(snp.sub)==0)
		stop("No SNPs are found in the genotype file.");

	r.time <- NULL;
	if( toupper(method )=="GLS")
	{
		r.time <- system.time( r <- snpscan_nocurve( obj.gen, obj.phe, snp.sub, options ) );
		if (r$error)
			stop(r$err.info);

		ret$ret.gls <- r;
	}

	if( toupper(method )=="MIXED")
	{
		r.time <- system.time( r <- fg_mixed_scan( obj.gen, obj.phe, snp.sub, degree=options$degree, ncores=options$ncores ) );
		if (r$error)
			stop(r$err.info);

		ret$ret.mixed <- r;
	}

	if( toupper(method )=="FAST" || toupper(method )=="FAST-NORM" || toupper(method )=="FGWAS" || toupper(method )=="OPTIM-FGWAS")
	{
		if(is.null(obj.phe$est.curve))
			obj.phe <- fg_dat_est( obj.phe,
					curve.type = ifelse( is.null(obj.phe$obj.curve), "auto", obj.phe$obj.curve@type),
					covariance.type = ifelse( is.null(obj.phe$obj.curve), "auto", obj.phe$obj.covar@type), options );
		obj.phe$h0 <- proc_est_h0( obj.phe, options );
	}

	if( toupper(method )=="FAST" || toupper(method )=="FAST-NORM")
	{
		r.time <- system.time( r <- snpsnp_curve( obj.gen, obj.phe, snp.sub, options ) );
		if (r$error)
			stop(r$err.info);

		ret$ret.fast <- r;
	}

	if( toupper(method )=="FGWAS" || toupper(method )=="OPTIM-FGWAS")
	{
		r.time <- system.time( r <- snpsnp_curve( obj.gen, obj.phe, snp.sub, options ) );
		if (r$error)
			stop(r$err.info);

		ret$ret.fgwas <- r;
	}

	ret$obj.gen <- obj.gen;
	ret$obj.phe <- obj.phe;
	ret$system.time <- r.time;
	ret$options <- options;

	class(ret) <- "fgwas.scan.obj";

	## dont load plink data, just scanning result object
	return( ret );
}


snpsnp_curve<-function(obj.gen, obj.phe, snp.idx=NULL, options )
{
	if(is.null(snp.idx)) snp.idx <- 1:obj.gen$n.snp
	snp.len <- length(snp.idx);
	snp.sect0 <- seq( 1, snp.len, options$piecewise );
	snp.sect1 <- c( snp.sect0[-1]-1, snp.len );

	intercept <- ifelse(is.null(obj.phe$intercept), FALSE, obj.phe$intercept );
	r.list <- list();

	for(i in 1:length(snp.sect0))
	{
		cpu.fun<-function( k )
		{
			## dont use requireNamespace method which doesnt call .onAttach() hook.
			##requireNamespace("fGWAS");
			require(fGWAS);

			i.start <- ifelse( (k-1)*snps.percore+1 > NROW(snp.idx.sub), NROW(snp.idx.sub), (k-1)*snps.percore+1);
			i.stop  <- ifelse( k*snps.percore>NROW(snp.idx.sub), NROW(snp.idx.sub), k*snps.percore );
			r.mle <- c()
			for(i in i.start:i.stop)
			{
				r <- try( proc_mle( snp.idx.sub[i], snp.info[i,], snp.mat$snpmat[,i], obj.phe, intercept, options ), silent =FALSE );
				if(class(r)!="try-error")
					r.mle <- rbind( r.mle, r);
			}
			return( r.mle );
		}

		cat("  Calculated SNP Range =", snp.idx[snp.sect0[i]], snp.idx[snp.sect1[i]], "\n" );

		snp.idx.sub <- snp.idx[ snp.sect0[i]:snp.sect1[i] ];
		snp.mat  <- obj.gen$reader$get_snpmat( snp.idx.sub);
		snp.info <- obj.gen$reader$get_snpinfo( snp.idx.sub);

		snps.percore <- ceiling( NROW(snp.idx.sub)/options$ncores );
		used.ncores <-  ceiling( NROW(snp.idx.sub)/snps.percore );

		r.cluster <- list();
		if( options$ncores>1 && options$use.snowfall==TRUE && require(snowfall) )
		{
			if(options$verbose) cat("Starting parallel computing(snowfall/snow)......\n");

			sfInit(parallel = TRUE, cpus = options$ncores, type = "SOCK")

			sfExport("snp.idx.sub", "snp.mat", "snp.info", "obj.phe", "options", "snps.percore" );

			r.cluster <- sfClusterApplyLB( 1:used.ncores, cpu.fun );

			sfStop();

			if(options$verbose) cat("Stopping parallel computing(snowfall/snow)\n");
		}
		else
		{
			if(options$verbose)	cat("Starting piecewise analysis(parallel package)......\n");

			r.cluster <- mclapply( 1:used.ncores, cpu.fun, mc.cores=options$ncores );

			if(options$verbose)	cat("Stopping piecewise.\n");
		}

		r.list[[i]] <- do.call( "rbind", r.cluster );
	}

	r.fgwas <- do.call( "rbind", r.list );

	curve.par.names <- get_param_info( obj.phe$obj.curve, obj.phe$pheT )$names;
	covar.par.names <- get_param_info( obj.phe$obj.covar, obj.phe$pheT )$names;
	par.curve.len   <- get_param_info( obj.phe$obj.curve, obj.phe$pheT )$count;
	par.covar.len   <- get_param_info( obj.phe$obj.covar, obj.phe$pheT )$count;

	re.names <- colnames(r.fgwas);
	re.names[1:4] <- c("INDEX", "NAME", "CHR", "POS");
	re.names[5:13] <- c("MAF", "NMISS", "SNP0", "SNP1", "SNP2", "GENO", "LR", "pv", "L0");

	if(intercept)
		par.X.len <- ifelse( is.null(obj.phe$pheX), 1, 1+NCOL(obj.phe$pheX) )
	else
		par.X.len <- ifelse( is.null(obj.phe$pheX), 0, NCOL(obj.phe$pheX) );

	idx.start <- 14;
	if(par.X.len>0)
	{
		re.names[idx.start:(idx.start+par.X.len-1)] <- paste("h0.X", c(1:par.X.len)-ifelse(intercept,1,0), sep="")
		idx.start <- idx.start + par.X.len
	}

	re.names[idx.start:(idx.start+par.curve.len-1)] <- paste("h0", curve.par.names, sep=".")
	idx.start <- idx.start+par.curve.len
	re.names[idx.start:(idx.start + par.covar.len-1)] <- paste("h0X", covar.par.names, sep=".")
	idx.start <- idx.start + par.covar.len
	re.names[idx.start] <- "L1";
	idx.start <- idx.start + 1
	if(par.X.len>0)
	{
		re.names[idx.start:(idx.start + par.X.len-1)] <- paste("h1.X", c(1:par.X.len)-ifelse(intercept,1,0), sep="")
		idx.start <- idx.start + par.X.len
	}

	re.names[idx.start:(idx.start+par.curve.len-1)] <- paste("h1.G0", curve.par.names, sep=".")
	idx.start <- idx.start+par.curve.len
	re.names[idx.start:(idx.start+par.curve.len-1)] <- paste("h1.G1", curve.par.names, sep=".")
	idx.start <- idx.start+par.curve.len
	re.names[idx.start:(idx.start+par.curve.len-1)] <- paste("h1.G2", curve.par.names, sep=".")
	idx.start <- idx.start+par.curve.len
	re.names[idx.start:(idx.start+par.covar.len-1)] <- paste("h1X", covar.par.names, sep=".")

	re.names[idx.start+par.covar.len] <- "h0.R2"
	re.names[idx.start+par.covar.len+1] <- "h1.R2"

	colnames(r.fgwas ) <- re.names;
	rownames(r.fgwas ) <- r.fgwas$NAME;

	return(list(error=F, result=r.fgwas));
}

snp_cluster <- function(obj.gen, snp.idx=NULL, dist=3 )
{
	get_min_dist<-function( snp0 )
	{
		min_dist_a0a1a2<-function(a0,a1,a2)
		{
			snpx <- snp0;
			snpx[which(snp0==0)]<-a0
			snpx[which(snp0==1)]<-a1;
			snpx[which(snp0==2)]<-a2;

			dists <- unlist( lapply(1:length(r.cluster), function(grp){
				sum( snpx - r.cluster[[grp]]$snp0 != 0, na.rm=T ) ;
			}));

			return(list(min.grp=which.min(dists), dist=min(dists, na.rm=T) ) );
		}

		d1 <- min_dist_a0a1a2(0, 1, 2)
		d2 <- min_dist_a0a1a2(0, 2, 1)
		d3 <- min_dist_a0a1a2(1, 0, 2)
		d4 <- min_dist_a0a1a2(1, 2, 0)
		d5 <- min_dist_a0a1a2(2, 0, 1)
		d6 <- min_dist_a0a1a2(2, 1, 0)

		min.grp=c(d1$min.grp, d2$min.grp, d3$min.grp, d4$min.grp, d5$min.grp, d6$min.grp);
		min.dist=c(d1$dist,    d2$dist,    d3$dist,    d4$dist,    d5$dist,    d6$dist);

		idx <- which.min( min.dist );

		return(list( min.grp=min.grp[idx], dist=min.dist[idx] ) );
	}

	if(is.null(snp.idx)) snp.idx <- 1:obj.gen$n.snp
	snp.len <- length(snp.idx);
	snp.sect0 <- seq( 1, snp.len, 1000 );
	snp.sect1 <- c( snp.sect0[-1]-1, snp.len );

	snp.idx.sub <- snp.idx[ snp.sect0[1]:snp.sect1[1] ];
	snp.mat  <- obj.gen$reader$get_snpmat( snp.idx.sub);

	r.cluster <- list();
	r.cluster[[1]] <- list(snp0=snp.mat$snpmat[,1], grp=c(snp.idx[ snp.sect0[1] ]) )

	for(i in 1:length(snp.sect0))
	{
		#if(options$verbose)
		cat("  Calculated SNP Range =", snp.idx[snp.sect0[i]], snp.idx[snp.sect1[i]], "\n" );

		snp.idx.sub <- snp.idx[ snp.sect0[i]:snp.sect1[i] ];
		snp.mat  <- obj.gen$reader$get_snpmat( snp.idx.sub);

		for(k in 1:NROW(snp.idx.sub))
		{
			min.dist <- get_min_dist( snp.mat$snpmat[,k] )
			cat(k, min.dist$dist, min.dist$min.grp, "\n" );
			if(min.dist$dist<=dist)
				r.cluster[[min.dist$min.grp]]$grp <- c(r.cluster[[min.dist$min.grp]]$grp, snp.idx.sub[k])
			else
			{
				r.cluster[[length(r.cluster) +1 ]] <- list(snp0=snp.mat$snpmat[,k], grp=c(snp.idx.sub[k]) )
			}
		}
	}

	return(r.cluster);
}