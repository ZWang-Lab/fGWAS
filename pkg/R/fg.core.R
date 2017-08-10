#------------------------------------------------------------------------------------------
#
# Functional Mapping Method
#
# !!!! fy[which(fy<=.Machine$double.eps)] <- .Machine$double.eps;
#
#------------------------------------------------------------------------------------------
proc_est_h0<-function( obj.phe, options )
{
	pheY <- obj.phe$pheY;
	pheT <- obj.phe$pheT;
	pheX <- obj.phe$pheX;

	if(obj.phe$intercept)
	{
		if ( !is.null(obj.phe$pheX) )
			pheX <- cbind(1, obj.phe$pheX )
		else
			pheX <- array(1, dim=c(NROW(obj.phe$pheY),1))
	}

	options$max.time <- max(pheT, na.rm=T);
	options$min.time <- min(pheT, na.rm=T);
	if(is.null(options$min.optim.success)) options$min.optim.success <- 5;
	options$b.permu  <- FALSE;

	h0 <- proc_full_h0( obj.phe, pheY, pheX, pheT, options, h0.ref=NULL );

	n.par.curve <- get_param_info(obj.phe$obj.curve, pheT, options)$count;

	parin.X     <- NULL;
	parin.curve <- h0$par;
	if (!is.null(pheX))
	{
		parin.X     <- h0$par[1:NCOL(pheX)];
		parin.curve <- h0$par[-c(1:NCOL(pheX))];
	}

	parin.covar <- parin.curve[-c(1:n.par.curve)];
	parin.curve <- parin.curve[c(1:n.par.curve)];

	X <- 0;
	if (!is.null(pheX))
		X  <- matrix( rep( pheX %*% parin.X , NCOL(pheY) ), byrow=F, ncol=NCOL(pheY) );

	mu <- get_curve( obj.phe$obj.curve, parin.curve, pheT, options=options);
    Y.delt <- pheY - mu - X;
	mat.cov <- get_matrix( obj.phe$obj.covar, parin.covar, pheT );

	h0$pv  <- try( dmvnorm_fast( Y.delt, rep(0,NROW(mat.cov)), mat.cov, nna.vec=NULL, log=T), .RR("try.silent") );
    h0$ssr <- sum(Y.delt^2, na.rm=T)

	return(h0);
}

proc_mle<-function( index, snp.info, snp.vec, obj.phe, intercept=F, options=list() )
{
	default_options = list(b.permu=F);
	default_options[names(options)]= options;
	options = default_options;

	gen.par <- max(snp.vec, na.rm=T) - min(snp.vec, na.rm=T) + 1;

	pheY <- obj.phe$pheY;
	pheT <- obj.phe$pheT;
	pheX <- obj.phe$pheX;

	if(intercept)
	{
		if ( !is.null(obj.phe$pheX) )
			pheX <- cbind(1, obj.phe$pheX )
		else
			pheX <- array(1, dim=c(NROW(obj.phe$pheY),1))
	}

	h0.ref  <- sum(-1 * obj.phe$h0$pv, na.rm=T);
	snp.miss <- which(is.na(snp.vec) | snp.vec<0 | snp.vec>2 );
	if (length(snp.miss)>0)
	{
		if (!is.null(pheX))
			pheX <- pheX[-snp.miss,,drop=F];
		pheY <- pheY[-snp.miss,,drop=F];
		pheT <- pheT[-snp.miss,,drop=F];
		snp.vec <- snp.vec[-snp.miss];
		h0.ref  <- sum( -1*obj.phe$h0$pv[-snp.miss], na.rm=T);
	}

	options$max.time=max(pheT, na.rm=T);
	options$min.time=min(pheT, na.rm=T);
	options$max.optim.failure <- 20;

	if(options$opt.method=="FAST" || options$opt.method=="FAST-NORM" )
	{
		if(is.null(options$min.optim.success)) options$min.optim.success <- 1;
		h0 <- proc_fast_h0( obj.phe, pheY, pheX, pheT, options );

		if(is.null(options$min.optim.success)) options$min.optim.success <- 2;
		h1 <- proc_fast_h1( obj.phe, pheY, pheX, pheT, snp.vec, options, h0);
	}
	else if(options$opt.method=="FGWAS")
	{
		if(is.null(options$min.optim.success)) options$min.optim.success <- 1;
		h0 <- proc_full_h0( obj.phe, pheY, pheX, pheT, options, h0.ref );
		if(is.null(options$min.optim.success)) options$min.optim.success <- 2;
		h1 <- proc_full_h1( obj.phe, pheY, pheX, pheT, snp.vec, options, h0 );
	}
	## for default method = 'OPTIM-FGWAS'
	else
	{
		if(is.null(options$min.optim.success)) options$min.optim.success <- 1;
		h0 <- proc_fast_h0( obj.phe, pheY, pheX, pheT, options );
		if(is.null(options$min.optim.success)) options$min.optim.success <- 2;
		h1 <- proc_fast_h1( obj.phe, pheY, pheX, pheT, snp.vec, options, h0  );

		r.val <- ( h0$value - h1$value )*2;
		r.pv <- pchisq( r.val, df=(gen.par-1)*get_param_info( obj.phe$obj.curve, pheT )$count, lower.tail=FALSE );
		if( r.pv <= 0.001 )
		{
			if(is.null(options$min.optim.success)) options$min.optim.success <- 1;
				h0 <- proc_full_h0( obj.phe, pheY, pheX, pheT, options, h0.ref );
			if(is.null(options$min.optim.success)) options$min.optim.success <- 2;
				h1 <- proc_full_h1( obj.phe, pheY, pheX, pheT, snp.vec, options, h0  );
		}
	}

	bOptim <- TRUE;
	if (is.null(h0) || is.na(h0$value) || is.infinite(h0$value) )
		bOptim <- FALSE;
	if (is.null(h1) || is.na(h1$value) || is.infinite(h1$value) )
		bOptim <- FALSE;

	if(bOptim)
	{
		r.val <- ( h0$value - h1$value )*2
		r.pv <- pchisq( r.val, df=(gen.par-1)*get_param_info( obj.phe$obj.curve, pheT )$count, lower.tail=FALSE );
	}
	else
		r.val <- r.pv <- NA;


	maf <- ifelse( mean(snp.vec)/2>0.5, 1-mean(snp.vec)/2, mean(snp.vec)/2 );

	re <- data.frame( index, snp.info[c(1:3)], matrix(c(maf, length(snp.miss), sum(snp.vec==0), sum(snp.vec==1), sum(snp.vec==2),
			gen.par, r.val, pv=r.pv, h0$value, h0$par, h1$value, h1$par, h0$R2, h1$R2 ) , nrow=1 ), stringsAsFactors=F);

	if( .RR("debug")  )
		cat("snp[", as.character(re[1,2]), "] loci=", re[1,3], re[1,4], "LR2=", re[1,11], re[1,12], "Param=", unlist( re[1, -c(1:12)] ), "\n");

	return(re);
}

proc_full_h0 <- function( obj.phe, pheY, pheX, pheT, options, h0.ref=NULL )
{
	gradient_h0_cov<-function( parin, pheY, pheT, pheX, obj.curve, obj.covariance, nna.vec, n.par.curve, options )
	{
		if (!is.null(pheX))
		{
			parin.X     <- parin[1:NCOL(pheX)];
			parin0      <- parin[-c(1:NCOL(pheX))];
			parin.curve <- parin0[1:n.par.curve];
			parin.covar <- parin0[-c(1:n.par.curve)];
		}
		else
		{
			parin.X     <- NULL;
			parin.curve <- parin[1:n.par.curve];
			parin.covar <- parin[-c(1:n.par.curve)];
		}

		d <- rep(0, length(parin));
		mat.cov <- get_matrix( obj.covariance, parin.covar, pheT );

		for(nna in unique(nna.vec))
		{
			nna.idx <- which(nna.vec == nna);
			col.idx <- !is.na(pheY[nna.idx[1],]);

			mat.cov.nna <- mat.cov[col.idx, col.idx, drop=F];
			mat.inverse <-  solve(mat.cov.nna);

			X <- rep(0, NROW(nna.idx));
			if (!is.null(pheX))
				X  <- matrix( rep(  pheX[nna.idx,] %*% parin.X , sum(col.idx) ), byrow=F, ncol=sum(col.idx) );

			mu <- get_curve( obj.curve, parin.curve, pheT[nna.idx,col.idx,drop=F], options=options);
		    Y.delt <- pheY[nna.idx,col.idx,drop=F]- mu - X;
		    d.u <- mat.inverse %*%t(Y.delt);

			d.X <- c ();
			if (!is.null(pheX))
			{
				for(i in 1:NCOL(pheX) )
					d.X <- c(d.X, sum( -1 * t(d.u) * pheX[nna.idx,i], na.rm=T ) );
			}

			d.curve <- c();
			d.list <- get_gradient( obj.curve, parin.curve , pheT[nna.idx,col.idx,drop=F], options=options);
			for(i in 1:length(parin.curve))
				d.curve <- c(d.curve, sum( -1 * d.u * t(d.list[[i]]), na.rm=T) );

			Y.delt[is.na(Y.delt)] <- 0;
			d.list <- get_gradient( obj.covariance, parin.covar, pheT[nna.idx,col.idx,drop=F] );
			d.sigma <- -0.5*(mat.inverse*NROW(Y.delt) - mat.inverse%*% ( t(Y.delt) %*% Y.delt) %*% mat.inverse);
			d.cov <- c();
			for(i in 1:length(parin.covar))
				d.cov <- c(d.cov, sum(d.sigma * d.list[[i]]) );

			d <- d + c(d.X, d.curve, -1*d.cov);
		}

		return(d);
	}

	proc_h0_cov<-function( parin, pheY, pheT, pheX, obj.curve, obj.covariance, nna.vec, n.par.curve, options )
	{
		if (!is.null(pheX))
		{
			parin.X     <- parin[1:NCOL(pheX)];
			parin       <- parin[-c(1:NCOL(pheX))];
			parin.curve <- parin[1:n.par.curve];
			parin.covar <- parin[-c(1:n.par.curve)];
		}
		else
		{
			parin.X     <- NULL;
			parin.curve <- parin[1:n.par.curve];
			parin.covar <- parin[-c(1:n.par.curve)];
		}

		X <- rep(0, NROW(pheY));
		if (!is.null(pheX))
			X  <- matrix( rep(  pheX %*% parin.X , NCOL(pheY) ), byrow=F, ncol=NCOL(pheY) );

		mu <- get_curve( obj.curve, parin.curve, pheT, options=options);
	    Y.delt <- pheY - mu - X;

		mat.cov <- get_matrix( obj.covariance, parin.covar, pheT );

		pv <- try( dmvnorm_fast( Y.delt, rep(0,NROW(mat.cov)), mat.cov, nna.vec, log=T), .RR("try.silent") );
		if ( class(pv)=="try-error" || any(is.na(pv)) )
			return (NaN);

		A <- sum(pv);
		return(-A);
	}

	nna.vec <- get_non_na_number(pheY);
	n.par.curve <- get_param_info(obj.phe$obj.curve, pheT, options)$count;

	## TEST gradient function
	if ( !is.null(options$use.gradient) &&  options$use.gradient && !.RR("graident.test", FALSE) && requireNamespace("numDeriv", quietly=T))
	{
		parin <- c(obj.phe$est.curve$parX, obj.phe$est.curve$param, obj.phe$est.covar$param)
		r0 <- gradient_h0_cov( parin, pheY, pheT, pheX, obj.phe$obj.curve, obj.phe$obj.covar, nna.vec, n.par.curve, options)
		r1 <- numDeriv::grad( function(u) proc_h0_cov( u, pheY, pheT, pheX, obj.phe$obj.curve, obj.phe$obj.covar, nna.vec, n.par.curve, options), parin);
		if ( max( abs(r0-r1) ) > 0.1 )
		{
			cat("Gradient=", r0, "\n");
			cat("numDeriv=", r1, "\n");
			stop("The gradient function doesn't work well.\n");
		}

		.RW("graident.test", TRUE);
	}

	h0.failure <- 0;
	h0.best <- NULL;
	while( h0.failure < 5 )
	{

		system.time( h0.best <- optim_BFGS( obj.phe$est.curve$parX, obj.phe$est.curve$param, obj.phe$est.covar$param,
					proc_h0_cov, proc_gr = if(options$use.gradient) gradient_h0_cov else NULL,
					pheY = pheY, pheT = pheT, pheX = pheX,
					obj.curve = obj.phe$obj.curve, obj.covariance = obj.phe$obj.covar,
					nna.vec = nna.vec, n.par.curve = n.par.curve, options = options));

		if( !is.null(h0.ref) && h0.best$value > h0.ref  * 1.01 )
		{
			h0.failure <- h0.failure + 1;
			options$min.optim.success <- options$min.optim.success * 2;
			h0.best <- NULL;
			reset_seed();
		}
		else
			break;
	}

	if( is.null(h0.best) || is.infinite(h0.best$value) || is.na(h0.best$value))
		h0.best <- list(par=c(obj.phe$est.curve$parX, obj.phe$est.curve$param, obj.phe$est.covar$param), value=NA, R2=NA)
	else
		h0.best$R2 <- get_R2 (pheY, pheX, pheT, obj.phe$obj.curve, h0.best$par );

	if( .RR("debug") && !options$b.permu )
		cat( "H0[F]", h0.best$value, "Param=",h0.best$par, "\n");

	return(h0.best);
}

proc_full_h1<-function( obj.phe, pheY, pheX, pheT, snp.vec, options, h0)
{
	proc_only_mu<-function( parin, pheY, pheT, pheX, snp.vec, obj.curve, sum.type )
	{
		if(!is.null(pheX))
			mu  <- get_curve( obj.curve, parin[-c(1:NCOL(pheX))], pheT, options=options )
		else
			mu  <- get_curve( obj.curve, parin, pheT, options=options )

		if(!all(is.na(mu) == is.na(pheT)) )
			return(NaN);

		X <- matrix( rep(  pheX %*% parin[1:NCOL(pheX)] , NCOL(pheY) ), byrow=F, ncol=NCOL(pheY) );
		Y.delt <- pheY - mu - X;

	    if(sum.type==0)
		    A <- sum(Y.delt^2, na.rm=T)
		else
		    A <- Y.delt[!is.na(Y.delt)];

		return(A);
	}

	gradient_h1_cov<-function( parin, pheY, pheT, pheX, obj.curve, obj.covar, nna.vec, snp.vec, n.par.curve, n.gentype=n.gentype, options  )
	{
		if (!is.null(pheX))
		{
			parin.X <- parin[1:NCOL(pheX)];
			parin0  <- parin[-c(1:NCOL(pheX))];
			parin.curve <- parin0[1:(n.par.curve * n.gentype)];
			parin.covar <- parin0[-c(1:(n.par.curve*n.gentype))];
		}
		else
		{
			parin.X     <- NULL;
			parin.curve <- parin[1:(n.par.curve * n.gentype)];
			parin.covar <- parin[-c(1:(n.par.curve*n.gentype))];
		}

		d <- rep(0, length(parin));
		mat.cov <- get_matrix( obj.covar, parin.covar, pheT );

		idx.k <- c();
		for(nna in unique(nna.vec))
		{
			nna.idx <- which(nna.vec == nna);
			col.idx <- !is.na(pheY[nna.idx[1],]);

			mat.cov.nna <- mat.cov[col.idx, col.idx, drop=F];
			mat.inverse <- solve(mat.cov.nna);

			if(!is.null(pheX))
				Y.delt <- pheY[nna.idx, col.idx, drop=F]  -  matrix( rep(  pheX[nna.idx,,drop=F] %*% parin.X , sum(col.idx) ), byrow=F, ncol=sum(col.idx) )
			else
				Y.delt <- pheY[nna.idx, col.idx, drop=F];

			snp.vec.nna <- snp.vec[nna.idx];
			parin.curve0 <- parin.curve;
			for(k in 0:2)
			{
				idx.k <- which(snp.vec==k & nna.vec == nna)
				if (length(idx.k ) >0 )
				{
					parin.curveX <- parin.curve0[1:n.par.curve];
					parin.curve0 <- parin.curve0[-c(1:n.par.curve)];

					Y.delt[snp.vec.nna==k, ] <- Y.delt[snp.vec.nna==k, , drop=F] - get_curve(obj.curve, parin.curveX, pheT[idx.k, col.idx, drop=F], options=options );
				}
			}

		    d.u <- mat.inverse %*% t(Y.delt);

			d.X <- c ();
			if (!is.null(pheX))
				for( i in 1:NCOL(pheX))
					d.X <-  c(d.X, sum( -1 * t(d.u) * pheX[nna.idx,i], na.rm=T) );

			d.curve <- c();
			parin.curve0 <- parin.curve;
			for(k in 0:2)
			{
				if( sum(snp.vec==k) > 0 )
				{
					parin.curveX <- parin.curve0[1:n.par.curve];
					parin.curve0 <- parin.curve0[-c(1:n.par.curve)];

					idx.k <- which(snp.vec==k & nna.vec == nna)
					if (length(idx.k ) >0 )
					{
						d.list <- get_gradient( obj.curve, parin.curveX, pheT[idx.k, col.idx, drop=F], options=options);
						for(i in 1:n.par.curve)
							d.curve <- c(d.curve, sum( -1 * t(d.u[, snp.vec.nna==k, drop=F]) * d.list[[i]], na.rm=T) )
					}
					else
						d.curve <- c(d.curve, rep(0, n.par.curve) )
				}
			}

			Y.delt[is.na(Y.delt)] <- 0;
			d.list <- get_gradient( obj.covar, parin.covar, pheT[nna.idx, col.idx, drop=F] );
			d.sigma <- -0.5*(mat.inverse*NROW(Y.delt) - mat.inverse%*% ( t(Y.delt) %*% Y.delt) %*% mat.inverse);
			d.cov <- c();
			for(i in 1:length(parin.covar))
				d.cov <- c(d.cov, sum(d.sigma * d.list[[i]]) );

			d <- d + c(d.X, d.curve, -1*d.cov);
		}


		return(d);
	}

	proc_h1_cov<-function( parin, pheY, pheT, pheX, obj.curve, obj.covar, nna.vec, snp.vec, n.par.curve, n.gentype=n.gentype, options )
	{
		if (!is.null(pheX))
		{
			parin.X <- parin[1:NCOL(pheX)];
			parin <- parin[-c(1:NCOL(pheX))];
			parin.curve <- parin[1:(n.par.curve * n.gentype)];
			parin.covar <- parin[-c(1:(n.par.curve*n.gentype))];
		}
		else
		{
			parin.curve <- parin[1:(n.par.curve * n.gentype)];
			parin.covar <- parin[-c(1:(n.par.curve*n.gentype))];
		}

		if(!is.null(pheX))
			Y.delt <- pheY -  matrix( rep(  pheX %*% parin.X , NCOL(pheY) ), byrow=F, ncol=NCOL(pheY) )
		else
			Y.delt <- pheY;

		idx.k <- c();
		for(k in 0:2)
		{
			idx.k <- which(snp.vec==k)
			if (length(idx.k ) >0 )
			{
				parin.curve0 <- parin.curve[1:n.par.curve];
				parin.curve  <- parin.curve[-c(1:n.par.curve)];

				Y.delt[idx.k,] <- Y.delt[idx.k,, drop=F] - get_curve(obj.curve, parin.curve0, pheT[idx.k,, drop=F], options=options );
			}
		}

		mat.cov <- get_matrix( obj.covar, parin.covar, pheT );
		pv <- try(dmvnorm_fast( Y.delt, rep(0, NCOL(mat.cov)), mat.cov, nna.vec=nna.vec, log=T), .RR("try.silent") );
		if ( class(pv)=="try-error" || is.na(pv) )
			return (NaN);

		return( - sum( pv ) );
	}

	nna.vec <- get_non_na_number( pheY );
	n.gentype <- length(unique(snp.vec));
	n.par.curve <- get_param_info(obj.phe$obj.curve, pheT, options)$count;
	parin.x     <- if( is.null(pheX)) NULL else h0$par[ 1:NCOL(pheX) ];
	X.ncol <- if( is.null(pheX)) 0 else NCOL(pheX);
	parin.curve <- rep( h0$par[ (X.ncol+1):(X.ncol+n.par.curve) ], n.gentype );
	parin.covar <- h0$par[ -(1:(X.ncol+n.par.curve)) ] ;

	if(0)
	{
	if ( !is.null(options$use.gradient) &&  options$use.gradient && !.RR("graident.test", FALSE) && requireNamespace("numDeriv", quietly=T) )
	{
		parin <- c(parin.x, parin.curve*1.0, parin.covar);
		r0 <- gradient_h1_cov( parin, pheY, pheT, pheX, obj.phe$obj.curve, obj.phe$obj.covar, nna.vec, snp.vec, n.par.curve, n.gentype, options)
		r1 <- numDeriv::grad( function(u) proc_h1_cov( u, pheY, pheT, pheX, obj.phe$obj.curve, obj.phe$obj.covar, nna.vec, snp.vec, n.par.curve, n.gentype, options), parin);
		if ( max( abs(r0-r1) ) > 0.1 )
		{
			cat("Gradient=", r0, "\n");
			cat("numDeriv=", r1, "\n");
			stop("The gradient function doesn't work well.\n");
		}

		.RW("graident.test", TRUE);
	}
	}


	if(options$use.gradient)
		h1.best <-  optim_BFGS ( parin.x, parin.curve*1.0, parin.covar,
					  proc_h1_cov, proc_gr= gradient_h1_cov,
					  pheY = pheY, pheT = pheT, pheX = pheX,
					  obj.curve = obj.phe$obj.curve, obj.covar = obj.phe$obj.covar,
					  nna.vec=nna.vec, snp.vec=snp.vec,
					  n.par.curve=n.par.curve, n.gentype=n.gentype, options=options )
	else
		h1.best <-  optim_BFGS ( parin.x, parin.curve*1.0, parin.covar,
					  proc_h1_cov, proc_gr= NULL,
					  pheY = pheY, pheT = pheT, pheX = pheX,
					  obj.curve = obj.phe$obj.curve, obj.covar = obj.phe$obj.covar,
					  nna.vec=nna.vec, snp.vec=snp.vec,
					  n.par.curve=n.par.curve, n.gentype=n.gentype, options=options );

	if( is.null(h1.best) || is.infinite(h1.best$value) || is.na(h1.best$value))
	{
		parin.curve <- rep( h0$par[ (X.ncol+1):(X.ncol+n.par.curve) ], 3);
		h1.best <- list(par=c(parin.x, parin.curve, parin.covar), value=NA, R2=NA)
	}
	else
	{
		if(!is.null(pheX))
		{
			parin.X <- h1.best$par[1:NCOL(pheX)];
			parin <- h1.best$par[-c(1:NCOL(pheX))];
		}
		else
		{
			parin.X <- NULL;
			parin <- h1.best$par;
		}

		idx.k <- c();
		parin.curve <- parin[1:(n.par.curve * n.gentype)];
		parin.covar <- parin[-c(1:(n.par.curve * n.gentype))];
		par.fin <- c();
		for(k in 0:2)
		{
			idx.k <- which(snp.vec==k)
			if ( length(idx.k) > 0 )
			{
				par.fin <- rbind(par.fin, parin.curve[1:n.par.curve]);
				parin.curve  <- parin.curve[-c(1:n.par.curve)];
			}
			else
				par.fin <- rbind(par.fin, rep(NA,n.par.curve)) ;
		}

		h1.best$par <- c(parin.X, t(par.fin), parin.covar);
		h1.best$R2 <- get_R2 (pheY, pheX, pheT, obj.phe$obj.curve, h1.best$par, snp.vec );
	}

	if( .RR("debug")  && !options$b.permu )
		cat( "H1[F]", h1.best$value, "Param=", h1.best$par, "\n");

	return(h1.best);

}

optim_least_square<-function ( parin, fn_mle, pheY, pheT, pheX, snp.vec=NULL, parX=NULL, obj.curve=NULL, method=NULL, ssr.ref=NULL )
{
	f_optim <- function()
	{
		parinx<-  parin* runif( length(parin), 0.9, 1.1 );
		control <- list(maxit = 500 + length(parin) * 200, reltol=1e-8);
		loop.optim <- 1;

		while(TRUE)
		{
			h0 <- try( optim( parinx, fn_mle, pheY = pheY, pheT = pheT, pheX = pheX, snp.vec = snp.vec, parX=parX,
							obj.curve = obj.curve, sum.type=0,
							method = ifelse(loop.optim%%2==0, "Nelder-Mead", "BFGS" ),
							control= control ),
							.RR("try.silent") );

			if (class(h0)=="try-error" || any(is.na(h0)) || h0$convergence!=0 )
			{
				if ( is.list(h0) )
				{
					if( h0$convergence == 1)
					{
						control$maxit <- control$maxit*2;
						if ( control$maxit > 500*4096 )
						{
							control$maxit <- 500*4096;
							control$reltol <- control$reltol*2;
						}
					}
					else
						cat("optim() convergence=", h0$convergence, " ", h0$message, "\n")
				}
				reset_seed();
				loop.optim <- loop.optim + 1;
				next;
			}
			else
				break;
		}
		return(h0);
	}

	f_minpack.lm <- function()
	{
		parinx<-  parin* runif( length(parin), 0.9, 1.1 );
		control <- list(maxiter = 500 + length(parin) * 100 );
		if ( control$maxiter > 1024 ) control$maxiter <- 1024;

		temp.ret <- list();

		while(TRUE)
		{
			h0.nls <-  try( nls.lm(par=parinx, fn = fn_mle,
						pheY = pheY, pheT = pheT, pheX = pheX, snp.vec = snp.vec, parX = parX,
						obj.curve = obj.curve, 	sum.type=1,
						control = nls.lm.control(nprint=0, maxiter=control$maxiter ) ),
						.RR("try.silent") );

			if (class(h0.nls)=="try-error" || any(is.na(h0.nls)) || !(h0.nls$info %in% c(1,2,3) ) )
			{
				if ( is.list(h0.nls) )
				{
					if( h0.nls$info == 9)
					{
						control$maxiter <- control$maxiter*2;
						if ( control$maxiter > 1024 ) control$maxiter <- 1024;
					}
					else
						cat("nls.lm(), info=", h0.nls$info, "", h0.nls$message, "\n")
				}
				reset_seed();
				next;
			}
			else
			{
				reset_seed();
				h0.nls$value <- sum(h0.nls$fvec^2);
				temp.ret[[ length(temp.ret) + 1 ]] <- h0.nls;

				if(is.null(ssr.ref))
					break
				else
				if ( h0.nls$value < ssr.ref * 1.1 )
				{
					index <- which.min(unlist(lapply(temp.ret, function(x){x$value})))
					h0.nls <- temp.ret[[index]];
					break;
				}
			}
		}
		return(h0.nls);
	}


	## !!! remove "fast-norm" method, use it by default.
	##
	#if(method != "FAST-NORM")
		A1 <- f_minpack.lm()
	#else
	#	A1 <- f_optim();

	return( A1 );

}

proc_fast_h0 <- function( obj.phe, pheY, pheX, pheT, options)
{
	proc_h0_curve<-function( parin, pheY, pheT, pheX, snp.vec, parX, obj.curve, sum.type )
	{
		if(!is.null(pheX))
		{
			parin.X <- parin[1:NCOL(pheX)];
			parin   <- parin[-c(1:NCOL(pheX))];
			X  <- matrix( rep(  pheX %*% parin.X , NCOL(pheY) ), byrow=F, ncol=NCOL(pheY) );
		}
		else
			X <- rep(0, NROW(pheY));

		mu <- get_curve( obj.curve, parin, pheT, options=options);
	    Y.delt <- pheY - mu - X;

	    if(sum.type==0)
		{
			Y.var <- colVars(Y.delt);
		    Y.delt[is.na(Y.delt)] <- 0;
		    pv <- try( dmvnorm( Y.delt, sigma=diag(Y.var), log=T) );
		    A <- -sum(pv);
		}
		else
		{
			A <- Y.delt[!is.na(Y.delt)];
		}

		return(A);
	}

	proc_h0_covmatrix<-function( parin, Y.delt, pheT, obj.covariance, nna.vec, options )
	{
		mat.cov <- get_matrix( obj.covariance, parin, pheT );

		pv <- try( dmvnorm_fast( Y.delt, rep(0,NROW(mat.cov)), mat.cov, nna.vec=nna.vec, log=T), .RR("try.silent") );
		if ( class(pv)=="try-error" || any(is.na(pv)) )
			return (NaN);

		A <- sum(pv);
		return(-A);
	}

	parinx <- c( obj.phe$est.curve$parX, obj.phe$est.curve$param );
	h0.curve <- optim_least_square( parinx, proc_h0_curve, pheY = pheY , pheT = pheT, pheX = pheX, snp.vec= NULL, parX=NULL, obj.curve = obj.phe$obj.curve, method=(options$opt.method=="FAST")  );
	if( is.null(h0.curve) || is.na(h0.curve$value) || is.infinite(h0.curve$value) )
		return(list(par=c(parinx, obj.phe$est.covar$param), value=NA, R2=NA));

	X <- 0;
	parin.curve <- h0.curve$par;
	if (!is.null(pheX))
	{
		parin.X <- h0.curve$par[1:NCOL(pheX)];
		parin.curve <- h0.curve$par[-c(1:NCOL(pheX))];
		X  <- matrix( rep(  pheX %*% parin.X , NCOL(pheY) ), byrow=F, ncol=NCOL(pheY) );
	}

	nna.vec = get_non_na_number( pheY );
	mu <- get_curve( obj.phe$obj.curve, parin.curve, pheT, options=options);
    Y.delt <- pheY - mu - X;
	h0.cm <- optim_BFGS( c(), c(), obj.phe$est.covar$param,
					proc_h0_covmatrix, Y.delt=Y.delt, pheT=pheT,
					obj.covariance = obj.phe$obj.covar,
					nna.vec = nna.vec,
					options=options);

	if( is.null(h0.cm) || is.na(h0.cm$value) || is.infinite(h0.cm$value) )
		h0.cm <- list(par=c( h0.curve$par, obj.phe$est.covar$param), value=NA, R2=NA)
	else
	{
		mat.cov <- get_matrix( obj.phe$obj.covar, h0.cm$par, pheT );
		h0.cm$pv  <- try( dmvnorm_fast( Y.delt, rep(0,NROW(mat.cov)), mat.cov, nna.vec=nna.vec, log=T), .RR("try.silent") );
	    h0.cm$par <- c( h0.curve$par, h0.cm$par);
	    h0.cm$ssr <- sum(Y.delt^2, na.rm=T)
		h0.cm$R2 <- get_R2 (pheY, pheX, pheT, obj.phe$obj.curve, h0.cm$par);
	}

	if( .RR("debug") && !options$b.permu )
		cat( "H0[F]", h0.cm$value, h0.cm$ssr, "Param=",h0.cm$par, "\n");

	return(h0.cm);
}

proc_fast_h1_comb<-function( obj.phe, pheY, pheX, pheT, snp.vec, options, h0, snp.comb=c(0,1,2), parX=NULL )
{
	proc_h1_curve<-function( parin, pheY, pheT, pheX, snp.vec, parX, obj.curve, sum.type)
	{
		if(!is.null(pheX) && is.null(parX) )
		{
			parin.X <- parin[1:NCOL(pheX)];
			parin.curve <- parin[-c(1:NCOL(pheX))];
		}
		else
		{
			parin.curve <- parin;
			parin.X <- parX;
		}

		n.par.curve	<- get_param_info( obj.curve, pheT, options)$count;

		Y.delt <- c();
		for(k in 0:2)
		{
			idx.k <- which(snp.vec==k)
			if (length(idx.k ) >0 )
			{
				if(!is.null(pheX))
				{
					pheX.k<- try( pheX[idx.k, , drop = F] )
					if(class(pheX.k)=="try-error") browser();
					X  <- matrix( rep(  pheX.k %*% parin.X , NCOL(pheY) ), byrow=F, ncol=NCOL(pheY) );
				}
				else
					X <- rep(0, NROW(idx.k));

				parin.curve0 <- parin.curve[1:n.par.curve];
				Y.delt <- rbind(Y.delt, pheY[idx.k,, drop=F] - get_curve( obj.curve, parin.curve0, pheT[idx.k,, drop=F], options=options ) - X );
			}

			parin.curve  <- parin.curve[- c(1:n.par.curve)];
		}

	    if(sum.type==0)
		{
			Y.var <- colVars(Y.delt);
		    Y.delt[is.na(Y.delt)] <- 0;
		    pv <- try( dmvnorm( Y.delt, sigma=diag(Y.var), log=T) );
		    A <- -sum(pv);

		}
		else
		{
			A <- Y.delt[!is.na(Y.delt)];
		}

		return(A);
	}

	proc_h1_cov<-function( parin, Y.delt, pheT, obj.covariance, nna.vec, options )
	{
		mat.cov <- get_matrix( obj.covariance, parin, pheT );

		pv <- try(dmvnorm_fast( Y.delt, rep(0, NCOL(mat.cov)), mat.cov, nna.vec, log=T), .RR("try.silent") );
		if ( class(pv)=="try-error" || is.na(pv) )
			return (NaN);

		return( - sum(pv) );
	}

	if( is.null(parX))
		parin    <- c( obj.phe$est.curve$parX, rep( obj.phe$est.curve$param, 3 ) )
	else
		parin    <- c( rep( obj.phe$est.curve$param, 3 ) );

	snp.vec <- snp.comb[snp.vec+1];
	h1.curve <- optim_least_square( parin, proc_h1_curve, pheY = pheY, pheT = pheT, pheX = pheX, snp.vec=snp.vec, parX=parX, obj.curve = obj.phe$obj.curve, method=(options$opt.method=="FAST"), ssr.ref=h0$ssr );
	if( is.null(h1.curve) || is.na(h1.curve$value) || is.infinite(h1.curve$value) )
		return(list(par=c(parin, obj.phe$est.covar$param), value=NA, R2=NA));

	Y.delt <- matrix( NA, nrow=NROW(pheY), ncol=NCOL(pheY) );
	if(is.null(pheX) || !is.null(parX) )
	{
		parin.curve <- matrix( h1.curve$par, nrow=3, byrow=T);
		parin.X <- c();
		if(!is.null(parX))
			parin.X <- parX;
	}
	else
	{
		parin.curve <- matrix( h1.curve$par[-c(1:NCOL(pheX))], nrow=3, byrow=T);
		parin.X <- h1.curve$par[ c(1:NCOL(pheX)) ]
	}

	for(k in 0:2)
	{
		idx.k <- which(snp.vec==k)
		if (length(idx.k ) >0 )
		{
			if(!is.null(pheX))
			{
				pheX.k<- try( pheX[idx.k, , drop = F] )
				if(class(pheX.k)=="try-error") browser();
				X  <- matrix( rep(  pheX.k %*% parin.X , NCOL(pheY) ), byrow=F, ncol=NCOL(pheY) );
			}
			else
				X <-0;

			Y.delt[idx.k,] <- pheY[idx.k,, drop=F] - get_curve(obj.phe$obj.curve, parin.curve[k+1,], pheT[idx.k,, drop=F], options=options ) - X;
		}
	}

	parin.covar <- h0$par[-c(1:(length(obj.phe$est.curve$parX)+length(obj.phe$est.curve$param)))];
	mat.cov <- get_matrix( obj.phe$obj.covar, parin.covar, pheT );
	pv <- try( dmvnorm_fast( Y.delt, rep(0,NROW(mat.cov)), mat.cov, NULL, log=T), .RR("try.silent") );

	h1.cm <- h1.curve;
	h1.cm$pv  <- pv;
	h1.cm$value <- -sum(pv);
	h1.cm$ssr <- sum(Y.delt^2, na.rm=T);
	h1.cm$par <- c(parin.X, c(t(parin.curve[snp.comb+1,])), parin.covar );
	h1.cm$R2 <- get_R2 (pheY, pheX, pheT, obj.phe$obj.curve, h1.cm$par, snp.vec);

	if( .RR("debug")  && !options$b.permu )
		cat( "H1[F]", h1.cm$value, h1.cm$ssr, "Param=", h1.cm$par, "\n");

	return(h1.cm);

}

proc_fast_h1<-function( obj.phe, pheY, pheX, pheT, snp.vec, options, h0)
{
	h1 <- proc_fast_h1_comb(obj.phe, pheY, pheX, pheT, snp.vec, options, h0, c(0,1,2) );

	pv.diff <- h1$pv - h0$pv;
	pv.neg <- which( c(sum(pv.diff[snp.vec==0]), sum(pv.diff[snp.vec==1]), sum(pv.diff[snp.vec==2]) ) < 0)
	if(length(pv.neg)==0)
		return(h1);

	h1x <- list(h1);
	h1x[[length(h1x)+1]] <- proc_fast_h1_comb(obj.phe, pheY, pheX, pheT, snp.vec, options, h0, c(0,1,1));
	h1x[[length(h1x)+1]] <- proc_fast_h1_comb(obj.phe, pheY, pheX, pheT, snp.vec, options, h0, c(0,0,2));
	h1x[[length(h1x)+1]] <- proc_fast_h1_comb(obj.phe, pheY, pheX, pheT, snp.vec, options, h0, c(0,1,0));

	idx.min <- which.min(unlist(lapply(h1x, function(h1){return(h1$value)})))
	if( h0$value - h1x[[idx.min]]$value >= 0 )
		return(h1x[[idx.min]]);

	parX <- NULL;
	if(length(obj.phe$est.curve$parX)>0)
		parX <- h0$par[1:length(obj.phe$est.curve$parX)];

	Y.delt <- matrix( NA, nrow=NROW(pheY), ncol=NCOL(pheY) );
	if( is.null(pheX) )
	{
		parin.X <- c();
		parin.curve <- h0$par[ 1:length(obj.phe$est.curve$param) ] ;
	}
	else
	{
		parin.X <- h0$par[ c(1:NCOL(pheX)) ]
		parin.curve <- h0$par[ (length(parin.X)+1):(length(parin.X) + length(obj.phe$est.curve$param) )] ;
	}

	if(!is.null(pheX))
		Y.delt<- pheY - matrix( rep(  pheX %*% parin.X , NCOL(pheY) ), byrow=F, ncol=NCOL(pheY) )
	else
		Y.delt<- pheY;

    nna.vec = get_non_na_number(pheY);

	h1.cm <- list();
	h1.cm$pv  <- rep(NA, NROW(pheY));
	h1.cm$par <- parin.X;
	parin.covar <- h0$par[-c(1:(length(obj.phe$est.curve$parX)+length(obj.phe$est.curve$param)))];
	n.par.curve	<- get_param_info( obj.phe$obj.curve, pheT, options)$count;

	optim.func <- function( parin, Y.delt, idx.k, pheT )
	{
		for(k in 0:2)
		{
			parin.curve <- parin[1:n.par.curve]
			idx.k <- which(snp.vec==k);
			if ( length(idx.k) >0 )
			{
				Y.delt[idx.k,] <- Y.delt[idx.k,,drop=F] - get_curve( obj.phe$obj.curve, parin.curve, pheT[idx.k,, drop=F], options=options );
			}
			parin <- parin[-c(1:n.par.curve)];
		}
		parin.covar <- parin;
		mat.cov <- get_matrix( obj.phe$obj.covar, parin.covar, pheT );
		pv <- try(dmvnorm_fast( Y.delt, rep(0, NCOL(mat.cov)), mat.cov, nna.vec, log=T), .RR("try.silent") );
		if ( class(pv)=="try-error" || is.na(pv) )
			return (NaN)
		else
			return(-sum(pv));
	}

	hx <- try( optim( c(rep(parin.curve,3), parin.covar), optim.func, Y.delt=Y.delt, idx.k=idx.k, pheT=pheT, method="BFGS"), .RR("try.silent") );
	if( is.null(hx) || class(hx)=="try-error" || is.na(hx$value) || is.infinite(hx$value) )
		return(list(par=c(h1.cm$par, c(rep(parin.curve,3), parin.covar )), value=NA, R2=NA));

	h1.cm$par <- c( h1.cm$par, hx$par );
	h1.cm$ssr   <- sum(Y.delt^2, na.rm=T);
	h1.cm$R2 <- get_R2 (pheY, pheX, pheT, obj.phe$obj.curve, h1.cm$par, snp.vec);
	h1.cm$value <- hx$value;

	return( h1.cm );
}
