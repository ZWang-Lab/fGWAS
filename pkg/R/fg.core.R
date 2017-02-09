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

	if(obj.phe$params$intercept)
	{
		if ( !is.null(obj.phe$pheX) )
			pheX <- cbind(1, obj.phe$pheX )
		else
			pheX <- array(1, dim=c(NROW(obj.phe$pheY),1))
	}

	options$max.time <- max(pheT, na.rm=T);
	options$min.time <- min(pheT, na.rm=T);
	options$max.loop <- 10;
	options$b.permu  <- FALSE;

	h0 <- proc_full_h0( obj.phe, pheY, pheX, pheT, options, pv.ref=NULL );

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

	h0$pv  <- try( dmvnorm_fast( Y.delt, rep(0,NROW(mat.cov)), mat.cov, nna.vec=NULL, log=T), .RR("try.slient") );
    h0$ssr <- sum(Y.delt^2, na.rm=T)

	return(h0);
}

proc_mle<-function( index, snp.info, snp.vec, obj.phe, intercept=F, options=list() )
{
	default_options = list(max.loop=5, b.permu=F);
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

	pv.ref  <- sum(-1 * obj.phe$h0$pv, na.rm=T);
	snp.miss <- which(is.na(snp.vec) | snp.vec<0 | snp.vec>2 );
	if (length(snp.miss)>0)
	{
		if (!is.null(pheX))
			pheX <- pheX[-snp.miss,,drop=F];
		pheY <- pheY[-snp.miss,,drop=F];
		pheT <- pheT[-snp.miss,,drop=F];
		snp.vec <- snp.vec[-snp.miss];
		pv.ref  <- sum( -1*obj.phe$h0$pv[-snp.miss], na.rm=T);
	}

	options$max.time=max(pheT, na.rm=T);
	options$min.time=min(pheT, na.rm=T);

	if(options$opt.method=="FAST" || options$opt.method=="FAST-NORM" )
	{
		options$max.loop<-1;
		h0 <- proc_fast_h0( obj.phe, pheY, pheX, pheT, options );
		options$max.loop<-10;
		h1 <- proc_fast_h1( obj.phe, pheY, pheX, pheT, snp.vec, options, h0);
	}
	else if(options$opt.method=="FGWAS")
	{
		options$max.loop<-1;
		h0 <- proc_full_h0( obj.phe, pheY, pheX, pheT, options, pv.ref );
		options$max.loop<-5;
		h1 <- proc_full_h1( obj.phe, pheY, pheX, pheT, snp.vec, options, h0  );
	}
	else
	{
		options$max.loop<-1;
		h0 <- proc_fast_h0( obj.phe, pheY, pheX, pheT, options );
		options$max.loop<-10;
		h1 <- proc_fast_h1( obj.phe, pheY, pheX, pheT, snp.vec, options, h0  );

		r.val <- ( h0$value - h1$value )*2;
		r.pv <- pchisq( r.val, df=(gen.par-1)*get_param_info( obj.phe$obj.curve, pheT )$count, lower.tail=FALSE );
		if( r.pv <= 0.001 )
		{
			options$max.loop<-1;
			h0 <- proc_full_h0( obj.phe, pheY, pheX, pheT, options, pv.ref );
			options$max.loop<-5;
			h1 <- proc_full_h1( obj.phe, pheY, pheX, pheT, snp.vec, options, h0  );
		}
	}

	r.val <- ( h0$value - h1$value )*2
	r.pv <- pchisq( r.val, df=(gen.par-1)*get_param_info( obj.phe$obj.curve, pheT )$count, lower.tail=FALSE );

	maf <- ifelse( mean(snp.vec)/2>0.5, 1-mean(snp.vec)/2, mean(snp.vec)/2 );

	re <- data.frame( index, snp.info[c(1:3)], matrix(c(maf, length(snp.miss), sum(snp.vec==0), sum(snp.vec==1), sum(snp.vec==2),
			gen.par, r.val, pv=r.pv, h0$value, h0$par, h1$value, h1$par ) , nrow=1 ), stringsAsFactors=F);

	if( .RR("log") >= LOG_INFO )
		cat("snp[", as.character(re[1,2]), "] loci=", re[1,3], re[1,4], "LR2=", re[1,11], re[1,12], "Param=", unlist( re[1, -c(1:12)] ), "\n");

	return(re);
}

proc_full_h0 <- function( obj.phe, pheY, pheX, pheT, options, pv.ref=NULL )
{
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

		pv <- try( dmvnorm_fast( Y.delt, rep(0,NROW(mat.cov)), mat.cov, nna.vec, log=T), .RR("try.slient") );
		if ( class(pv)=="try-error" || any(is.na(pv)) )
			return (NaN);

		A <- sum(pv);
		return(-A);
	}


	nna.vec <- get_non_na_number(pheY);
	n.par.curve <- get_param_info(obj.phe$obj.curve, pheT, options)$count;

	while( TRUE )
	{
		h0.best <- optim_BFGS( obj.phe$est.curve$parX, obj.phe$est.curve$param, obj.phe$est.covar$param, proc_h0_cov,
						pheY = pheY, pheT = pheT, pheX = pheX,
						obj.curve = obj.phe$obj.curve, obj.covariance = obj.phe$obj.covar,
						nna.vec = nna.vec, n.par.curve = n.par.curve, options = options);

		if( !is.null(pv.ref) && h0.best$value > pv.ref  * 1.01 )
		{
#cat("H0=",h0.best$value, "REF=", pv.ref, "\n");
			options$max.loop <- options$max.loop * 2;
			reset_seed();
		}
		else
			break;
	}

	if( .RR("log")>= LOG_INFO && !options$b.permu )
		cat( "H0[F]", h0.best$value, h0.best$par, "\n");

	return(h0.best);
}

proc_full_h1<-function( obj.phe, pheY, pheX, pheT, snp.vec, options, h0)
{
	proc_only_mu<-function( parin, pheY, pheT, pheX, snp.vec, obj.curve, sum.type )
	{
		mu  <- get_curve( obj.curve, parin[-c(1:NCOL(pheX))], pheT, options=options )
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

		Y.delt <- matrix( NA, nrow=NROW(pheY), ncol=NCOL(pheY) );
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
					X <- 0;

				parin.curve0 <- parin.curve[1:n.par.curve];
				parin.curve  <- parin.curve[-c(1:n.par.curve)];

				Y.delt[idx.k,] <- pheY[idx.k,, drop=F] - get_curve(obj.curve, parin.curve0, pheT[idx.k,, drop=F], options=options ) - X;
			}

		}

		mat.cov <- get_matrix( obj.covar, parin.covar, pheT );
		pv <- try(dmvnorm_fast( Y.delt, rep(0, NCOL(mat.cov)), mat.cov, nna.vec=nna.vec, log=T), .RR("try.slient") );
		if ( class(pv)=="try-error" || is.na(pv) )
			return (NaN);

		return( - sum( pv ) );
	}

	#par.fin <- c();
	#n.gentype <- 0;
	#for(k in 0:2)
	#{
	#	if (length( which(snp.vec==k)) >0 )
	#	{
	#		parin <- c(obj.phe$est.curve$parX, obj.phe$est.curve$param );
	#		h1c   <- optim_least_square( parin, proc_only_mu,
	#		                    pheY = pheY[which(snp.vec==k),,drop=F],
	#							pheT = pheT[which(snp.vec==k),,drop=F],
	#							pheX = pheX[which(snp.vec==k),,drop=F],
	#							snp.vec = NULL,
	#						  	obj.curve = obj.phe$obj.curve );
	#
	#		par.fin <- rbind( par.fin, h1c$par[-(1:length(obj.phe$est.curve$parX))] );
	#		n.gentype <- n.gentype + 1;
	#	}
	#}
	#
	#if( .RR("log") >= LOG_INFO  && !options$b.permu ) show(par.fin);
	#parin   <- c( obj.phe$est.curve$parX, c(t(par.fin)), obj.phe$est.covar$param );
	#n.par.curve <- get_param_info(obj.phe$obj.curve, pheT, options)$count;
	#nna.vec <- get_non_na_number( pheY );

	nna.vec <- get_non_na_number( pheY );
	n.gentype <- length(unique(snp.vec));
	n.par.curve <- get_param_info(obj.phe$obj.curve, pheT, options)$count;
	parin.x     <- h0$par[ 1:NCOL(pheX) ]
	parin.curve <- rep(h0$par[ (NCOL(pheX)+1):(NCOL(pheX)+n.par.curve) ], n.gentype );
	parin.covar <- h0$par[ -(1:(NCOL(pheX)+n.par.curve)) ] ;

	h1.best <-  optim_BFGS ( parin.x, parin.curve, parin.covar,
					  proc_h1_cov, pheY = pheY, pheT = pheT, pheX = pheX,
					  obj.curve = obj.phe$obj.curve, obj.covar = obj.phe$obj.covar,
					  nna.vec=nna.vec, snp.vec=snp.vec,
					  n.par.curve=n.par.curve, n.gentype=n.gentype, options=options );

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
	if( .RR("log") >= LOG_INFO  && !options$b.permu )
		cat( "H1[F]", h1.best$value, h1.best$par, "\n");

	return(h1.best);

}

optim_least_square<-function ( parin, fn_mle, pheY, pheT, pheX, snp.vec=NULL, obj.curve=NULL, method=NULL, ssr.ref=NULL )
{
	f_optim <- function()
	{
		parinx<-  parin* runif( length(parin), 0.9, 1.1 );
		control <- list(maxit = 500 + length(parin) * 200, reltol=1e-8);
		loop.optim <- 1;

		while(TRUE)
		{
			h0 <- try( optim( parinx, fn_mle, pheY = pheY, pheT = pheT, pheX = pheX, snp.vec = snp.vec,
							obj.curve = obj.curve, sum.type=0,
							method = ifelse(loop.optim%%2==0, "Nelder-Mead", "BFGS" ),
							control= control ),
							.RR("try.slient") );

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
						cat("H0_1, convergence=", h0$convergence, "", h0$message, "\n")
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
						pheY = pheY, pheT = pheT, pheX = pheX, snp.vec = snp.vec,
						obj.curve = obj.curve, 	sum.type=1,
						control = nls.lm.control(nprint=0, maxiter=control$maxiter ) ),
						.RR("try.slient") );

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
						cat("nls.lm, info=", h0.nls$info, "", h0.nls$message, "\n")
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

	if(method!="FAST-NORM")
		A1 <- f_minpack.lm()
	else
		A1 <- f_optim();

	return( A1 );

}

colVar <- function(M )
{
	apply( M, 2, function(mc) {var(mc, na.rm=T)} )
}

proc_fast_h0 <- function( obj.phe, pheY, pheX, pheT, options)
{
	proc_h0_curve<-function( parin, pheY, pheT, pheX, snp.vec, obj.curve, sum.type )
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
		    #Y.var <- colVar(Y.delt);
		    #A <- NROW(Y.delt)*sum(log(sqrt(1/Y.var))) - 0.5*sum( rowSums((Y.delt^2)/matrix(Y.var, ncol=NCOL(Y.delt), nrow=NROW(Y.delt), byrow=T), na.rm=T) );
		    #A <- -A;

			Y.var <- colVar(Y.delt);
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

		pv <- try( dmvnorm_fast( Y.delt, rep(0,NROW(mat.cov)), mat.cov, nna.vec=nna.vec, log=T), .RR("try.slient") );
		if ( class(pv)=="try-error" || any(is.na(pv)) )
			return (NaN);

		A <- sum(pv);
		return(-A);
	}

	parinx <- c( obj.phe$est.curve$parX, obj.phe$est.curve$param );
	h0.curve <- optim_least_square( parinx, proc_h0_curve, pheY = pheY , pheT = pheT, pheX = pheX, snp.vec= NULL, obj.curve = obj.phe$obj.curve, method=(options$opt.method=="FAST")  );

	X <- 0;
	parin.curve <- h0.curve$par;
	if (!is.null(pheX))
	{
		parin.X <- h0.curve$par[1:NCOL(pheX)];
		parin.curve <- h0.curve$par[-c(1:NCOL(pheX))];
		X  <- matrix( rep(  pheX %*% parin.X , NCOL(pheY) ), byrow=F, ncol=NCOL(pheY) );
	}

	mu <- get_curve( obj.phe$obj.curve, parin.curve, pheT, options=options);
    Y.delt <- pheY - mu - X;

	nna.vec = get_non_na_number( pheY );
	h0.cm <- optim_BFGS( c(), c(), obj.phe$est.covar$param,
					proc_h0_covmatrix, Y.delt=Y.delt, pheT=pheT, obj.covariance = obj.phe$obj.covar, nna.vec = nna.vec, options=options);
	mat.cov <- get_matrix( obj.phe$obj.covar, h0.cm$par, pheT );
	h0.cm$pv  <- try( dmvnorm_fast( Y.delt, rep(0,NROW(mat.cov)), mat.cov, nna.vec=nna.vec, log=T), .RR("try.slient") );
    h0.cm$par <- c( h0.curve$par, h0.cm$par);
    h0.cm$ssr <- sum(Y.delt^2, na.rm=T)

	if( .RR("log")>= LOG_INFO && !options$b.permu )
		cat( "H0[F]", h0.cm$ssr, h0.cm$value, h0.cm$par, "\n");

	return(h0.cm);
}

proc_fast_h1_comb<-function( obj.phe, pheY, pheX, pheT, snp.vec, options, h0, snp.comb=c(0,1,2) )
{
	proc_h1_curve<-function( parin, pheY, pheT, pheX, snp.vec, obj.curve, sum.type)
	{
		if(!is.null(pheX))
		{
			parin.X <- parin[1:NCOL(pheX)];
			parin.curve <- parin[-c(1:NCOL(pheX))];
		}
		else
			parin.curve <- parin;

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
		    #Y.var <- colVar(Y.delt);
		    #A <- NROW(Y.delt)*sum(log(sqrt(1/Y.var))) - 0.5*sum( rowSums((Y.delt^2)/matrix(Y.var, ncol=NCOL(Y.delt), nrow=NROW(Y.delt), byrow=T), na.rm=T) );
		    #A <- -A;

			Y.var <- colVar(Y.delt);
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

		pv <- try(dmvnorm_fast( Y.delt, rep(0, NCOL(mat.cov)), mat.cov, nna.vec, log=T), .RR("try.slient") );
		if ( class(pv)=="try-error" || is.na(pv) )
			return (NaN);

		return( - sum(pv) );
	}

	parin    <- c( obj.phe$est.curve$parX, rep( obj.phe$est.curve$param, 3 ) );

	snp.vec <- snp.comb[snp.vec+1];
	h1.curve <- optim_least_square( parin, proc_h1_curve, pheY = pheY, pheT = pheT, pheX = pheX, snp.vec=snp.vec, obj.curve = obj.phe$obj.curve, method=(options$opt.method=="FAST"), ssr.ref=h0$ssr );

	Y.delt <- matrix( NA, nrow=NROW(pheY), ncol=NCOL(pheY) );
	parin.curve <- matrix( h1.curve$par, nrow=3, byrow=T);
	parin.X <- c();
	if(!is.null(pheX))
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

	#nna.vec = get_non_na_number(pheY);
	#h1.cm <- optim_BFGS( c(), c(), obj.phe$est.covar$param, proc_h1_cov, Y.delt = Y.delt, pheT = pheT, obj.covariance = obj.phe$obj.covar, nna.vec = nna.vec, options = options);
    #h1.cm$par <- c( h1.curve$par, h1.cm$par);

	parin.covar <- h0$par[-c(1:(length(obj.phe$est.curve$parX)+length(obj.phe$est.curve$param)))];
	mat.cov <- get_matrix( obj.phe$obj.covar, parin.covar, pheT );
	pv <- try( dmvnorm_fast( Y.delt, rep(0,NROW(mat.cov)), mat.cov, NULL, log=T), .RR("try.slient") );

	h1.cm <- h1.curve;
	h1.cm$pv  <- pv;
	h1.cm$value <- -sum(pv);
	h1.cm$ssr <- sum(Y.delt^2, na.rm=T);
	h1.cm$par <- c(parin.X, c(t(parin.curve[snp.comb+1,])), parin.covar );

	if( .RR("log") >= LOG_INFO  && !options$b.permu )
		cat( "H1[F]", h1.cm$ssr, h1.cm$value, h1.cm$par, "\n");

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
	return(h1x[[idx.min]]);
}
