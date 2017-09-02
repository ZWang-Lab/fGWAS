fg_dat_est<-function( obj.phe, curve.type="auto", covariance.type="auto", file.plot.pdf=NULL, options=list() )
{
	default_options <- list( max.optim.failure=100, min.optim.success=20, R2.loop=5, verbose=FALSE);
	default_options[names(options)] <- options;
	options <- default_options;

	intercept <- ifelse( is.null(obj.phe$intercept),  FALSE, obj.phe$intercept );
	pheX <- obj.phe$pheX;
	if(intercept)
	{
		if(is.null(pheX))
			pheX <- matrix(1, nrow=NROW(obj.phe$pheY), ncol=1)
		else
			pheX <- cbind(1, pheX );
	}

	r.test <- NULL;
	if( is.null( obj.phe$obj.curve ) || toupper(obj.phe$obj.curve@type) != toupper(curve.type))
	{
		cat("  No curve or new curve type is specified, curve fitting is being performed.\n");
		r.test <- fg_fit_curve( obj.phe$pheY, pheX, obj.phe$pheT, curve.type=curve.type, file.plot.pdf=file.plot.pdf, options=options );

		if(r.test$error)
			stop("? No curve is fitted to this data.");

		obj.phe$obj.curve <- fg.getCurve( r.test$type );
		obj.phe$summary.curve <- r.test;
	}

	cat("  Curve Type==> ", obj.phe$obj.curve@type, " <==\n");
	par.init <- if(!is.null(r.test)) r.test$par else NULL;
	r <- proc_est_curve( obj.phe$pheY, pheX, obj.phe$pheT, obj.phe$obj.curve, par.init= par.init, options=options )
	if( r$error )
		stop("? Can't estimate the parameters of mean vector according to the curve function" );

	parX.len <- NCOL(pheX)
	if(is.null(pheX)) parX.len <- 0;

	cat(" Parameter range estimation ...... \n")
	options$max.optim.failure <- 10;
	options$min.optim.success <- 2;
	range <- proc_est_curve_range(obj.phe$pheY, pheX, obj.phe$pheT, obj.phe$obj.curve, par.init = r$par, options=options);
	obj.phe$est.curve <- list( type = obj.phe$obj.curve@type,
	                           intercept = intercept,
							   param = if(parX.len==0) r$par else r$par[-c(1:parX.len)],
							   param.lower = if(parX.len==0) r$par else range$lower[-c(1:parX.len)],
							   param.upper = if(parX.len==0) r$par else range$upper[-c(1:parX.len)],
							   parX  = if(parX.len==0) c() else r$par[1:parX.len],
							   parX.lower = if(parX.len==0) c() else range$lower[1:parX.len],
							   parX.upper = if(parX.len==0) c() else range$upper[1:parX.len],
							   R2 = r$R2);

	if( r$R2 < - 0.1 )
		cat("! The R2 value for curve fitting indicates the curve type is inappropriate for the phenotype data.(R2=", r$R2,")\n")

	if( is.null(obj.phe$obj.covar) || ( toupper(obj.phe$obj.covar@type) != toupper(covariance.type) ) )
	{
		r.test <- fg_fit_covar( obj.phe$pheY, pheX, obj.phe$pheT, r$y.resd, obj.phe$obj.curve, covariance.type, options=options );

		if(r.test$error)
			stop("? No covariance is fitted to this data.")
		else
		{
			obj.phe$summary.covariance <- r.test;
			obj.phe$obj.covar <- fg.getCovariance( r.test$type );
		}
	}

	cat("  Covariance Type==> ", obj.phe$obj.covar@type, " <==\n");

	options$max.optim.failure <- 50;
	options$min.optim.success <- 5;
	r.est <- proc_est_covar( r$y.resd, NULL, obj.phe$pheT, obj.phe$obj.curve, obj.phe$obj.covar, options=options );
	if ( r.est$error )
		stop("? Can't estimate the parameters of covariance structure according to the curve function" );

	obj.phe$est.covar<- list( type = obj.phe$obj.covar@type, param = r.est$par);
	obj.phe$error <- F;

	return(obj.phe);
}

fn_get_resd<-function(pheY, pheX, pheT, obj.curve, parin  )
{
	par_X <- c()
	par_c <- parin;
	if ( !is.null(pheX) )
	{
		par_X <- parin[ 1:NCOL(pheX)];
		par_c <- parin[ -c(1:NCOL(pheX)) ]
	}

	mu_gen <- get_curve( obj.curve, par_c, pheT, options=list(max.time=max(pheT, na.rm=T), min.time=min(pheT, na.rm=T)) )
	if(all(is.na( mu_gen )))
		return(NULL);

	if(!is.null(pheX) )
		X <- matrix( rep( pheX %*% par_X, NCOL(pheY)), ncol=NCOL(pheY), byrow=F)
	else
		X <- 0;

	if( is.vector( mu_gen ) )
		y_resd <- t(t(pheY) - mu_gen) -  X
	else
		y_resd <- pheY - mu_gen  - X;

	return( y_resd );
}

get_R2<-function(pheY, pheX, pheT, obj.curve, h.par, snp.vec=NULL )
{
	pheY_hat <- matrix(NA, nrow=NROW(pheY), ncol=NCOL(pheY) );
	pheY_hat_vec <- c(pheY_hat);

	pheT_vec <- c(pheT);
	pheY_vec <- c(pheY);
	for (ti in unique(pheT_vec) )
	{
		if (!is.na(ti))
			pheY_hat_vec[which(pheT_vec==ti)] <- mean(pheY_vec[which(pheT_vec==ti)], na.rm=T)
	}

	pheY_hat <- matrix(pheY_hat_vec, nrow=NROW(pheY), ncol=NCOL(pheY) );

	if(is.null(snp.vec))
		y.resd <- fn_get_resd( pheY, pheX, pheT, obj.curve, h.par )
	else
	{
		y.resd <- c();
		par.X <- if(is.null(pheX)) NULL else h.par[1:NCOL(pheX)];
		par.curve <- if(is.null(pheX)) h.par else h.par[-c(1:NCOL(pheX))];
		len.curve <- get_param_info( obj.curve, pheT)$count;

		for(i in 0:2)
		{
			idx <- which(snp.vec==i);
			if(length(idx)>0)
			{
				hx.par <- c(par.X, par.curve[1:len.curve]);
				par.curve <- par.curve[-c(1:len.curve)];
				y.resd <- rbind(y.resd, fn_get_resd( pheY[idx,,drop=F], pheX[idx,,drop=F], pheT[idx,,drop=F], obj.curve, hx.par ) );
			}
		}
	}

	#R2   <- sum(((pheY-y.resd)-pheY_hat)^2, na.rm=T)/sum((pheY-pheY_hat)^2, na.rm=T)
	R2   <- 1 - sum(y.resd^2, na.rm=T)/sum((pheY-pheY_hat)^2, na.rm=T)

	return(R2);
}

proc_est_curve<-function(  pheY, pheX, pheT, obj.curve, par.init=NULL, options=NULL  )
{
	options$max.time <- max(pheT, na.rm=T)
	options$min.time <- min(pheT, na.rm=T)

	get_init_covariate_par <- function(pheY, pheX)
	{
		colnames(pheX) <- paste("x", 1:NCOL(pheX), sep="");
		fit <- lm(Y ~.-1, data=data.frame(Y=rowMeans(pheY), pheX) );
		if(!is.null(fit$coefficients))
			return(fit$coefficients)
		else
			return(mean(pheY, na.rm=T)/colMeans(pheX, na.rm=T));
	}

	get_init_curve_par<-function( pheY, pheX, pheT, f.obj )
	{
		par.X <- c(); YX <- 0;
		if( !is.null(pheX) )
		{
			par.X <- get_init_covariate_par(pheY, pheX)
			YX <- matrix( rep( pheX %*% par.X, NCOL(pheY)), ncol=NCOL(pheY), byrow=F);
		}

		par.curve <- est_init_param( f.obj, pheY, pheX, pheT, options=options );
		return( c( par.X, par.curve)  );
	}

	get_rand_curve_par<-function( pheY, pheX, pheT, f.obj, parin )
	{
		uc <- check_fittness( pheY, pheX, pheT, f.obj, parin )

		if(uc$fit)
		{
			return( parin*runif(length(parin), 0.9, 1.1));
		}
		else
		{
			if ( uc$over0.05 < NCOL(pheY)/4 )
				return( parin * runif(length(parin), 0.5, 1.5))
			else
			{
				par.curve <- get_init_curve_par( pheY, pheX, pheT, f.obj );
				return( par.curve * runif( length(par.curve) , 0.5, 1.5 ) );
			}
		}
	}

	# parin : par.X[1..n], par.curve.
	fn_mle_est<-function( parin, extra_par )
	{
		y_resd <- fn_get_resd( extra_par$pheY, extra_par$pheX, extra_par$pheT, extra_par$obj.curve, parin );
		if(is.null(y_resd))
			return(NaN);

		A <- sum( y_resd^2, na.rm=T );

		return(A);
	}


	if( is.null( par.init) )
		par.init <- get_init_curve_par(pheY, pheX, pheT, obj.curve);

	if(.RR("debug")) cat("Initial parameters for curve: ", par.init, "\n")

	R2.max <- -Inf;
	R2.success <- 0;
	R2.failed <- 0
	R2.list <- NULL;
	while ( R2.success < ifelse( is.null(options$R2.loop), 5, options$R2.loop ) && R2.failed <= 500  )
	{
		reset_seed();
		h0 <- proc_mle_loop(  pheY, pheX, pheT,
				obj.curve,
				fn_mle_est,
				mle_extra_par = list( pheY=pheY, pheX=pheX, pheT=pheT, obj.curve=obj.curve ),
				parin = par.init,
				fn_init_par = get_init_curve_par,
				fn_rand_par = get_rand_curve_par,
				options=options )

		if(is.null(h0) || is.na(h0$value) || is.infinite(h0$value))
		{
			R2.failed <- R2.failed + 1;
			next;
		}
		else
		{
			R2 <- get_R2(pheY, pheX, pheT,  obj.curve, h0$par );
if(.RR("debug")) cat( "MLE results for curve: R2=", R2, h0$value, h0$par, "\n");

			if( R2 > 1 || R2 < -0.1 )
				R2.failed <- R2.failed + 1
			else
				R2.success <- R2.success + 1;

			if( R2 <= R2.max ) next;

			R2.max <- R2;
			par.init <- h0$par*runif(length(h0$par), 0.99, 1);

			y.resd <- fn_get_resd( pheY, pheX, pheT,  obj.curve, h0$par );
			## The sum of the squared differences between each observation and its predicted value.
			y.SSE <- sum(y.resd^2,na.rm=T) ;
			##Gives the average of the squares of the errors of each value.
			y.MSE <- y.SSE/length(which(!is.na(pheY)));
			##The square root of the MSE that estimates the standard deviation of the random error.
			y.RMSE <- sqrt(y.MSE)
			## pramater count
			K <- get_param_info(obj.curve, pheT)$count;
			n.sample <- NROW(pheY);

			AIC  <- 2*K +n.sample*log(y.SSE/n.sample)
			AICc <- log(y.SSE/n.sample) + (n.sample+K)/(n.sample-K-2)
			BIC  <- n.sample*log(y.SSE/n.sample) + K*log(n.sample)

			R2.list <- list(error=F, par.count = K, AIC = AIC, AICc = AICc, BIC = BIC,
							  SSE = y.SSE, MSE = y.MSE, RMSE = y.RMSE, R2 = R2,
							  par=h0$par, y.resd=y.resd)
		}
	}

	if (is.null(R2.list))
		return(list(error=T, par=NA, val=NA, R2=NA))
	else
		return(R2.list);
}

proc_est_curve_range<-function( pheY, pheX, pheT, f.curve, par.init, options=NULL )
{
	n.obs <- NROW( pheY );
	mu.pars <- c();

	loop <- 0;
	while( loop < .RR("mu.range.loop")  )
	{
		y.sub <- sample(n.obs)[1:round( n.obs * runif(1,0.5,0.9) )]

		pheX0 <- NULL;
		pheY0 <- pheY[y.sub,,drop=F];
		if( !is.null(pheX) ) pheX0 <- pheX[y.sub,,drop=F];
		if( !is.null(dim(pheT)) )  pheT0 <- pheT[y.sub,,drop=F] else pheT0 <- pheT;

		r <- proc_est_curve( pheY0, pheX0, pheT0, f.curve, par.init, options=options );
		if (r$error) next;

		mu.pars <- rbind(mu.pars, r$par);
		loop <- loop + 1;
	}

	mu.lower <-c();
	mu.upper <-c();

	for(i in 1:NCOL(mu.pars) )
	{
		mu.lower <- c(mu.lower, min(mu.pars[,i]) )
		mu.upper <- c(mu.upper, max(mu.pars[,i]) )
	}

	return(list(lower=mu.lower, upper=mu.upper))
}


proc_est_covar<-function( Y.resd, pheX, pheT, obj.curve, obj.covar, par.init=NULL, options=NULL  )
{
	get_init_covar_par<-function( Y.resd, pheX, pheT, f.covar)
	{
		return( est_init_param( f.covar, Y.resd, pheX, pheT, options=options ) );
	}

	get_rand_covar_par<-function( Y.resd, pheX, pheT, f.covar, parin )
	{
		return( parin* runif(length(parin), 0.9, 1.1) );
	}

	#parin:
	# phi1, s1, phi2, s2
	fn_mle_est<-function( parin, extra_par)
	{
		y.resd  <- extra_par$Y.resd;
		pheX    <- extra_par$pheX;
		pheT    <- extra_par$pheT;
		f.covar <- extra_par$f.covar;

		cov.mat  <- get_matrix(f.covar, parin, pheT);
		if(any(is.na(cov.mat))) return(NaN);

		pv <- dmvnorm_fast( y.resd, rep(0, NCOL(y.resd)), cov.mat, log=T);
		if(any(is.infinite(pv)))
			return(NaN);

		A <- sum( pv );
		return( -A );
	}

	if( is.null( par.init) )
		par.init <- get_init_covar_par(  Y.resd, pheX, pheT, obj.covar );

	if(.RR("debug"))  cat( "Initial parameters for covariance: ", par.init, "\n");

	h0 <- proc_mle_loop( Y.resd, pheX, pheT,
			obj.covar,
			fn_mle_est,
			mle_extra_par=list(Y.resd=Y.resd, pheX=pheX, pheT=pheT, f.covar=obj.covar),
			parin = par.init,
			fn_init_par = get_init_covar_par,
			fn_rand_par = get_rand_covar_par,
			options=options )

	if(.RR("debug"))  cat( "MLE results for curve covariance:", h0$value, h0$par, "\n");

	if( is.finite( h0$value ) )
	{
		K <- get_param_info(obj.covar, pheT)$count;
		###Note: L=-log(L)
		AIC  <- 2*K + 2*h0$value;
		BIC  <- 2*h0$value + K*log(NROW(Y.resd))
		return(list(error=F, AIC=AIC, BIC=BIC, par=h0$par, logL = -h0$value ))
	}
	else
		return(list(error=T, par=NA, err.info="Failed to estimate covariance structure."));
}

proc_mle_loop<-function(  pheY, pheX, pheT, f.obj, fn_mle, mle_extra_par, parin, fn_init_par, fn_rand_par, options=list(min.optim.success=5, max.optim.failure=100)  )
{
	control <- list(maxit = 50000, reltol=1e-8);
	h0      <- list( value=Inf, par=parin );

	init.optim <- 1
	while( is.infinite(h0$value) )
	{
		try( h0 <- optim( parin, fn_mle, extra_par=mle_extra_par,
				method = ifelse(runif(1)>0.75, "Nelder-Mead", "BFGS"),
				control= control ), .RR("try.silent", FALSE)  );

		if ( class(h0)=="try-error" || is.na(h0) || is.infinite(h0$value) || is.null(h0$convergence) || (h0$convergence !=0) )
		{
			if ( !is.na(h0) && is.list(h0) && !is.null(h0$convergence) )
			{
				if (h0$convergence == 1)
				{
					control$maxit <- control$maxit*2;
					if ( control$maxit > 500*4096 )
					{
						control$maxit <- 500*4096;
						control$reltol <- control$reltol*2;
					}
				}
				else
					if(.RR("debug")) cat("h0$convergence =", h0$convergence , "\n");
			}

			reset_seed();
			parin <- fn_init_par( pheY, pheX, pheT, f.obj );

			init.optim <- init.optim + 1;
			if(init.optim > 100) return( list(error=T, par=NA, value=NA) );

			next;
		}

		reset_seed();
	}

	if(.RR("debug")) cat( "MLE[0]", h0$value, h0$par, "\n");

	parin0 <- h0$par;
	n.optim.failure <- 0;
	n.optim.success <- 0;

	while ( n.optim.failure < options$max.optim.failure && n.optim.success < options$min.optim.success )
	{
		parinx<- fn_rand_par( pheY, pheX, pheT, f.obj, parin0 );
		h2 <- NULL;
		try( h2 <- optim( parinx, fn_mle, extra_par=mle_extra_par,
				method = ifelse((n.optim.failure+n.optim.success)%%2==1, "Nelder-Mead", "BFGS"),
				control=control ), .RR("try.silent", FALSE)  );

		if (is.null(h2) || class(h2)=="try-error" || any(is.na(h2)) || (h2$convergence !=0)  )
		{
			if ( is.list(h2) && h2$convergence == 1)
			{
				control$maxit <- control$maxit*2;
				if ( control$maxit > 500*4096 )
				{
					control$maxit <- 500*4096;
					control$reltol <- control$reltol*2;
				}
			}
			n.optim.failure <- n.optim.failure+1;
		}
		else
		{
			if( h2$value < h0$value)
				h0 <- h2;

			n.optim.success <- n.optim.success + 1;
		}

		reset_seed();
	}

	if(.RR("debug")) cat( "MLE[F]", h0$value, h0$par, "\n");
	return(list(error=F, par=h0$par, value=h0$value))
}

check_fittness<-function( pheY, pheX, pheT, f.obj, parin )
{
	if( !inherits(f.obj, "fg.curve.base") )
		return(list(fit=TRUE));

	y_resd <- fn_get_resd( pheY, pheX, pheT, f.obj, parin );
    if( is.null(y_resd) )
   		return(list(fit=TRUE));

if(NCOL(y_resd)==1) browser();

	y_sd   <- apply( pheY, 2, function(x) sd(x) )

	py.prob  <- pnorm( colMeans(y_resd), mean=0, sd=y_sd, lower.tail=F );
	py.05 <- length( which(py.prob<0.05) );

	return(list(fit = ifelse( py.05 > 1, F, T), over0.05=py.05) );
}

fg_fit_curve<-function( pheY, pheX, pheT, curve.type="auto", file.plot.pdf=NULL, options=options )
{
	cat(" Searching curve type ......\n" );

	obj.curves <- list();
	if(curve.type=="auto" || curve.type=="" || is.null(curve.type) )
	{
		for ( i in 1:fg_get_curve_count())
			obj.curves[[i]] <- fg.getCurve( i );
	}
	else
	{
		for(i in 1:length(curve.type))
			obj.curves[[i]] <- fg.getCurve(curve.type[i]);
	}

	if(!is.null(file.plot.pdf)) pdf(file.plot.pdf);

	est.curve <- list();
	index = 1
	for ( obj.curve in obj.curves)
	{
		cat( "* [", index, "/", length(obj.curves), "] try curve: ", obj.curve@type, "\n" );
		r <- proc_est_curve( pheY, pheX, pheT, obj.curve, options=options )
		if( r$error )
			warning("Can't estimate the parameters of mean vector according to the curve function[ curve.type=", obj.curve@type, "]" )
		else
		{
			r$type <- obj.curve@type;
			est.curve[[index]] <- r;
			index<-index+1;

			if(!is.null(file.plot.pdf))
			{
				plot(1,1, type="n", xlim=c(min(pheT, na.rm=T), max(pheT, na.rm=T)), ylim=c(min(pheY, na.rm=T), max(pheY, na.rm=T)), main=r$type, xlab="Time", ylab="Phenotype");
				for(i in 1:NROW(pheY))
					lines(pheT[i,], pheY[i,], col="gray", lwd=0.2);

				pheT_vec <- c(pheT)
				pheY_vec <- c(pheY)
				ti <- unique(c(pheT));
				ti <- sort( ti [ !is.na(ti) ] );
				y.mu <- rep(NA, length(ti));
				for(i in 1:length(ti) )
				   y.mu[i] <-  mean(pheY_vec[pheT_vec==ti[i]], na.rm=T);

				lines(ti, y.mu, col="black", lty="22", lwd=1);

				mu_X <- 0 ;
				par_c <- r$par;
				if ( !is.null(pheX) )
				{
				   par_X <- c( par_X, r$par[ 1:NCOL(pheX)]);
				   par_c <- r$par[ -c(1:NCOL(pheX)) ]
				   mu_X <- mean( pheX %*% par_X, na.rm=T );
				}

				ti <- seq( min(pheT, na.rm=T), max(pheT, na.rm=T), 1);
				mu_gen <- get_curve( obj.curve, par_c, ti, options=list(max.time=max(pheT, na.rm=T), min.time=min(pheT, na.rm=T)))

				lines(ti, mu_gen + mu_X, col="black", lwd=1.5);
			}
		}
	}

	if(!is.null(file.plot.pdf)) dev.off();

	if(length(est.curve)>0)
		est.summary = data.frame(type = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$type);     }) ),
							 parm = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$par.count);}) ),
							 AIC  = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$AIC);      }) ),
							 AICc = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$AICc);     }) ),
							 BIC  = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$BIC);      }) ),
							 SSE  = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$SSE);      }) ),
							 MSE  = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$MSE);      }) ),
							 RMSE = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$RMSE);     }) ),
							 R2   = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$R2);       }) ) )
	else
	{
		est.summary = NULL;
		stop("Failed to do curve fitting.\n");
	}

	## use AIC to determine the best curve type.
	cat("  Curve Fitting Summary:\n");
	old.width <- options(width=132)
	show(est.summary);
	options(width=old.width$width)

	fit.idx <- which.min( est.summary[,3] );
	return(list(error=F, type=est.curve[[fit.idx]]$type, y.resd=est.curve[[fit.idx]]$y.resd, par = est.curve[[fit.idx]]$par, summary=est.summary, est.curve=est.curve) );
}

fg_fit_covar<-function( pheY, pheX, pheT, y.resd, obj.curve, covariance.type="auto", options=NULL)
{
	cat(" Searching covariance matrix .......\n" );

	obj.covars <- list();
	if(covariance.type=="auto" || covariance.type=="" || is.null(covariance.type) )
	{
		for ( i in 1:fg_get_covar_count()) obj.covars[[i]] <- fg.getCovariance( i );
	}
	else
	{
		for(i in 1:length(covariance.type))
			obj.covars[[i]] <- fg.getCovariance(covariance.type[i]);
	}

	est.covar <- list();
	index = 1
	for ( obj.covar in obj.covars)
	{
		cat(" *[", index, "/", fg_get_covar_count(), "]: try covariance matrix: ", obj.covar@type, "\n" );

		r.est <- proc_est_covar( y.resd, NULL, pheT, obj.curve, obj.covar, options = options );
		if( r.est$error )
			warning("Can't estimate the parameters of covariance structure according to the covariance function[ covariance.type=", obj.covar@type, "]" )
		else
		{
			r.est$type <- obj.covar@type;
			est.covar[[index]] <- r.est;
			index <- index+1;
		}
	}

	if(length(est.covar)>=1)
	{
		est.summary = data.frame(type = unlist( lapply(1:length(est.covar), function(i){return(est.covar[[i]]$type); })),
							 L    = unlist( lapply(1:length(est.covar), function(i){return(est.covar[[i]]$logL[1]); })),
							 AIC  = unlist( lapply(1:length(est.covar), function(i){return(est.covar[[i]]$AIC); })),
							 BIC  = unlist( lapply(1:length(est.covar), function(i){return(est.covar[[i]]$BIC); })));


		cat("  Covariance Estimation Summary:\n");
		old.width <- options(width=132)
		show(est.summary);
		options(width=old.width$width)

		fit.idx <- which.min(est.summary[,3]);
		return(list(error=F, type=as.character(est.covar[[fit.idx]]$type), par = est.covar[[fit.idx]]$par, summary=est.summary));
	}
	else
		return(list(error=T));
}
