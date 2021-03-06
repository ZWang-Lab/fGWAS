## This common library for Funmap and fGWAS
## The newest source is in fGWAS library
## !!!! ALL modifies in fGWAS !!!!

##-----------------------------------------------------------
##
##
##-----------------------------------------------------------

setClass("fg.curve.base",
	representation(
		type           = "character",    #curve_type
		description    = "character"
  )
)

setMethod("show", signature(object="fg.curve.base"), function(object){
   cat("     Class :", class(object), "\n");
   cat("Curve Type :", object@type, "\n");

   info <- get_param_info(object, NULL );
   cat("Parameters :", info$names, "\n");
   cat("   Formula :", info$formula, "\n");
});

#------------------------------------------------------------------------------------------
# mufunc
#------------------------------------------------------------------------------------------

fg.getCurve<-function(type)
{
	obj.curve <- NULL;
	if(is.character(type))
	{
		for(obj in fg.allCurves())
			if(toupper(obj@type)==toupper(type))
				obj.curve <- obj;
	}

	if(is.numeric(type))
		obj.curve <- fg.allCurves()[[type]];

	return(obj.curve);
}

fg_get_curve_count<-function()
{
	return(length(fg.allCurves()));
}

fg.allCurves<-function()
{
	list.curve <- list()
	cn<-0
	list.curve[[cn<-cn+1]] <- new("fg.curve.log",   type = "Logistic",   description = "logistic curve");
	list.curve[[cn<-cn+1]] <- new("fg.curve.log2",  type = "Bi-Logistic",  description = "Double logistic curve");
	list.curve[[cn<-cn+1]] <- new("fg.curve.abrk",  type = "ABRK",  description = "ABRK model");

	list.curve[[cn<-cn+1]] <- new("fg.curve.pharmacology",  type = "Pharmacology",  description = "Pharmacology curve");
	list.curve[[cn<-cn+1]] <- new("fg.curve.exponential",  type = "Exponential",  description = "Exponential curve");
	list.curve[[cn<-cn+1]] <- new("fg.curve.bi.exponential",  type = "Bi-Exponential",  description = "Bi-exponential curve");

	list.curve[[cn<-cn+1]] <- new("fg.curve.power",  type = "Power",  description = "power curve");

	list.curve[[cn<-cn+1]] <- new("fg.curve.legendre2",  type = "Legendre2",  description = "Legendre Polynomial(2nd-order)");
	list.curve[[cn<-cn+1]] <- new("fg.curve.legendre3",  type = "Legendre3",  description = "Legendre Polynomial(3rd-order)");
	list.curve[[cn<-cn+1]] <- new("fg.curve.legendre4",  type = "Legendre4",  description = "Legendre Polynomial(4th-order)");

	#list.curve[[cn<-cn+1]] <- new("fg.curve.ChapmanRichard",  type = "ChapmanRichard",  description = "Chapman-Richard");

	return(list.curve);
}

fg.addCurve<-function()
{

}

get_mean_vector<-function(pheY, pheT)
{
	t.count <- length(unique(pheT));
	if(t.count>=20)
	{
		t.min <- min(pheT, na.rm=T)
		t.max <- max(pheT, na.rm=T)
		pheT <- round((pheT - t.min)/(t.max-t.min)*20)/20*(t.max-t.min) + t.min;
	}

	t.all <- sort(unique(pheT));
	y.all <- c();
	for(t in t.all)
	{
		y.all <- c( y.all, mean(pheY[which(pheT==t)], na.rm=T) );
	}

	return(list(t=t.all, y=y.all));
}


##-----------------------------------------------------------
## Logistic curve
##
##    y = a/(1+b*exp(-r*t))
##
##-----------------------------------------------------------
log_get_cueve <- function(object, par, times, options=list())
{
	y <- par[1]/(1+par[2]*exp(-1 * par[3]*times) );

	return(y);
}

log_get_gradient <- function(object, par, times, options=list())
{
	d.a <- 1/(1+par[2]*exp(-1 * par[3]*times) );
	d.b <- (-1)* par[1] / ( (1+par[2]*exp(-1 * par[3]*times) )^2 ) * exp(-1 * par[3]*times);
	d.r <- (-1)* par[1] / ( (1+par[2]*exp(-1 * par[3]*times) )^2 ) * par[2] * exp(-1 * par[3]*times) * (-1*times);

	return(list(d.a, d.b, d.r));
}

log_get_param_info<-function(object, times, options=list())
{
	return(list(count=3, names=c("a", "b", "r"), formula="y = a/(1+b*exp(-r*t))" ));
}

log_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

log_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(18.18, 9.98, 0.99), c(17.08, 9.78, 0.97), c(15.95, 9.88, 0.98)	) );
}

log_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	mc <- get_mean_vector(pheY, pheT);
	mc$t <- mc$t[mc$y>0]
	mc$y <- mc$y[mc$y>0]
	m <- length(mc$t);
	if(m==0)
	{
		mc <- get_mean_vector(pheY, pheT);
		mc$y <- mc$y - min(mc$y)*1.01;
	}

	par <- c();
	ls.i <- ls.max <- Inf;
	a.rate <- mean(mc$y[-1]/mc$y[-m]);
	if(a.rate==1) a.rate <- 0.99
	for(i in 1:10)
	{
		par.a <- mc$y[m] * a.rate^i;
		par.r <- try( (log(par.a/mc$y[1]-1) - log(par.a/mc$y[m]-1))/(mc$t[1]-mc$t[m]));
		if(class(par.r)=="try-error" || is.infinite(par.r) )
			next;

		par.b <- (par.a / mc$y[m] -1)/exp(-par.r*mc$t[m]);

		y.ls <- sum(abs(mc$y - par.a/(1+par.b*exp(-par.r*mc$t)))^2, na.rm=T);

		if (y.ls < ls.max)
		{
			ls.i <- i;
			ls.max <- y.ls;
			par <- c(par.a, par.b, par.r);
		}
	}

	return(par*runif(length(par),0.95, 1.05) );
}

##-----------------------------------------------------------
## S4 class:
##           fg.curve.log
##-----------------------------------------------------------

setClass("fg.curve.log",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.log" ), log_get_cueve )

setMethod("get_gradient",  signature( object = "fg.curve.log" ), log_get_gradient )

setMethod("get_param_info",  signature( object = "fg.curve.log" ), log_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.log" ), log_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.log"), log_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.log"), log_est_init_param)


##-----------------------------------------------------------
## Double Logistic
##
##    y = a1/(1+b1*exp(-r1*t)) +  a2/(1+b2*exp(-r2*t))
##
##-----------------------------------------------------------

log2_show<-function(object)
{

}

log2_get_cueve <- function(object, par, times, options=list())
{
	y <- par[1]/(1+par[2]*exp(-par[3]*times) ) + par[4]/(1+par[5]*exp(-par[6]*times) )

	return(y);
}

log2_get_gradient <- function(object, par, times, options=list())
{
	d.a1 <- 1/(1+par[2]*exp(-1 * par[3]*times) );
	d.b1 <- (-1)* par[1] / ( (1+par[2]*exp(-1 * par[3]*times) )^2 ) * exp(-1 * par[3]*times);
	d.r1 <- (-1)* par[1] / ( (1+par[2]*exp(-1 * par[3]*times) )^2 ) * par[2] * exp(-1 * par[3]*times) * (-1*times);

	d.a2 <- 1/(1+par[5]*exp(-1 * par[6]*times) );
	d.b2 <- (-1)* par[4] / ( (1+par[5]*exp(-1 * par[6]*times) )^2 ) * exp(-1 * par[6]*times);
	d.r2 <- (-1)* par[4] / ( (1+par[5]*exp(-1 * par[6]*times) )^2 ) * par[5] * exp(-1 * par[6]*times) * (-1*times);

	return(list(d.a1, d.b1, d.r1, d.a2, d.b2, d.r2));
}


log2_get_param_info<-function(object, times, options=list())
{
	return(list(count=6,
	        names=c("a1", "b1", "r1","a2", "b2", "r2"),
	        formula="y = a1/(1+b1*exp(-r1*t)) +  a2/(1+b2*exp(-r2*t))" ));
}

log2_check_param<-function(object, times, options=list())
{
	return(TRUE);

}

log2_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(18.18, 9.98, 0.99), c(17.08, 9.78, 0.97), c(15.95, 9.88, 0.98)	) );
}

log2_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	par <- log_est_init_param(object, pheY, pheX, pheT, options)
	par[1] <- par[1]/2;
	return( c(par,par) );
}

##-----------------------------------------------------------
## S4 class:
##           fg.curve.log2
##-----------------------------------------------------------

setClass("fg.curve.log2",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.log2" ), log2_get_cueve )

setMethod("get_gradient",  signature( object = "fg.curve.log2" ), log2_get_gradient )

setMethod("get_param_info",  signature( object = "fg.curve.log2" ), log2_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.log2" ), log2_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.log2"), log2_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.log2"), log2_est_init_param)


##-----------------------------------------------------------
## ABRK:
##
##    y = a*(1+b*exp(-r*t))^(1/(1-k))
##
## Reference:<no>
##-----------------------------------------------------------

abrk_show<-function(object)
{

}

abrk_get_cueve <- function(object, par, times, options=list())
{
	y<- par[1]*( 1 + par[2]*exp(-par[3]*times) )^(1/(1-par[4]) );
	return(y);
}

abrk_get_gradient <- function(object, par, times, options=list())
{
	d.a <- ( 1 + par[2]*exp(-par[3]*times) )^(1/(1-par[4]) );
	d.b <- par[1] * (1/(1-par[4])) * ( 1 + par[2]*exp(-par[3]*times) )^( 1/(1-par[4]) -1 ) * exp( -par[3]*times );
	d.r <- par[1] * (1/(1-par[4])) * ( 1 + par[2]*exp(-par[3]*times) )^( 1/(1-par[4]) -1 ) * par[2]*exp(-par[3]*times) * (-1*times);
	d.k <- par[1] * (1 + par[2]*exp(-par[3]*times))^( 1/(1-par[4])) * log((1 + par[2]*exp(-par[3]*times))) * (1-par[4])^(-2);

	return(list(d.a, d.b, d.r, d.k));
}


abrk_get_param_info<-function(object, times, options=list())
{
	return(list(count=4, names=c("a","b","r","k"),
	       formula="y = a*(1+b*exp(-r*t))^(1/(1-k))" ) );
}

abrk_check_param<-function(object, times, options=list())
{
	return(TRUE);
}

abrk_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(18.18, 9.98, 0.99, 2.6 ), c(17.08, 9.78, 0.97, 2.5), c(15.95, 9.88, 0.98, 2.4)	) );
}

abrk_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	par <- log_est_init_param(object, pheY, pheX, pheT, options)
	par.k <- 2;

	return( c(par, par.k) * runif( 4, 0.95, 1.05) );
}

##-----------------------------------------------------------
## S4 class:
##           fg.curve.abrk
##-----------------------------------------------------------

setClass("fg.curve.abrk",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.abrk" ), abrk_get_cueve )

setMethod("get_gradient",  signature( object = "fg.curve.abrk" ), abrk_get_gradient )

setMethod("get_param_info",  signature( object = "fg.curve.abrk" ), abrk_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.abrk" ), abrk_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.abrk"), abrk_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.abrk"), abrk_est_init_param)


#-----------------------------------------------------------------
# Pharmacology Curve
#
#    y = E0 + Emax*t/(E50+t)
#
#-----------------------------------------------------------------

pc_show <- function(object)
{

}

pc_get_cueve <- function(object, par, times, options=list())
{
	return(  par[1] + par[3]*times/(par[2] + times) );
}

pc_get_gradient <- function(object, par, times, options=list())
{
	d.E0 <-  times/times;
	d.E50 <-  (-1) * par[3]*times/(par[2] + times)^2 ;
	d.Emax <-  times/(par[2] + times) ;

	return(list(d.E0, d.E50, d.Emax));
}

pc_get_param_info<-function(object, times, options=list())
{
	return(list(count=3, names=c("E0", "E50", "Emax"),
	       formula="y = E0 + Emax*t/(E50+t)" ));
}

pc_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

pc_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(10.9824, 15.909, 20.7768), c( 8.9824, 16.098, 20.7768), c(6.9507, 12.090,18.5737) ) );
}

pc_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	mc <- get_mean_vector(pheY, pheT);
	m <- length(mc$t);

	par <- c();
	ls.max <- Inf;
	a.rate <- mc$y[m]/mc$y[m-1];
	for(i in 1:10)
	{
		par.Emax <- mc$y[m]*a.rate^i;
		par.E50 <- mc$y[which(mc$t==median(mc$t))];
		par.E0 <- mc$y[m] - par.Emax*mc$t[m]/(par.E50 + mc$t[m] );

		y.ls <- sum(abs(mc$y - par.E0 - par.Emax*mc$t/(par.E50+mc$t))^2, na.rm=T);
		if (y.ls < ls.max)
		{
			ls.max <- y.ls;
			par <- c( par.E0, par.E50, par.Emax);
		}
	}


	return(par*runif(3, 0.95, 1.05))
}

##-----------------------------------------------------------
## S4 class:
##           fg.curve.pharmacology
##-----------------------------------------------------------

setClass("fg.curve.pharmacology",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.pharmacology" ), pc_get_cueve )

setMethod("get_gradient",  signature( object = "fg.curve.pharmacology" ), pc_get_gradient )

setMethod("get_param_info",  signature( object = "fg.curve.pharmacology" ), pc_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.pharmacology" ), pc_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.pharmacology"), pc_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.pharmacology"), pc_est_init_param)


#-----------------------------------------------------------------
# Exponential  Curve
#
#    y = a*exp(r*t)
#
#-----------------------------------------------------------------

exp_show <- function(object)
{
}

exp_get_cueve <- function(object, par, times, options=list())
{
	return(  par[[1]]*exp( par[[2]]*times ) );
}

exp_get_gradient <- function(object, par, times, options=list())
{
	d.a <- exp( par[[2]]*times );
	d.r <- par[[1]] * exp( par[[2]]*times ) * times;

	return(list(d.a, d.r));
}

exp_get_param_info<-function(object, times, options=list())
{
	return(list(count=2, names=c("a", "r"), formula="y = a*exp(r*t)"));
}

exp_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

exp_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c( 2, 0.0128), c( 1.8, 0.02), c(1.6, 0.024) ) );
}

exp_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	suppressWarnings( r <-  mean( (log( pheY[, 2]) - log(pheY[, 1]) )/(pheT[,2]-pheT[,1]), na.rm=T) )
	w0 <- mean( pheY[, 1]/exp(r*pheT[,1]), na.rm=T);

	return(c(w0,r ));
}

##-----------------------------------------------------------
## S4 class:
##           fg.curve.exponential
##-----------------------------------------------------------

setClass("fg.curve.exponential",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.exponential" ), exp_get_cueve )

setMethod("get_gradient",  signature( object = "fg.curve.exponential" ), exp_get_gradient )

setMethod("get_param_info",  signature( object = "fg.curve.exponential" ), exp_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.exponential" ), exp_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.exponential"), exp_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.exponential"), exp_est_init_param)



#-----------------------------------------------------------------
# Bi-exponential Curve
#
#    y = a1*exp(-r1*t) + a2*exp(-r2*t)
#
#-----------------------------------------------------------------

biexp_show <- function(object)
{

}

biexp_get_cueve <- function(object, par, times, options=list())
{
	return(  par[1]*exp(-par[2]*times) + par[3]*exp(-par[4]*times) );
}

biexp_get_gradient <- function(object, par, times, options=list())
{
	d.a1 <- exp(-par[2]*times);
	d.r1 <- par[1]*exp(-par[2]*times) *(-times) ;
	d.a2 <- exp(-par[4]*times) ;
	d.r2 <- par[3]*exp(-par[4]*times) *(-times) ;

	return(list(d.a1, d.r1, d.a2, d.r2));
}


biexp_get_param_info<-function(object, times, options=list())
{
	return(list(count=4, names=c("a1", "r1", "a2", "r2"),
	          formula="y = a1*exp(-r1*t) + a2*exp(-r2*t)"));
}

biexp_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

biexp_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(19.9824, 0.4699, 8.7768,  1.4699), c( 17.9824, 0.0699, 9.7768, 1.0699), c(15.9507, 0.1836, 10.5737, 1.8836) ) );
}

biexp_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	suppressWarnings( r <-  -1* mean( (log( pheY[, 2]) - log(pheY[, 1]))/(pheT[,2]-pheT[,1]), na.rm=T) )
	a_double <- mean(pheY[, 1]/exp(-1*r*pheT[,1]), na.rm=T);

	return(c(a_double/2, r, a_double/2, r))
}

##-----------------------------------------------------------
## S4 class:
##           fg.curve.bi.exponential
##-----------------------------------------------------------

setClass("fg.curve.bi.exponential",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.bi.exponential" ), biexp_get_cueve )

setMethod("get_gradient",  signature( object = "fg.curve.bi.exponential" ), biexp_get_gradient )

setMethod("get_param_info",  signature( object = "fg.curve.bi.exponential" ), biexp_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.bi.exponential" ), biexp_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.bi.exponential"), biexp_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.bi.exponential"), biexp_est_init_param)



##-----------------------------------------------------------------
## Power Curve
##
##         y = a*t^b
##
##-----------------------------------------------------------------
power_show <- function(object)
{

}

power_get_cueve <- function(object, par, times, options=list())
{
	return( par[1]*( times^par[2]) );
}

power_get_gradient <- function(object, par, times, options=list())
{
	return( list( times^par[2] , par[1]*( times^par[2])*log(times) )  );
}


power_get_param_info<-function(object, times, options=list())
{
	return(list(count=2, names=c("a","b"), formula="y = a*t^b"));
}

power_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

power_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(simu_a = 11.049, simu_b = 1.151), c(simu_a = 9.049, simu_b = 1.251), c(simu_a = 7.148,  simu_b = 1.359) ) );
}

power_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	suppressWarnings ( b <- mean( (log(pheY[, NCOL(pheY)])-log(pheY[, 1])) / (log(pheT[, NCOL(pheT)])-log(pheT[, 1])), na.rm=T) );
	suppressWarnings ( a <- exp( mean( log(pheY[, NCOL(pheY)]) - b*log(pheT[, NCOL(pheY)]), na.rm=T )) );
	return(c(a, b));
}


##-----------------------------------------------------------
## S4 class:
##           fg.curve.power
##-----------------------------------------------------------

setClass("fg.curve.power",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.power" ), power_get_cueve )

setMethod("get_gradient",  signature( object = "fg.curve.power" ), power_get_gradient )

setMethod("get_param_info",  signature( object = "fg.curve.power" ), power_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.power" ), power_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.power"), power_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.power"), power_est_init_param)


##-----------------------------------------------------------------
## Legendre:  Legendre Polynomial(2nd-order)
##
##        y = u0 + u1 *t + u2*1/2*(3*t^2-1)
##
##-----------------------------------------------------------------

Legendre2_show <- function(object)
{
}

Legendre2_get_cueve <- function(object, par, times, options=list())
{
	ti <- -1 + 2*(times - options$min.time)/( options$max.time - options$min.time );
	return( par[1] + ti*par[2] + 0.5*(3*ti*ti-1)* par[3]  );
}

Legendre2_get_gradient <- function(object, par, times, options=list())
{
	ti <- -1 + 2*(times - options$min.time)/( options$max.time - options$min.time );
	return( list( array(1, dim=dim(ti)), ti, 0.5*(3*ti*ti-1) )  );
}


Legendre2_get_param_info<-function(object, times, options=list())
{
	return(list(count=3, names=c("u0","u1","u2"), formula="y = u0 + u1*t + u2*1/2*(3*t^2-1)"));
}

Legendre2_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

Legendre2_get_simu_param<-function(object, times, options=list())
{
	QQ2 = c(simu_u0 = 11.049, simu_u1 = 1.551, simu_u2 = -8.019, simu_u3 = 3.151, simu_u4 = 0.652, simu_u5 = -0.597, simu_u6 = 0.821);
	Qq1 = c( simu_u0 = 9.049, simu_u1 = 1.151, simu_u2 = -6.019, simu_u3 = 2.651, simu_u4 = 0.652, simu_u5 = -0.797, simu_u6 = 0.621);
	qq0 = c( simu_u0 = 7.148, simu_u1 = 1.379, simu_u2 = -4.489, simu_u3 = 2.004, simu_u4 = 0.662, simu_u5 = -0.836, simu_u6 = 0.432)

	return( rbind( QQ2, Qq1, qq0)[, c(1:3),drop=F ] );
}

Legendre2_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	ti <- -1 + 2*(pheT-options$min.time)/( options$max.time - options$min.time );
    y1 <- pheY[,1];
    y2 <- pheY[,2];
    y3 <- pheY[,3];

	u2 <- mean( ( (y1-y2)/(ti[,1]-ti[,2]) - (y1-y3)/(ti[,1]-ti[,3]) ) / ( (ti[,1]+ti[,2])-(ti[,1]+ti[,3]) ) / 1.5, na.rm = T);
	u1 <- mean( (y1-y2) /(ti[,1]-ti[,2]) - u2 * 1.5 * (ti[,1]+ti[,2]), na.rm=T);
	u0 <- mean( y3 + u2*0.5 - u1*ti[,3] - u2 * 1.5 * ti[,3]^2, na.rm=T);

	return(c(u0, u1, u2));
}


##-----------------------------------------------------------
## S4 class:
##           fg.curve.legendre2
##-----------------------------------------------------------

setClass("fg.curve.legendre2",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.legendre2" ), Legendre2_get_cueve )

setMethod("get_gradient",  signature( object = "fg.curve.legendre2" ), Legendre2_get_gradient )

setMethod("get_param_info",  signature( object = "fg.curve.legendre2" ), Legendre2_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.legendre2" ), Legendre2_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.legendre2"), Legendre2_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.legendre2"), Legendre2_est_init_param)



##-----------------------------------------------------------------
## Legendre:  Legendre Polynomial(3rd-order)
##
##        y = u0 + u1 *t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t)
##
##-----------------------------------------------------------------

Legendre3_show <- function(object)
{
}

Legendre3_get_cueve <- function(object, par, times, options=list())
{
	ti <- -1 + 2*(times-options$min.time)/( options$max.time - options$min.time );
	return( par[1] + ti*par[2] + 0.5*(3*ti*ti-1)* par[3] +  0.5*(5*ti^3-3*ti)*par[4] );
}

Legendre3_get_gradient <- function(object, par, times, options=list())
{
	ti <- -1 + 2*(times-options$min.time)/( options$max.time - options$min.time );
	return( list( array(1, dim=dim(ti)), ti, 0.5*(3*ti*ti-1), 0.5*(5*ti^3-3*ti) ) );
}


Legendre3_get_param_info<-function(object, times, options=list())
{
	return(list(count=4, names=c("u0","u1","u2","u3"), formula="y = u0 + u1*t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) "));
}

Legendre3_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

Legendre3_get_simu_param<-function(object, times, options=list())
{
	QQ2 = c(simu_u0 = 11.049, simu_u1 = 1.551, simu_u2 = -8.019, simu_u3 = 3.151, simu_u4 = 0.652, simu_u5 = -0.597, simu_u6 = 0.821);
	Qq1 = c( simu_u0 = 9.049, simu_u1 = 1.151, simu_u2 = -6.019, simu_u3 = 2.651, simu_u4 = 0.652, simu_u5 = -0.797, simu_u6 = 0.621);
	qq0 = c( simu_u0 = 7.148, simu_u1 = 1.379, simu_u2 = -4.489, simu_u3 = 2.004, simu_u4 = 0.662, simu_u5 = -0.836, simu_u6 = 0.432)

	return( rbind( QQ2, Qq1, qq0)[, c(1:4),drop=F] );
}

Legendre3_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	u3 <- Legendre2_est_init_param(object, pheY, pheX, pheT, options);
	return(c(u3, 0.0001));
}


##-----------------------------------------------------------
## S4 class:
##           fg.curve.legendre3
##-----------------------------------------------------------

setClass("fg.curve.legendre3",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.legendre3" ), Legendre3_get_cueve )

setMethod("get_gradient",  signature( object = "fg.curve.legendre3" ), Legendre3_get_gradient )

setMethod("get_param_info",  signature( object = "fg.curve.legendre3" ), Legendre3_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.legendre3" ), Legendre3_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.legendre3"), Legendre3_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.legendre3"), Legendre3_est_init_param)



##-----------------------------------------------------------------
## Legendre:  Legendre Polynomial(4th-order)
##
##        y = u0 + u1 *t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) + u4*1/8*(35*ti^4-30*ti^2+3)
##
##-----------------------------------------------------------------

Legendre4_show <- function(object)
{
}

Legendre4_get_cueve <- function(object, par, times, options=list())
{
	ti <- -1 + 2*(times-options$min.time)/( options$max.time - options$min.time );
	return( par[1] + ti*par[2] + 0.5*(3*ti*ti-1)* par[3] +  0.5*(5*ti^3-3*ti)*par[4] +  0.125*(35*ti^4-30*ti^2+3)* par[5] );
}

Legendre4_get_gradient <- function(object, par, times, options=list())
{
	ti <- -1 + 2*(times-options$min.time)/( options$max.time - options$min.time );
	return( list( array(1, dim=dim(ti)), ti, 0.5*(3*ti*ti-1), 0.5*(5*ti^3-3*ti),  0.125*(35*ti^4-30*ti^2+3) )  );
}


Legendre4_get_param_info<-function(object, times, options=list())
{
	return(list(count=5, names=c("u0","u1","u2","u3", "u4"), formula="y = u0 + u1*t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) + u4*1/8*(35*ti^4-30*ti^2+3)"));
}

Legendre4_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

Legendre4_get_simu_param<-function(object, times, options=list())
{
	QQ2 = c(simu_u0 = 11.049, simu_u1 = 1.551, simu_u2 = -8.019, simu_u3 = 3.151, simu_u4 = 0.652, simu_u5 = -0.597, simu_u6 = 0.821);
	Qq1 = c( simu_u0 = 9.049, simu_u1 = 1.151, simu_u2 = -6.019, simu_u3 = 2.651, simu_u4 = 0.652, simu_u5 = -0.797, simu_u6 = 0.621);
	qq0 = c( simu_u0 = 7.148, simu_u1 = 1.379, simu_u2 = -4.489, simu_u3 = 2.004, simu_u4 = 0.662, simu_u5 = -0.836, simu_u6 = 0.432)

	return( rbind( QQ2, Qq1, qq0)[, c(1:5),drop=F] );
}

Legendre4_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	u3 <- Legendre2_est_init_param(object, pheY, pheX, pheT, options);
	return(c(u3, 0.0001, 0.0001));
}


##-----------------------------------------------------------
## S4 class:
##           fg.curve.legendre4
##-----------------------------------------------------------

setClass("fg.curve.legendre4",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.legendre4" ), Legendre4_get_cueve )

setMethod("get_gradient",  signature( object = "fg.curve.legendre4" ), Legendre4_get_gradient )

setMethod("get_param_info",  signature( object = "fg.curve.legendre4" ), Legendre4_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.legendre4" ), Legendre4_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.legendre4"), Legendre4_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.legendre4"), Legendre4_est_init_param)


##-----------------------------------------------------------
## Chapman-Richard
##
##    y = a*(1-exp(-rt))^b
##
##-----------------------------------------------------------

cr_show <- function(object)
{
}

cr_get_cueve <- function(object, par, times, options=list())
{
	return ( par[1]*(1 - exp(-par[3]*times))^par[2] );
}

cr_get_gradient <- function(object, par, times, options=list())
{
	d.a <- (1 - exp(-par[3]*times))^par[2] ;
	d.b <- par[1] * (1 - exp(-par[3]*times))^par[2] * log( 1 - exp(-par[3]*times) );
	d.r <- par[1] * (1 - exp(-par[3]*times))^(par[2]-1) * par[2]* exp(-par[3]*times) *(-1*times);

	return ( list(d.a, d.b, d.r) );
}


cr_get_param_info<-function(object, times, options=list())
{
	return(list(count=3, names=c("a","b","r"), formula="y = a*(1-exp(-rt))^b"));
}

cr_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

cr_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(21.98, 0.47, 9.78), c(19.98, 0.47, 8.77), c(15.95, 0.48, 7.58)	) );
}

cr_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	a <- max(pheY, na.rm=T);
	b <- 1
	r <- mean( log(1-pheY/a)/pheT/(-1), trim=0.2, na.rm=T)
	if( is.infinite(r)) r<-0;

	return(c(a, b, r));
}


##-----------------------------------------------------------
## S4 class:
##           fg.curve.ChapmanRichard
##-----------------------------------------------------------

setClass("fg.curve.ChapmanRichard",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.ChapmanRichard" ), cr_get_cueve )

setMethod("get_gradient",  signature( object = "fg.curve.ChapmanRichard" ), cr_get_gradient )

setMethod("get_param_info",  signature( object = "fg.curve.ChapmanRichard" ), cr_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.ChapmanRichard" ), cr_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.ChapmanRichard"), cr_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.ChapmanRichard"), cr_est_init_param)
