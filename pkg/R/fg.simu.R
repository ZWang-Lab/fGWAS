#---------------------
# fg_simulate
#
#   return fg.obj
#
#	# str(fg.obj)
#	# $ obj.gen : 'fgwas.gen.obj':
#	# $ obj.phe : 'fgwas.phe.obj':
#	# $ error   : logi FALSE
#---------------------

fg_simulate<-function( curve.type, covariance.type, n.obs, n.snp, time.points, par0=NULL, par1=NULL, par2=NULL, par.covar=NULL, par.X=NULL,
		phe.missing=0.03, snp.missing=0.03, sig.pos=NULL, plink.format=FALSE, file.prefix=NULL )
{
	## check 'curve.type'
	if ( missing(curve.type) || is.null(curve.type) )
		stop("No curve type is specified in the paramater 'curve.type'.")
	else
		if( is.character(curve.type) && length(curve.type) > 1 )
			stop( paste( "'curve.type' should be a string indicating curve type." ) );

	## check 'covariance.type'
	if ( missing(covariance.type) || is.null(covariance.type) )
		stop("No covariance type is specified in the paramater 'covariance.type'.")
	else
		if( is.character(covariance.type) && length(covariance.type) > 1 )
			stop( paste( "'covariance.type' should be a string indicating covariance type." ) );

	## check 'n.obs'
	if ( missing(n.obs) || is.null(n.obs) )
		stop("Individual count is not specified in the paramater 'n.obs'.")
	else
		if ( !(is.numeric( n.obs ) && length(n.obs)) )
			stop("Not integer in parameter 'n.obs'.");

	## check 'n.snp'
	if ( missing(n.snp) || is.null(n.snp) )
		stop("SNP count is not specified in the paramater 'n.snp'.")
	else
		if ( !(is.numeric( n.snp ) && length(n.snp)) )
			stop("Not integer in parameter 'n.snp'.");

	## check 'time.points'
	if ( missing(time.points) || is.null(time.points) )
		stop("Time points are not specified in the paramater 'time.points'.")
	else
	{
		if ( !(is.numeric( time.points ) && length(time.points)) )
			stop("Not integer in parameter 'time.points'.");

		if(length(time.points)==1)
			time.points <- 1:time.points
	}

	## check 'par.X'
	if( !missing(par.X)  && !is.null(par.X) )
		if( !all(is.numeric(par.X)) )
			stop( paste( "Covariate paramater should be numeric values." ) );

	## check 'file.prefix'
	if ( !missing(file.prefix) && !is.null(file.prefix))
		if( !is.character(file.prefix) || length(file.prefix)>1 )
			stop( paste( "'file.prefix' should be a string." ) );

	if(plink.format && (missing(file.prefix) || is.null(file.prefix)))
		stop("'file.prefix' is NULL")

	## check 'curve.type'
	if( class(curve.type)=="fgwas.curve")
		fg_curve <- curve.type
	else
		fg_curve <- fg.getCurve( curve.type );

	simu_parm <- get_simu_param( fg_curve, time.points );
	simu_len <- NCOL(simu_parm);

	## check 'par0'
	if( missing(par0) || is.null(par0))
		par0 <- simu_parm[1,]
	else
	{
		if( !( all(is.numeric(par0)) && length(par0)==simu_len) )
			stop( paste( "Curve paramater should be", simu_len, "numeric values." ) );
	}

	## check 'par1'
	if( missing(par1) || is.null(par1))
		par1 <- simu_parm[2,]
	else
	{
		if( !( all(is.numeric(par1)) && length(par1)==simu_len) )
			stop( paste( "Curve paramater should be", simu_len, "numeric values." ) );
	}

	## check 'par2'
	if( missing(par2) || is.null(par2))
		par2 <- simu_parm[3,]
	else
	{
		if( !( all(is.numeric(par2)) && length(par2)==simu_len) )
			stop( paste( "Curve paramater should be", simu_len, "numeric values." ) );
	}

	## check 'covariance.type'
	if(class(covariance.type)=="fgwas.covar")
		fg_covar <- covariance.type
	else
		fg_covar <- fg.getCovariance( covariance.type );

	simu_covar <- get_simu_param(fg_covar, time.points);
	simu_len <- NROW(simu_covar);

	## check 'par.covar'
	if(missing(par.covar) || is.null(par.covar) )
		par.covar <- simu_covar
	else
	{
		if( !( all(is.numeric(par.covar)) && length(par.covar)==simu_len) )
			stop( paste( "Covariance paramater should be", simu_len, "numeric numbers." ) );
	}

	## check 'sig.pos'
	if ( missing(sig.pos) || is.null(sig.pos))
	{
		sig.pos <- round( runif(1, n.snp*0.25, n.snp*0.75) );
		cat(" * A significant SNP is randomly specified to location(", sig.pos, ")\n" );
	}
	else
	{
		if ( !(is.numeric( sig.pos ) && length(sig.pos)) )
			stop("Not integer in parameter 'n.snp'.");
	}

	fg.obj <- proc_dat_simu( n.obs, n.snp, par.X, par0, par1, par2, par.covar, fg_curve, fg_covar, time.points, sig.pos, snp.missing, phe.missing );

	if(!is.null(file.prefix))
	{
		fg.obj$obj.phe$file.pheX <- paste( file.prefix, ".pheX.csv", sep="" );
		fg.obj$obj.phe$file.pheY <- paste( file.prefix, ".pheY.csv", sep="" );
		fg.obj$obj.phe$file.pheT <- paste( file.prefix, ".pheT.csv", sep="" );

		write.csv(data.frame(ID = fg.obj$obj.phe$ids, fg.obj$obj.phe$pheY), file = fg.obj$obj.phe$file.pheY, quote=F, row.names=F );
		write.csv(data.frame(ID = fg.obj$obj.phe$ids, fg.obj$obj.phe$pheT), file = fg.obj$obj.phe$file.pheT, quote=F, row.names=F );
		if(!is.null(fg.obj$obj.phe$pheX))
			write.csv(data.frame(ID = fg.obj$obj.phe$ids, fg.obj$obj.phe$pheX), file = fg.obj$obj.phe$file.pheX, quote=F, row.names=F )
		else
			fg.obj$obj.phe$file.pheX=NULL;

		if ( !plink.format )
		{
			file.gen.dat  <- paste(file.prefix, ".geno.tab", sep="");
			tb.gen <- data.frame(fg.obj$obj.gen$reader$get_snpinfo( NULL ),
								 fg.obj$obj.gen$reader$get_snpmat( NULL, impute=F, allel=T )$snpmat );
			colnames(tb.gen) <- c(colnames(tb.gen), rownames(fg.obj$obj.phe$file.pheX) );
			write.table( tb.gen, file=file.gen.dat, quote=F, row.names=F, col.names=T, sep="\t" );

			fg.obj$obj.gen$files = list(file.gen.dat);
		}
		else
		{
			snp.mat <- fg.obj$obj.gen$reader$get_snpmat( NULL, impute=F, allel=F)$snpmat;
			snp.info <- fg.obj$obj.gen$reader$get_snpinfo(NULL );

			r <- convert_simpe_to_plink( data.frame(snp.info[,c(2,1)], 0, snp.info[,c(3:5)]),  snp.mat, paste(file.prefix, ".geno", sep="") );

			fg.obj$obj.gen$files = list(
				file.plink.bed = r$file.plink.bed,
   		    	file.plink.bim = r$file.plink.bim,
   		    	file.plink.fam = r$file.plink.fam);
		}
	}

	return(fg.obj);
}

proc_dat_simu<-function( n.obs, n.snp, par.X, par0, par1, par2, par.covar, f.curve, f.covar, times, sig.idx, snp.missing, phe.missing )
{
	obj.gen <- list();
	obj.phe <- list(options=list(), params=list(), intercept=F);

	#generate SNPs
	tmp.matrix <- proc_simu_geno( n.obs, n.snp, sig.idx, prob.miss = snp.missing )
	tmp.snpinfo <- data.frame(SNPNAME=paste("SNP", 1:n.snp, sep="_"), CHR=1, POS=1:n.snp, RefBase="A", AltBase="B");
	obj.gen$n.snp          = n.snp;
	obj.gen$n.ind.total    = n.obs;
    obj.gen$n.ind.used     = n.obs;

	obj.gen$reader <- new("fg.dm.simple",
			type           = "SIMPLE",
			description    = "Simple geno table",
			n.snp          = n.snp,
			n.ind.total    = n.obs,
			n.ind.used     = n.obs,
			ids.used       = colnames(tmp.matrix),
			file.simple.snp= "",
			rawdata        = data.frame(tmp.snpinfo , tmp.matrix),
			snpdata        = list(snp.info=tmp.snpinfo , snp.mat=tmp.matrix) );

	pheX <- NULL;
	if(length(par.X)>0)
	{
		pheX <- matrix(0, nrow=n.obs, ncol=0 );
		rownames(pheX) <- paste("N", 1:n.obs, sep="_");

		for(i in 1:length(par.X) )
		{
			if(i==1)
				pheX <- cbind( pheX,round( runif(n.obs, 1, 2) ) )
			else
				pheX <- cbind( pheX, runif(n.obs, -1, 1) );
		}
		colnames(pheX) <- paste("X", 1:length(par.X), sep="_");
	}


	pheY <- array( 0, dim=c(n.obs, length(times)));

	colnames( pheY ) <- paste("Y", times, sep="_");
	rownames( pheY ) <- paste("N", 1:n.obs, sep="_");

	#generate traits
	options <- list( max.time=max(times, na.rm=T), min.time=min(times, na.rm=T) );
	sim.mu   <-  get_curve( f.curve, par0, times, options=options );
	sim.mu   <-  rbind(sim.mu, get_curve( f.curve, par1, times, options=options ) );
	sim.mu   <-  rbind(sim.mu, get_curve( f.curve, par2, times, options=options ) );

	get_gencode <- function(d.gen, d.snpinfo)
	{
		d.g <- array( 9, n.obs );
		d.gen2 <- as.character(unlist(d.gen));
		snpB <- as.character(unlist(d.snpinfo[5]));
		snpA <- as.character(unlist(d.snpinfo[4]));

		QQ2<- paste(snpB, snpB, sep="");
		qq0<- paste(snpA, snpA, sep="");
		Qq1<- c(paste(snpA, snpB, sep=""), paste( snpB, snpA, sep="") ) ;

		d.g[which(d.gen2==QQ2)]<-2;
		d.g[which(d.gen2==qq0)]<-0;
		d.g[which(d.gen2==Qq1[1])]<-1;
		d.g[which(d.gen2==Qq1[2])]<-1;

		return(d.g);
	}

	d.g <- get_gencode(tmp.matrix[sig.idx,],  tmp.snpinfo[sig.idx,])
	d.g1 <- get_gencode(tmp.matrix[1,],  tmp.snpinfo[1,])
	d.g2 <- get_gencode(tmp.matrix[2,],  tmp.snpinfo[2,])

if(.RR("debug")) cat("COR(g, g1)=", cor(d.g, d.g1), "\n");
if(.RR("debug")) cat("COR(g, g2)=", cor(d.g, d.g2), "\n");

	sim.covar <- get_matrix( f.covar, par.covar,times );
	for (i in 1:n.obs)
	{
		 if (d.g[i]==9) d.g[i] <- round(runif(1, 0, 2));

		 pheY[i, ] <- rmvnorm(1, sim.mu[ d.g[i] + 1, ], sim.covar );
		 if( !is.null(pheX) )
		 	pheY[i, ]  <- pheY[i, ] + sum( pheX[i,] * par.X );
	}

	obj.phe$obj.curve <- f.curve;
	obj.phe$obj.covar <- f.covar;
	obj.phe$ids       <- rownames( pheY );
	obj.phe$pheY      <- pheY;
	obj.phe$pheX      <- pheX;
	obj.phe$pheT      <- matrix( rep(times, n.obs), nrow=n.obs, byrow=T)

	colnames(obj.phe$pheT) <- paste("T", times, sep="_");
	rownames(obj.phe$pheT) <- rownames( pheY );

	class(obj.gen) <- "fgwas.gen.obj"
	class(obj.phe) <- "fgwas.phe.obj"

	fg.obj <- list(obj.gen=obj.gen, obj.phe=obj.phe, error=F);

	cat(" Data simulation is done![Sig=", sig.idx, "]\n");

	return(fg.obj);
}

#--------------------------------------------------------------
# private: fin.generate_bc_marker;
#
# genarate N Backcross Markers from marker disttance (cM): dist.
#
# input
#     dist : the vector for the marker distances
#   samp_N : the sample size
#--------------------------------------------------------------
fin.generate_bc_marker<-function( n.obs, dist )
{
	if (dist[1] != 0)
		cm=c(0, dist)/100
	else
		cm=dist/100;

	n  <- length(cm);
	rs <- 1/2*( exp(2*cm) - exp(-2*cm) ) / (exp(2*cm)+exp(-2*cm));
	mk <- array( 0, dim=c( n.obs, n ) );

	for (j in 1:n.obs)
		mk[j,1] <- ( runif(1)>0.5 );

	for (i in 2:n)
		for ( j in 1:n.obs )
		{
			if (mk[j,i-1]==1)
				mk[j,i] <- ( runif(1)>rs[i] )
			else
				mk[j,i] <- ( runif(1)<rs[i] );
    	}

	return(mk);
}

proc_simu_geno<-function( n.obs, n.snp, sig.idx, prob.miss=0.03 )
{
	dist <- cumsum( runif(n.snp, 0.05, 0.12 ) );
	snp1 <- t( fin.generate_bc_marker( n.obs,  dist) );
	snp2 <- t( fin.generate_bc_marker( n.obs,  dist) );
	for(i in 1:n.snp)
	{
		nmiss <- runif( 1, 0, n.obs*prob.miss );
		if (nmiss>=1)
			snp1[ i, sample(n.obs)[1:round(nmiss)] ] <- 9;
	}

	for(i in 1:n.snp)
	{
		nmiss <- runif( 1, 0, n.obs*prob.miss );
		if (nmiss>=1)
			snp2[ i, sample(n.obs)[1:round(nmiss)] ] <- 9;
	}

	cors <- c();
	snpx <- snp1+snp2;
	for(i in 1:n.snp)
		cors <- c(cors, cor(snpx[i,], snpx[sig.idx, ]) );

	cor.high <- which( abs(cors)>0.75 );
	if ( length( cor.high )>=2)
	{
		for( i in 1:length(cor.high) )
			if( cor.high[i]!=sig.idx)
			{
				snp1[ cor.high[i], ] <- snp1[ cor.high[i] -1, ];
				snp2[ cor.high[i], ] <- snp2[ cor.high[i] -1, ];
			}
	}

	snp1.s <- c(as.character(snp1))
	snp2.s <- c(as.character(snp2))

	if( length( which(snp1.s=="0") ) >0 ) snp1.s[ which(snp1.s=="0") ] <- "A";
	if( length( which(snp1.s=="1") ) >0 ) snp1.s[ which(snp1.s=="1") ] <- "T";
	if( length( which(snp1.s=="9") ) >0 ) snp1.s[ which(snp1.s=="9") ] <- ".";

	if( length( which(snp2.s=="0") ) >0 ) snp2.s[ which(snp2.s=="0") ] <- "A";
	if( length( which(snp2.s=="1") ) >0 ) snp2.s[ which(snp2.s=="1") ] <- "T";
	if( length( which(snp2.s=="9") ) >0 ) snp2.s[ which(snp2.s=="9") ] <- ".";

	snp.s <- paste(snp1.s, snp2.s, sep="");
	gen<-array(snp.s, dim=c(n.snp, n.obs))
	colnames(gen)<-paste("N",1:n.obs, sep="_");
	rownames(gen)<-paste("SNP",1:n.snp, sep="_");

	return(gen);
}



#--------------------------------------------------------------
# private: fin.generate_bc_marker;
#
# genarate N Backcross Markers from marker disttance (cM): dist.
#
# input
#     dist : the vector for the marker distances
#   samp_N : the sample size
#--------------------------------------------------------------
fin.generate_bc_marker<-function( n.obs, dist )
{
	if (dist[1] != 0)
		cm=c(0, dist)/100
	else
		cm=dist/100;

	n  <- length(cm);
	rs <- 1/2*( exp(2*cm) - exp(-2*cm) ) / (exp(2*cm)+exp(-2*cm));
	mk <- array( 0, dim=c( n.obs, n ) );

	for (j in 1:n.obs)
		mk[j,1] <- ( runif(1)>0.5 );

	for (i in 2:n)
		for ( j in 1:n.obs )
		{
			if (mk[j,i-1]==1)
				mk[j,i] <- ( runif(1)>rs[i] )
			else
				mk[j,i] <- ( runif(1)<rs[i] );
    	}

	return(mk);
}

proc_simu_geno<-function( n.obs, n.snp, sig.idx, prob.miss=0.02 )
{
	dist <- cumsum( runif(n.snp, 0.05, 0.12 ) );
	snp1 <- t( fin.generate_bc_marker( n.obs,  dist) );
	snp2 <- t( fin.generate_bc_marker( n.obs,  dist) );
	for(i in 1:n.snp)
	{
		nmiss <- runif( 1, 0, n.obs*prob.miss );
		if (nmiss>=1)
			snp1[ i, sample(n.obs)[1:round(nmiss)] ] <- 9;
	}

	for(i in 1:n.snp)
	{
		nmiss <- runif( 1, 0, n.obs*prob.miss );
		if (nmiss>=1)
			snp2[ i, sample(n.obs)[1:round(nmiss)] ] <- 9;
	}

	cors <- c();
	snpx <- snp1+snp2;
	for(i in 1:n.snp)
		cors <- c(cors, cor(snpx[i,], snpx[sig.idx, ]) );

	cor.high <- which( abs(cors)>0.75 );
	if ( length( cor.high )>=2)
	{
		for( i in 1:length(cor.high) )
			if( cor.high[i]!=sig.idx)
			{
				snp1[ cor.high[i], ] <- snp1[ cor.high[i] -1, ];
				snp2[ cor.high[i], ] <- snp2[ cor.high[i] -1, ];
			}
	}

	snp1.s <- c(as.character(snp1))
	snp2.s <- c(as.character(snp2))

	if( length( which(snp1.s=="0") ) >0 ) snp1.s[ which(snp1.s=="0") ] <- "A";
	if( length( which(snp1.s=="1") ) >0 ) snp1.s[ which(snp1.s=="1") ] <- "B";
	if( length( which(snp1.s=="9") ) >0 ) snp1.s[ which(snp1.s=="9") ] <- ".";

	if( length( which(snp2.s=="0") ) >0 ) snp2.s[ which(snp2.s=="0") ] <- "A";
	if( length( which(snp2.s=="1") ) >0 ) snp2.s[ which(snp2.s=="1") ] <- "B";
	if( length( which(snp2.s=="9") ) >0 ) snp2.s[ which(snp2.s=="9") ] <- ".";

	snp.s <- paste(snp1.s, snp2.s, sep="");
	gen<-array(snp.s, dim=c(n.snp, n.obs))
	colnames(gen)<-paste("N",1:n.obs, sep="_");
	rownames(gen)<-paste("SNP",1:n.snp, sep="_");

	return(gen);
}

