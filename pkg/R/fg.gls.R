snpscan_nocurve<-function(obj.gen, obj.phe, snp.idx=NULL, options=list() )
{
	default_options <- list( ncores=1, method="GLS", verbose=TRUE);
	default_options[names(options)] <- options;
	options <- default_options;

	n.ind <- obj.phe$n.ind;
	n.snp <- obj.gen$n.snp;

	if(options$verbose)
	{
		cat("* Method =", options$method, "\n" );
		cat("* SNP Count =", obj.gen$n.snp, "\n" );
		cat("* Sample Count =", obj.phe$n.ind, "\n" );
		cat("* Parallel CPU Count =", options$ncores, "\n" );
		cat("* Verbose =", options$verbose, "\n" );
	}

	intercept <- ifelse ( !is.null(obj.phe$intercept), obj.phe$intercept, FALSE);

	if(is.null(snp.idx)) snp.idx <- 1:obj.gen$n.snp
	snp.len <- length(snp.idx);
	snp.sect0 <- seq( 1, snp.len, 10000 );
	snp.sect1 <- c( snp.sect0[-1]-1, snp.len );

	r.gls <- c();
	snp.mat <- c();

	r.fgwas <- c();
	for(i in 1:length(snp.sect0))
	{
		snp.idx.sub <- snp.idx[ snp.sect0[i]:snp.sect1[i] ];
		snp.mat  <- obj.gen$reader$get_snpmat( snp.idx.sub, impute=F );
		snp.info <- obj.gen$reader$get_snpinfo( snp.idx.sub);

		if(options$verbose)	cat("  Calculated SNP Range =", snp.sect0[i], snp.sect1[i], "\n" );

		snp.matx <- data.frame(snp.info[,c(1:3)], t(snp.mat$snpmat));
		rownames(snp.matx) <- colnames(snp.mat$snpmat)

		r.fgwas0 <- list();
		if( NCOL(obj.phe$pheY) == 1)
			r.fgwas0 <- bls.fgwas( i, snp.idx.sub, snp.matx, obj.phe$pheY, obj.phe$pheX, intercept=intercept, options$ncores )
		else
			r.fgwas0 <- gls.fgwas( i, snp.idx.sub, snp.matx, obj.phe$pheY, obj.phe$pheT, obj.phe$pheX, intercept=intercept, options$ncores )

		if( r.fgwas0$error )
			return(r.fgwas0);

		#adjust SNP.IDX field.
		#r.fgwas0$r[,1] <- r.fgwas0$r[,1] + snp.sect0[i]-1;

		r.fgwas <- rbind( r.fgwas, r.fgwas0$r);
	}

	return(list(error=F, result = r.fgwas) );
}

gls.fgwas <- function( sect.idx, snp.idx, snp.mat, pheY, pheZ, pheX=NULL, intercept=FALSE, ncores=1)
{
	if ( NCOL(pheY) != NCOL(pheZ) )
	{
		cat("! The Z columns are not matched with Y columns.\n");
		return(list(error=T, err.info="The Z columns are not matched with Y columns."))
	}

	covar.names <- c();
	if(!is.null(pheX)) covar.names <- colnames(pheX);

	phe.gls.mat <- array(0,dim=c(0, 1 + 2 + ifelse(is.null(pheX), 0, NCOL(pheX)) ))
	for(i in 1:NCOL(pheY))
	{
		phe.tmp <- data.frame( ID = 1:NROW(pheY), Y=pheY[,i,drop=F], Z=pheZ[,i,drop=F])
		if (!is.null(pheX)) phe.tmp <- data.frame(phe.tmp, pheX);

		phe.tmp.col <- colnames(phe.tmp);
		phe.tmp.col[ c(1:3) ] <- c("ID", "Y", "Z");
		colnames(phe.tmp) <- phe.tmp.col;

		phe.tmp.missing <- which( is.na(phe.tmp$Z) | is.na(phe.tmp$Y) );
		if(length(phe.tmp.missing)>0)
			phe.tmp <- phe.tmp[-phe.tmp.missing, , drop=F ];
		phe.gls.mat <- rbind( phe.gls.mat, phe.tmp );
	}

	reg.str0 <- ifelse( length(covar.names)>0, paste("Y ~ ", ifelse(intercept, "1+", ""), paste(covar.names,collapse= "+")),  "Y ~ 1" );
	reg.str1 <- ifelse( length(covar.names)>0, paste("Y ~ ", ifelse(intercept, "1+", ""), paste(covar.names,collapse= "+"), "+ as.factor(SNP)") ,  "Y ~ 1 + as.factor(SNP)" );

	if( sect.idx ==1)
	{
		cat("* H0 =", as.character(reg.str0), "\n" );
		cat("* H1 =", as.character(reg.str1), "\n" );
	}

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

		r.gls <- data.frame( INDEX=numeric(), NAME=character(), CHR=numeric(), POS=numeric(), NMISS=numeric(), MAF=numeric(), LR=numeric(), pv=numeric(),stringsAsFactors=F);

		reg.f0 <- as.formula(reg.str0);
		reg.f1 <- as.formula(reg.str1);

		for(i in c(range.fr:range.to) )
		{
			snp <- unlist( snp.mat[i, c(4:NCOL(snp.mat)) ] );
			snp.miss <- which( is.na(snp) | snp<0 | snp>2 )

			r.gls[i - range.fr + 1, "NAME"] <- as.character(snp.mat[i,1]);

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
				r.gls[i-range.fr+1, -c(2)] <- c( snp.idx[i], unlist(snp.mat[i,2:3]), length(snp.miss), maf,  NA, pv.max );
				next;
			}

			r <- try( anova( r0,r1 ) );
			if(any(class(r)=="try-error"))
			{
				r.gls[i-range.fr+1, -c(2)] <- c( snp.idx[i], unlist(snp.mat[i,2:3]), length(snp.miss), maf, NA, pv.max );
				next;
			}

			r.gls[i-range.fr+1, -c(2)] <- c( snp.idx[i], unlist(snp.mat[i,2:3]), length(snp.miss), maf, r[2,8], r[2,9]);
		}

		return(r.gls);
	}


	r.gls <- c();
	if( ncores > 1 && require("snowfall") )
	{
		if( .RR("debug") )
			cat("Starting parallel computing, snowfall/snow......\n");
		snowfall::sfInit(parallel = TRUE, cpus = ncores, type = "SOCK")

		n.percpu <- ceiling( NROW(snp.mat) / ncores );
		snowfall::sfExport("n.percpu", "phe.gls.mat", "snp.mat", "covar.names", "reg.str0", "reg.str1","snp.idx" );

		gls.cluster <- snowfall::sfClusterApplyLB( 1:ncores, cpu.fun);
		snowfall::sfStop();

		if( .RR("debug") )
			cat("Stopping parallel computing......\n");

		r.gls <- do.call( rbind, gls.cluster );
	}
	else
	{
		cat("Starting the fGWAS estimate for each SNP......\n");
		n.percpu <- NROW(snp.mat);
		r.gls <- cpu.fun(1);
	}

	colnames(r.gls) <- c("INDEX", "NAME", "CHR", "POS", "NMISS", "MAF", "LR", "pv");
	rownames(r.gls) <- r.gls$name;

	return(list(error=F, r=r.gls));
}

bls.fgwas <- function( sect.idx, snp.idx, snp.mat, pheY, pheX = NULL, intercept=FALSE, ncores = 1)
{
	sample.ids <- intersect( rownames(pheY), colnames(snp.mat)[-c(1:3)] );
	if(length(sample.ids)==0)
	{
		cat("! No same population data in the phenotypic data and genotypic data.\n")
		return(list(error=T, err.info="No same population data in the phenotypic data and genotypic data."))
	}

	snp.mat <- cbind( snp.mat[, c(1:3)], snp.mat[, sample.ids, drop=F] );
	if(is.null(pheX))
	{
		phe.mat <- data.frame(ID = c(1:NROW(pheY)), rowMeans(pheY, na.rm=T) );
		colnames(phe.mat) <- c("ID",  "Y");
	}
	else
	{
		phe.mat <- data.frame(ID = c(1:NROW(pheY)), pheX, rowMeans(pheY, na.rm=T) );
		colnames(phe.mat) <- c("ID", colnames(pheX), "Y");
	}

	phe.missing <- which( is.na(phe.mat[, NCOL(phe.mat)] ));
	if(length(phe.missing)>0)
		phe.mat <- phe.mat[-phe.missing, , drop=F ];

	reg.str0 <- ifelse( is.null(colnames(pheX)), "Y ~ 1", paste( "Y ~", ifelse(intercept, "1+", ""), paste(colnames(pheX),collapse= "+") )  );
	reg.str1 <- ifelse( is.null(colnames(pheX)), "Y ~ 1 + as.factor(SNP)", paste( "Y ~", ifelse(intercept, "1+", ""), paste(colnames(pheX),collapse= "+"), "+ as.factor(SNP)") );

	if( sect.idx ==1)
	{
		cat("* H0 =", as.character(reg.str0), "\n" );
		cat("* H1 =", as.character(reg.str1), "\n" );
	}

	r0 <- try( do.call("gls", args = list(as.formula(reg.str0), phe.mat, method="ML" ) ) );
	if(class(r0)=="try-error")
	{
		cat("! Failed to call gls() method.\n")
		return(list( error=T, err.info="Failed to call gls() method.") )
	}

	cpu.fun<-function( sect )
	{
		if( (sect-1)*n.percpu+1 > NROW(snp.mat) )
			return(NULL);

		range.fr <- (sect-1)*n.percpu+1;
		range.to <- sect*n.percpu;
		range.to <- ifelse( range.to > NROW(snp.mat),  NROW(snp.mat), range.to );

		r.bls <- data.frame( INDEX=numeric(), NAME=character(),CHR=numeric(), POS=numeric(), NMISS=numeric(), MAF=numeric(), LR=numeric(), pv=numeric(),stringsAsFactors=F);

		reg.f0 <- as.formula(reg.str0);
		reg.f1 <- as.formula(reg.str1);

		for(i in c(range.fr:range.to) )
		{
			snp <- unlist( snp.mat[i, c(4:NCOL(snp.mat)) ] );
			snp.bls <- snp[ phe.mat$ID ];
			bls.mat <- cbind( phe.mat, SNP=snp.bls );

			r.gls[i - range.fr + 1, "NAME"] <- as.character(snp.mat[i,1]);

			maf <- mean(snp.bls, na.rm=T)/2;
			if(maf>0.5) maf <- 1-maf;

			# if maf closes to 0, it will lead to an error: contrasts not defined for 0 degrees of freedom
			pv.max <- NA;
			if( maf <= 10^(-4)) pv.max <- 1.0;

			# remove NA row
			na.row <- unique( which(is.na(bls.mat))%% NROW(bls.mat) );
			if(length(na.row)>0) bls.mat <- bls.mat[-na.row,,drop=F];


			r0 <- try( do.call("gls", args = list( reg.f0, bls.mat, method="ML" ) ) );
			r1 <- try( do.call("gls", args = list( reg.f1, bls.mat, method="ML" ) ) );

			if(any(class(r0)=="try-error") || any(class(r1)=="try-error"))
			{
				r.bls[i-range.fr+1,-c(2)] <- c( snp.idx[i], unlist(snp.mat[i,2:3]), length(na.row), maf, NA, pv.max );
				next;
			}

			r <- try( anova( r0,r1 ) );
			if(any(class(r)=="try-error"))
			{
				r.bls[i-range.fr+1,-c(2)] <- c( snp.idx[i], unlist(snp.mat[i,2:3]), length(na.row), maf, NA, pv.max );
				next;
			}

			r.bls[(i-range.fr+1),-c(2)] <- c( snp.idx[i], unlist(snp.mat[i,2:3]), length(na.row), maf, r[2,8], r[2,9]);
		}

		return(r.bls);
	}

	r.bls <- c();
	if( ncores>1 && require("snowfall") )
	{
		if( .RR("debug") )
			cat("\n  Starting parallel computing, snowfall/snow......\n");

		snowfall::sfInit(parallel = TRUE, cpus = ncores, type = "SOCK")

		n.percpu <- ceiling( NROW(snp.mat) / ncores );
		snowfall::sfExport("n.percpu", "phe.mat", "snp.mat", "reg.str0", "reg.str1", "snp.idx" );

		bls.cluster <- snowfall::sfClusterApplyLB( 1:ncores, cpu.fun);
		snowfall::sfStop();

		if( .RR("debug") )
			cat("  Stopping parallel computing......\n\n");

		r.bls <- do.call( rbind, bls.cluster );

	}
	else
	{
		cat("  Starting the fGWAS estimate for each SNP......\n");
		n.percpu <- NROW(snp.mat);
		r.bls <- cpu.fun(1);
	}

	colnames(r.bls) <- c("INDEX", "NAME", "CHR", "POS", "NMISS", "MAF", "LR", "pv");
	rownames(r.bls) <- r.bls$name;

	return(list(error=F, r=r.bls));

}

