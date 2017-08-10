fg_mixed_scan<- function( obj.gen, obj.phe, order=3, snp.idx = NULL, ncores=1 )
{
	if( !requireNamespace("lme4", quietly=T) )
		stop("Package lme4 is not installed, please use 'install.packages' command to install it.");

	F.fisher<-function(p)
	{
		return( pchisq( -2 * sum(log(p)), df = 2 * length(p), lower.tail=F) );
	}

	Fast1 <- function( b.mat, snp.i )
	{
		p<-unlist(lapply(0:order, function(i){
			if(order>=i)
			{
				data = data.frame(Y = b.mat[,i+1], SNP = factor(snp.i));
				colnames(data) <- c("Y", "SNP");
				r <- lm(Y~factor(SNP), data = data);
				f <- summary(r)$fstatistic;
				p1 <- pf(f[1],f[2],f[3],lower.tail=F);
				attributes(p1) <- NULL;
				return(p1);
			}
			else
				return(1);
		} ) );

		return( c( F.fisher(p), min(p), exp(sum(log(p))), p ));
	}

	if(missing(snp.idx) || is.null(snp.idx))
		snp.idx <- c(1:obj.gen$n.snp);

	pheY = obj.phe$pheY;
	pheX = obj.phe$pheX;
	pheT = obj.phe$pheT;

	intercept <- ifelse ( !is.null(obj.phe$intercept), obj.phe$intercept, FALSE);

	if(!is.null(pheX))
	{
		colnames(pheX) <- paste("X", 1:NCOL(pheX), sep="");
		df <- cbind( UID = c(1:NROW(pheX)), pheX );
	}
	else
		df <- data.frame( UID = c(1:NROW(pheY)) );

	phe.one<- do.call("rbind", lapply(1:NROW(pheY), function(i){
		phe.indv <- c();
		for(j in 1:NCOL(pheY)){
			if(!is.na(pheY[i,j]))
				phe.indv <-rbind( phe.indv, data.frame(df[i,,drop=F], pheT[i,j], pheY[i,j]));
		}
		return(phe.indv);
	}));

	if(!is.null(pheX))
		colnames(phe.one)<-c("UID", colnames(pheX), "T", "Y")
	else
		colnames(phe.one)<-c("UID", "T", "Y");

	## formulae = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + 1 + I(T^1)+I(T^2)+I(T^3)+I(T^4)+I(T^5) +
	##                         ( 1 + I(T^1)+I(T^2)+I(T^3)+I(T^4)+I(T^5) | UID )
	##
	str.form = "Y ~ "
	if(!is.null(pheX))
		str.form = paste( str.form,  paste("X", 1:NCOL(pheX), collapse="+", sep=""));

	str.form.T = "1";
	for(i in 1:order)
		str.form.T = paste(str.form.T, " + I(T^", i, ") ", sep="");

	reg.form <- as.formula( paste(str.form, "+", str.form.T, "+(", str.form.T, "|UID)", sep=""));
	est1 <- try( do.call("lme4::lmer", args = list( reg.form, data=phe.one) ) );
	if(any(class(est1)=="try-error") )
		list(error=TRUE, err.info="Failed to call lmer.")	;

	est1.beta <- coef(est1)[[1]]
	b0.mat <- eval(scale(est1.beta$"(Intercept)"));
	for(i in 1:order)
		b0.mat <- cbind(b0.mat, eval(parse(text=paste("scale(est1.beta$\"I(T^", i, ")\")", sep="")) ) );

	time <- seq(min(obj.phe$pheT, na.rm=T), max(obj.phe$pheT, na.rm=T), 0.1);
	curve <- b0.mat;

	ret <- mclapply(snp.idx, function(k){
		snp.mat <- obj.gen$reader$get_snpmat( k, impute = T );
		snp.info <- obj.gen$reader$get_snpinfo( k );
		ret <- Fast1( b0.mat, snp.mat$snpmat);

		y.pv <- cbind( Index=k, snp.info[1:5], MAF=snp.mat$MAF[1], NMISS=snp.mat$NMISS[1], p.fisher=ret[1], p.min=ret[2], p.join=ret[3]);
		for(i in 0:order)
			y.pv <- cbind(y.pv, ret[3+i+1]);
		return(y.pv);

	}, mc.cores = ncores );

	ret <- as.data.frame(do.call("rbind", ret ) );

	colnames(ret) <- c("INDEX", "NAME","CHR", "POS", "Allel1", "Allel2", "MAF", "NMISS", "pv","p.min","p.join",paste("p", 0:order, sep="_"));
	rownames(ret) <- ret$NAME;

	return( list(lmer=est1, beta=b0.mat, result = ret, error=FALSE) );
}