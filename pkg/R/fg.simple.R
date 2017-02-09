create.simple.obj  <- function( file.simple.snp, options )
{
	if ( missing( file.simple.snp ) )
		stop("! file.simple.snp  must be assigned with the valid values.");

	fg.dm.simple <- getRefClass("fg.dm.simple");
	objref <- fg.dm.simple(
		file.simple.snp   = file.simple.snp,
		rawdata           = data.frame(Characters=character()),
		snpdata           = list());

	tb.gen <- read.table( objref$file.simple.snp, header=T );
	objref$n.snp <- NROW(tb.gen)
	objref$n.ind.count <- NCOL(tb.gen)-5
	objref$idsused <- colnames(tb.gen)[-c(1:5)];
	objref$n.ind.used <- length(objref$idsused);

	load_simple_snp(objref);

	return(objref);
}

simple.shrink<-function(objref)
{
	return(objref);
}

simple.get.snpmat<-function(objref, snp.idx=NULL, impute=F, allel=F  )
{
	allel.mat  <- NULL
	allel.info <- NULL
	if (is.null(snp.idx))
	{
		allel.mat <- objref$snpdata$snp.mat[, objref$ids.used,drop=F];
		allel.info <- objref$snpdata$snp.info;
	}
	else
	{
		allel.mat  <- objref$snpdata$snp.mat[ snp.idx, objref$ids.used, drop=F];
		allel.info <- objref$snpdata$snp.info[ snp.idx, , drop=F];
	}

	covert_geno <- function( allel.info, allel.geno)
	{
		snpB <- as.character(unlist(allel.info[5]));
		snpA <- as.character(unlist(allel.info[4]));
		QQ2<- paste(snpB, snpB, sep="");
		qq0<- paste(snpA, snpA, sep="");
		Qq1<- c(paste(snpA, snpB, sep=""), paste( snpB, snpA, sep="") ) ;

		d.g <- rep( NA, length(allel.geno) );
		d.g[which(allel.geno == QQ2)]<-2;
		d.g[which(allel.geno == qq0)]<-0;
		d.g[which(allel.geno == Qq1[1])]<-1;
		d.g[which(allel.geno == Qq1[2])]<-1;

		if( mean(d.g, na.rm=T)/2 > 0.5)
			d.g <- 2 - d.g;

		return(d.g);
	}

	snpmat <- do.call( "cbind", lapply(1:NROW(allel.mat), function(i){
		return( covert_geno ( allel.info[i,], as.character(unlist(allel.mat[i,])) ) );
	}));

	rownames(snpmat) <- colnames(allel.mat);
	colnames(snpmat) <- rownames(allel.mat);

	if(!allel)
	{
		if(impute)
		{
			snpmat.imp <- impute_simple_snp( snpmat );
			return(list(snpmat = snpmat.imp$snpmat, NMISS=snpmat.imp$NMISS, MAF=snpmat.imp$MAF))
		}
		else
		{
			NMISS <- unlist(apply(snpmat, 2, function(snp.i){length(which(is.na(snp.i)))}));
			return(list(snpmat = snpmat, NMISS=NMISS, MAF=colMeans(snpmat, na.rm=T)/2));
		}
	}
	else
		return(list(snpmat = allel.mat, NMISS=NA, MAF=NA))
}

simple.get.snpinfo<-function(objref, snp.idx=NULL )
{
	if (is.null(snp.idx))
		return(objref$snpdata$snp.info)
	else
		return(objref$snpdata$snp.info[snp.idx,,drop=F]);
}

simple.get.individuals<-function(objref)
{
	return( colnames(objref$snpdata$snp.mat) );
}

simple.select.individuals<-function(objref, ids.used)
{
	objref$ids.used <- ids.used;

	load_simple_snp(objref);

	ids.all <- colnames(objref$snpdata$snp.mat)
	if(any(is.na(match(ids.used, ids.all))))
		stop("Some IDs are not matached in the SNP data file");

	objref$n.ind.used <- length(ids.used);
	return(objref);
}

simple.get.used.individuals<-function(objref)
{
	return( objref$ids.used );
}

simple.get.snpindex <- function( objref, snp.names )
{
	return( match(snp.names, rownames(objref$snpdata$snp_info) ) );
}


#scaffoldId, Loci, RefBase, AltBase, P1, P2, P3,....
load_simple_snp<-function( objref )
{
	if(objref$file.simple.snp != "" )
		tb.gen <- read.table( objref$file.simple.snp, header=T)
	else
		tb.gen <- objref$rawdata;
		
	cat("gen.table[", dim(tb.gen), "]\n");

	snp_info <- tb.gen[, c(1:5)]
	snp_mat  <- tb.gen[, -c(1:5)]
	rownames(snp_info) <- snp_info[,1];
	objref$snpdata$snp_info <- snp_info

	ids.idx <- match( objref$ids.used, colnames(snp_mat) );
	objref$snpdata$snp_mat <- snp_mat[, ids.idx];
}

fg_load_simple<-function( file.simple.snp, options )
{
	cat( "[ Loading Simple ] \n");
	cat( "Checking the parameters ......\n");

	if ( missing(file.simple.snp) )
		stop("! file.simple.snp must be assigned with the valid file name..");

	cat("* SIMPLE SNP File = ",  file.simple.snp, "\n");

	if(!file.exists( file.simple.snp))
		stop("Failed to open Simple data files.");

	params <- list( file.simple.snp = file.simple.snp );

	obj.gen <- create.simple.obj ( file.simple.snp, options );

	ret <- list( reader=obj.gen, options=options, params=params );
	ret$n.snp <- obj.gen$n.snp
	ret$n.ind.total <- obj.gen$n.ind.count
	ret$n.ind.used <- obj.gen$n.ind.used
	class(ret) <- "fgwas.gen.obj";

	return(ret);

}


