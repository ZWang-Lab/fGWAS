## PLINK format
#
## BIM file:
# Chromosome
# Marker ID
# Genetic distance
# Physical position
# Allele 1
# Allele 2
#
#
#
# Example of a BIM file of the binary PLINK format:
#
# 21	rs11511647	0	26765	A	T
# X	     rs3883674	0	32380	C	G
# X	    rs12218882	0	48172	T	T
# 9	    rs10904045	0	48426	A	T
# 9	    rs10751931	0	49949	C	T
# 8	    rs11252127	0	52087	A	C
# 10	rs12775203	0	52277	A	A
# 8	    rs12255619	0	52481	G	T
#
##
## FAM files
#
# Family ID
# Sample ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Affection (0=unknown; 1=unaffected; 2=affected)
#
#
# FAM1	NA06985	0	0	1	1
# FAM1	NA06991	0	0	1	1
# 0		NA06993	0	0	1	1
# 0		NA06994	0	0	1	1
# 0		NA07000	0	0	2	1
# 0		NA07019	0	0	1	1
# 0		NA07022	0	0	2	1
# 0		NA07029	0	0	1	1
# FAM2	NA07056	0	0	0	2
# FAM2	NA07345	0	0	1	1


fg_load_plink<-function(file.plink.bed, file.plink.bim, file.plink.fam, plink.command=NULL, chr=NULL, options )
{
	default_options <- list( force.split=TRUE, verbose=FALSE);
	default_options[names(options)] <- options;
	options <- default_options;

	if ( missing(file.plink.bed) || missing(file.plink.bim) || missing(file.plink.fam) )
		stop("! file.plink.bed, file.plink.bim, file.plink.fam must be assigned with the valid values.");

	if ( !(is.logical(options$force.split) && length(options$force.split)==1 ) )
		stop("! The parameter of force.split should be a logical value(TRUE or FALSE).");

	if(options$verbose)
	{
		cat( "[ Loading PLINK ] \n");
		cat( "Checking the parameters ......\n");

		cat("* PLINK BED File = ",  file.plink.bed, "\n");
		cat("* PLINK BIM File = ",  file.plink.bim, "\n");
		cat("* PLINK FAM File = ",  file.plink.fam, "\n");
		cat("* Chromosome     = ",  chr, "\n");
		cat("* PLINK Command = ",   plink.command, "\n");
		cat("* Force Split by PLINK Command = ", options$force.split, "\n")
	}


	chk <- check_plink_file(file.plink.bed, file.plink.bim, file.plink.fam, options$verbose);
	if(chk$error)
		stop("! Failed to check PLINK data files.");

	obj.gen <- list();
	obj.gen$reader <- create.plink.obj ( file.plink.bed, file.plink.bim, file.plink.fam, plink.command, chromosome=chr, verbose=options$verbose);

	obj.gen$n.snp  <- obj.gen$reader$n.snp;
    obj.gen$n.ind.total <- obj.gen$reader$n.ind.total;
    obj.gen$n.ind.used  <- length(obj.gen$reader$ids.used);

	obj.gen$options = options;
	obj.gen$params  = list(file.plink.bed = file.plink.bed,
				file.plink.bim = file.plink.bim,
				file.plink.fam = file.plink.fam,
				plink.command = plink.command );

	class(obj.gen) <- "fgwas.gen.obj";

	return(obj.gen);

}

check_plink_file<-function( file.plink.bed, file.plink.bim, file.plink.fam, verbose=FALSE )
{
	if(verbose)
	{
		cat("Checking PLINK file......\n");
		cat("* BED file =", file.plink.bed, "\n");
		cat("* BIM file =", file.plink.bim, "\n");
		cat("* FAM file =", file.plink.fam, "\n");
	}

	bigdata  <- FALSE;
	error <- FALSE;

	snp.mat <- try( read.plink( file.plink.bed,  file.plink.bim, file.plink.fam) );
	if(class(snp.mat)=="try-error")
	{
		if(all( file.exists( file.plink.bed,  file.plink.bim, file.plink.fam ) ) )
			bigdata <- TRUE
		else
			error <- TRUE;
	}

	if(!error)
	{
		tb.fam <- read.table(file.plink.fam, header=F);
		tb.bim <- read.table(file.plink.bim, header=F);
		n.idv <- NROW(tb.fam);
		n.snp <- NROW(tb.bim);

		if(verbose)
		{
			cat("* Individuals =", n.idv, "SNPs=", n.snp, "\n")
			cat("* PLINK testing successfully.\n")
		}
		return(list(error=F, bigdata=bigdata, family=tb.fam[,2]));
	}
	else
		return(list(error=T));

}


create.plink.obj  <- function( file.plink.bed, file.plink.bim, file.plink.fam, plink.command, chromosome=NULL, verbose=FALSE )
{
	if ( missing( file.plink.bed) || missing(file.plink.bim) || missing(file.plink.fam) )
		stop("! file.plink.bed, file.plink.bim and file.plink.fam must be assigned with the valid values.");

	if(is.null(plink.command)) plink.command <- "plink";

	fg.dm.plink <- getRefClass("fg.dm.plink");
	objref <- fg.dm.plink(
		file.plink.bed    = file.plink.bed,  ##
		file.plink.bim    = file.plink.bim,  ##
		file.plink.fam    = file.plink.fam,  ##
		plink.command     = ifelse(is.null(plink.command),"", plink.command),   ##
		chromosome        = ifelse(is.null(chromosome),-1, chromosome),
		snp.blocksize     = 1000,
		snpdata           = list());

	t <- try(system( paste( plink.command, "--noweb", sep=" "), ignore.stdout=TRUE, ignore.stderr=TRUE ));
	if(class(t)=="try-error")
	{
		cat("! No PLINK command can be found in your environment( plink.command=",plink.command, ")\n")
		return(NULL);
	}

	plink.checkFiles( objref, verbose );

	return(objref);
}

plink.checkFiles <- function( refobj, verbose=FALSE )
{
	if(verbose)
	{
		cat( "Checking PLINK files ......\n");
		cat("* PLINK BED File = ",  refobj$file.plink.bed, "\n");
		cat("* PLINK BIM File = ",  refobj$file.plink.bim, "\n");
		cat("* PLINK FAM File = ",  refobj$file.plink.fam, "\n")
		cat("* PLINK Command = ",   refobj$plink.command, "\n")
	}

	if( all( file.exists( refobj$file.plink.bed,  refobj$file.plink.bim, refobj$file.plink.fam ) ) )
	{
		refobj$snpdata$plink.bim <- read.table(refobj$file.plink.bim, header=F, stringsAsFactors=F)
		refobj$snpdata$plink.fam <- read.table(refobj$file.plink.fam, header=F, stringsAsFactors=F)

		if(refobj$chromosome != -1)
		{
			snp.idx <- which(refobj$snpdata$plink.bim[,1] %in% refobj$chromosome);
			refobj$snpdata$plink.bim <- refobj$snpdata$plink.bim[snp.idx, ]
		}

		refobj$n.snp        <- NROW(refobj$snpdata$plink.bim);
		refobj$n.ind.total  <- NROW(refobj$snpdata$plink.fam);
		refobj$n.ind.used   <- refobj$n.ind.total ;
		refobj$ids.used     <- as.character(refobj$snpdata$plink.fam[,2])
	}
	else
		stop("PLInK files can not be found!");
}

plink.get.snpmat<-function(objref, snp.idx=NULL, impute=F, allel=F )
{
	snp.idx.ord <- order( snp.idx, decreasing=F);
	snp.idx.new <- snp.idx[ snp.idx.ord ];

	check_local_snpmat<-function(snp.k)
	{
		if(is.null(objref$snpdata$local.idx))
			return(FALSE);

		return(snp.k %in% objref$snpdata$local.idx);
	}

	load_local_snpmat<-function(snp.k)
	{
		select.subjects <- objref$snpdata$plink.fam[ match( objref$ids.used, objref$snpdata$plink.fam[,2] ), c(1,2)]
		loading.idx <- seq( snp.k, ifelse(snp.k + objref$snp.blocksize > objref$n.snp, objref$n.snp, snp.k + objref$snp.blocksize), 1);
		plink.obj <- plink.cmd.select( objref$plink.command, objref$file.plink.bed, objref$file.plink.bim, objref$file.plink.fam,
						select.snps = objref$snpdata$plink.bim[loading.idx, 2],
						select.subjects = select.subjects )

		if(class(plink.obj)=="try-error")
			stop("Error in PLINK data extraction.");

		objref$snpdata$local.idx   <- loading.idx;
		objref$snpdata$local.fam   <- plink.obj$fam;
		objref$snpdata$local.map   <- plink.obj$map;
		objref$snpdata$local.snpmat<- plink.obj$genotypes;
	}

	get_local_snpmat<-function(snp.k)
	{
		k.idx <- which(snp.k==objref$snpdata$local.idx);
		return(as.numeric(objref$snpdata$local.snpmat[,k.idx]) -1 );
	}

	snpmat <- matrix(NA, nrow = objref$n.ind.used, ncol=NROW(snp.idx.new));
	for(i in 1:NROW(snp.idx.new))
	{
		if ( !check_local_snpmat( snp.idx.new[i] ) )
			load_local_snpmat( snp.idx.new[i] )

		snpmat[,i] <- get_local_snpmat( snp.idx.new[i] )
		snpmat[ snpmat[,i]==-1 , i ]  <- NA;
		if (mean(snpmat[,i], na.rm=T)/2>0.5 )
			snpmat[,i] <- 2 - snpmat[,i];
	}

	snpmat <- snpmat[, order(snp.idx.ord), drop=F];

	if( impute )
	{
		snpmat.imp <- impute_simple_snp(snpmat);
		return(list(snpmat = snpmat.imp$snpmat, NMISS=snpmat.imp$NMISS, MAF=snpmat.imp$MAF))

	}
	else
	{
		NMISS <- unlist(apply(snpmat, 2, function(snp.k){ length(which(is.na(snp.k)));}));
		return(list(snpmat = snpmat, NMISS=NMISS, MAF=colMeans(snpmat, na.rm=T)/2));
	}
}

plink.get.snpinfo<-function(objref, snp.idx)
{
	# remove genetic distance at 3rd postion
	return( objref$snpdata$plink.bim[snp.idx,c(2,1,4,5,6), drop=F]);
}

plink.get.individuals<-function(objref)
{
	return( objref$snpdata$plink.fam[,2]);
}

plink.select.individuals<-function(objref, ids.used)
{
	objref$ids.used <- as.character(ids.used);
}

plink.get.used.individuals<-function(objref)
{
	return(objref$ids.used);
}

plink.get.snpindex <- function( objref, snp.names )
{
	return( match(snp.names, objref$snpdata$plink.bim[,2]) );
}

plink.shrink <-function(objref)
{
	objref$snpdata$plink.snpmat <- NULL;
	objref$snpdata$plink.fam <- NULL;
	objref$snpdata$plink.map <- NULL;
}

test_plink_func<-function(plink)
{
	n.snp <- NCOL( plink$genotypes );

	cat("...Get 100 SNPs\n");
	r.snp <- get_plink_subsnp(plink, 1:100);

	cat("...Get 1000 SNPs\n");
	r.snp <- get_plink_subsnp(plink, 1:1000);

	cat("...Get 10000 SNPs\n");
	r.snp <- get_plink_subsnp(plink, 1:10000);

	cat("...Get 50000 SNPs\n");
	r.snp <- get_plink_subsnp(plink, 1:50000);

	cat("...Get 100000 SNPs\n");
	r.snp <- get_plink_subsnp( plink, sample(n.snp)[1:100000] );
}

get_plink_subsnp<-function(snp.mat, snp.set.idx)
{
	s.mat <- as( snp.mat$genotypes[, snp.set.idx, drop=F ], "numeric");
	snp.imp <-c();
	snp.maf <- c();
	snp.names <- c();

	f.impute<-function(s.mat.i )
	{
		s.miss <- which( is.na(s.mat.i) );

		if (length(s.miss)>0)
		{
			n.AA <- length( which( s.mat.i == 0 ) );
			n.Aa <- length( which( s.mat.i == 1 ) );
			n.aa <- length( which( s.mat.i == 2 ) );
			n.s  <- n.AA + n.Aa + n.aa;

			r.miss <- runif( length(s.miss) );
			r.snp  <- rep(2, length(s.miss));
			r.snp [ r.miss <= n.AA/n.s ]<-0;
			r.snp [ r.miss <= (n.AA + n.Aa)/n.s ]<-1;
			s.mat.i[s.miss] <- r.snp;
		}

		if (mean(s.mat.i)/2>0.5) s.mat.i <- 2 - s.mat.i;

		return( s.mat.i );
	}


	snp.imp <- t( apply(s.mat, 2, f.impute) );
	snp.maf <- rowMeans(snp.imp) /2;

	map <- snp.mat$map[snp.set.idx, ,drop=F];
	rownames(snp.imp) <-rownames( map );

	return(list(maf=snp.maf, snp=snp.imp, info=map[,c(2,1,4)]) );
}

impute_simple_snp <- function( snpmat )
{
	f_imputed <- function( snp )
	{
		s.miss <- which( is.na( snp ) );
		if ( length(s.miss)>0 )
		{
			n.AA <- length(which( snp == 0 ) );
			n.Aa <- length(which( snp == 1 ) );
			n.aa <- length(which( snp == 2 ) );
			n.s  <- n.AA + n.Aa + n.aa;

			r.miss <- runif( length(s.miss) );
			r.snp  <- rep(2, length( s.miss ));
			r.snp [ r.miss <= n.AA/n.s ] <- 0;
			r.snp [ r.miss <= (n.AA + n.Aa)/n.s ] <-1;
			snp [ s.miss ] <- r.snp;
		}

		if (mean(  snp )/2>0.5)
			snp <- 2 - snp;

		return(list(snp=snp, NMISS=length(s.miss), MAF=mean(  snp )/2 ));
	}

	total_miss <- length(which(is.na(snpmat)));

	r.imp  <- apply(snpmat, 2, f_imputed);
	snpmat <- do.call("cbind", lapply(1:NCOL(snpmat), function(i){return(r.imp[[i]]$snp);}));
	NMISS  <- do.call("unlist", lapply(1:NCOL(snpmat), function(i){return(r.imp[[i]]$NMISS);}));
	MAF    <- do.call("unlist", lapply(1:NCOL(snpmat), function(i){return(r.imp[[i]]$MAF);}));


	cat("* Missing SNPs are imputed(", total_miss, "SNPs).\n");

	return(list(snpmat=snpmat, NMISS=NMISS, MAF=MAF));
}

plink.cmd.load.chromosome<-function(plink.command,file.plink.bed, file.plink.bim, file.plink.fam, chr)
{
	tmp <- tempfile();
	str.cmd <- paste( plink.command, "--noweb",
					"--bed", file.plink.bed,
					"--bim", file.plink.bim,
					"--fam", file.plink.fam,
					"--chr", chr,
					"--make-bed",
					"--out", tmp,
					sep=" ");

	t <- try(system( str.cmd, ignore.stdout=TRUE, ignore.stderr=TRUE) );
	if(class(t)=="try-error")
	{
		cat("! Error in PLINK command.\n")
		return(list(error=T, err.info="Error in PLINK command."));
	}

	plink <- try( read.plink( paste(tmp, ".bed", sep=""),  paste(tmp, ".bim", sep=""), paste(tmp, ".fam", sep="")) );
	if(class(plink) == "try-error")
	{
		cat("! Package snpStats can not open PLINK data.\n")
		return(list(error=T, err.info="Package snpStats can not open PLINK data."));
	}

	return(plink);
}


plink.cmd.select<-function(plink.command, file.plink.bed, file.plink.bim, file.plink.fam, select.snps, select.subjects)
{
	file.keep.snp <- tempfile()
	file.keep.subject <- tempfile()
	write.table(select.snps, file=file.keep.snp, quote=F, row.names=F, col.names=F)
	write.table(select.subjects, file=file.keep.subject, quote=F, row.names=F, col.names=F)

	tmp <- tempfile();
	str.cmd <- paste( plink.command, "--noweb",
					"--bed", file.plink.bed,
					"--bim", file.plink.bim,
					"--fam", file.plink.fam,
					"--keep", file.keep.subject,
					"--extract", file.keep.snp,
					"--make-bed",
					"--out", tmp,
					sep=" ");

	t <- try(system( str.cmd, ignore.stdout=TRUE, ignore.stderr=FALSE) );
	if(class(t)=="try-error")
	{
		cat("! Error in PLINK command.\n")
		return(list(error=T, err.info="Error in PLINK command."));
	}

	plink <- try( read.plink( paste(tmp, ".bed", sep=""),  paste(tmp, ".bim", sep=""), paste(tmp, ".fam", sep="")) );
	if(class(plink) == "try-error")
	{
		cat("! Package snpStats can not open PLINK data.\n")
		return(list(error=T, err.info="Package snpStats can not open PLINK data."));
	}

	unlink(c(file.keep.snp, file.keep.subject, paste(tmp, ".bed", sep=""),  paste(tmp, ".bim", sep=""), paste(tmp, ".fam", sep="") ));

	return(plink);
}


convert_simpe_to_plink <- function( snp.info, snp.mat, snp.file.base )
{
	# snp,mat : 0/1/2/NA
	# PLINK raw data: 1/2/3==> AA,AB,BB, 0==>NA
	snp.mat <- t(snp.mat + 1);
	snp.mat[is.na(snp.mat)] <- 0;

	sub.name <- colnames(snp.mat);
	snp.name <- rownames(snp.mat);

	###snps
	dim.snps <- dim(snp.mat);

	snps <- as.raw( as.matrix(snp.mat ) );
	snps <- array(snps, dim=dim.snps);
	colnames(snps) <- sub.name;
	rownames(snps) <- snp.name;
	class(snps) <- "SnpMatrix";

	r <- write.plink( file.base=snp.file.base, snp.major = F, snps=t(snps),
	    	id=sub.name,
	    	father=rep(0,dim.snps[2]),
	    	mother=rep(0,dim.snps[2]),
	    	sex=rep(0,dim.snps[2]),
	    	phenotype = rep(-9,dim.snps[2]),
			chromosome = as.character(snp.info[,1]),
			genetic.distance = as.numeric(snp.info[,3]),
			position= as.numeric(snp.info[,4]),
			allele.1 = as.character(snp.info[,5]),,
			allele.2 = as.character(snp.info[,6]),
			na.code=0);

	cat(" Genotype files have been converted into PLINK binary format(bed/bim/fam)\n");

	return(list(file.plink.bed = paste(snp.file.base, ".bed", sep=""),
   	    	file.plink.bim = paste(snp.file.base, ".bim", sep=""),
   	    	file.plink.fam = paste(snp.file.base, ".fam", sep="")));
}