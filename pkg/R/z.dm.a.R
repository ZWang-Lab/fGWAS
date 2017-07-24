#' S4 class for storing TFBS.
#' sequence preferences of a set of TFs.


setRefClass("fg.dm.genodata",

	fields = list(
		type           = "character",    ## plink, cvf, simple.
		description    = "character",    ## plink, cvf, simple.
		n.snp          = "numeric",      ## Number of SNP in the data file.
		n.ind.total    = "numeric",      ## Number of individual in the data file.
		n.ind.used     = "numeric",      ## Number of individual used.
		ids.used       = "character"),   ## individual used.

    methods = list(
		show = function()
		{
       		cat("  Reference matrix object of class",	classLabel(class(.self)), "\n")
			cat("  Data type: ",    type, "\n");
			cat("  Description: ",  description, "\n");
			cat("  SNP Count: ",    n.snp, "\n");
			cat("  Individual Count: ", n.ind.total, "\n");
			cat("  Individual Used: ",  n.ind.used, "\n");
	    },

		shrink = function()
		{
	        NA;
		},

		# output data.frame
		#
		# colnames:
		# SNPName(#rs or ##)
		# chromosome
		# position
		# Allel 1
		# Allel 2
		get_snpinfo = function(snp.idx )
		{
	        NA;
		},

		# output list with snpmat(Matrix[n_sub, n_snp]), NMISS(vector), MAF(vector)
		get_snpmat = function( snp.idx, impute=F, allel=F)
		{
	        NA;
		},

		#return object
		select_individuals = function( ids.used)
		{
	        NA;
		},

		#return all individuals
		get_individuals = function()
		{
	        NA;
		},

		#return used individuals
		get_used_individuals = function()
		{
	        NA;
		},

		# output index
		#
		get_snpindex = function( snp.names  )
		{
	        NA;
		})
)

setRefClass("fg.dm.plink",
	fields = list(
		file.plink.bed    = "character",  ##
		file.plink.bim    = "character",  ##
		file.plink.fam    = "character",  ##
		plink.command  	  = "character",  ##
		chromosome        = "numeric",
		snp.blocksize     = "numeric",
		snpdata           = "list"
	),
	contains = "fg.dm.genodata",
    methods = list(
		show = function()
		{
			if (chromosome!=-1) 
				cat("  Chromosome: ", chromosome, "\n")
			else	
				cat("  Chromosome: ", "all", "\n")

			callSuper();
			cat("  Plink Command: ",    plink.command, "\n");
	    },

		shrink = function()
		{
	        return(plink.shrink(.self));
		},

		get_snpinfo = function(snp.idx )
		{
			return(plink.get.snpinfo (.self, snp.idx));
		},

		get_snpmat = function( snp.idx, impute=F, allel=F)
		{
			return(plink.get.snpmat(.self, snp.idx, impute=impute, allel=allel));
		},

		select_individuals = function( ids.used)
		{
			return(plink.select.individuals(.self, ids.used));
		},

		get_used_individuals = function()
		{
			return(plink.get.used.individuals(.self));
		},

		get_individuals = function()
		{
			return(plink.get.individuals(.self));
		},

		get_snpindex = function( snp.names  )
		{
			return(plink.get.snpindex( .self, snp.names ));
		} )
		

)


setRefClass("fg.dm.simple",

	fields=list(
		file.simple.snp   = "character",  ##
		rawdata           = "data.frame",
		snpdata           = "list"
	),
	contains = "fg.dm.genodata",

    methods = list(
		show = function()
		{
			callSuper();
	    },

		shrink = function()
		{
	        return(simple.shrink(.self));
		},

		get_snpinfo = function(snp.idx )
		{
			return(simple.get.snpinfo (.self, snp.idx));
		},

		get_snpmat = function( snp.idx, impute=F, allel=F)
		{
			return(simple.get.snpmat(.self, snp.idx, impute=impute, allel=allel));
		},

		select_individuals = function( ids.used)
		{
			return(simple.select.individuals(.self, ids.used));
		},

		get_used_individuals = function()
		{
			return(simple.get.used.individuals(.self));
		},

		get_individuals = function()
		{
			return(simple.get.individuals(.self));
		},

		get_snpindex = function( snp.names  )
		{
			return(simple.get.snpindex( .self, snp.names ));
		} )

)


setGeneric("select_individuals",
	def = function( object,  ids.used  ){
		standardGeneric("select_individuals");
	})



fgwas.phe.select.individuals<-function( object, ids.used )
{
	idx.match <- match(as.character(ids.used), object$ids);
	if(any(is.na(idx.match)))
		stop("Some IDs are not matached in the SNP data file");

	object$n.ind <- length(ids.used);
	object$ids   <- as.character(ids.used);
	object$pheY  <- object$pheY[idx.match,,drop=F];
	object$pheX  <- object$pheX[idx.match,,drop=F];
	object$pheZ  <- object$pheZ[idx.match,,drop=F];

	return(object);
}

setMethod("select_individuals",  signature( object = "fgwas.phe.obj" ), fgwas.phe.select.individuals )



