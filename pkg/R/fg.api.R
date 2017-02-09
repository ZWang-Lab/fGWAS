cur<-function()
{
	return(FG_ENV);
}

fg.simulate<-function( curve.type, covariance.type, n.obs, n.snp, time.points, par0=NULL, par1=NULL, par2=NULL, par.covar=NULL, par.X=NULL,
		sig.pos=NULL, phe.missing=0.03, snp.missing=0.03, plink.format=FALSE, file.prefix=NULL )
{
	obj <- fg_simulate( curve.type, covariance.type, n.obs, n.snp, time.points, par0, par1, par2, par.covar, par.X,
		phe.missing, snp.missing, sig.pos, plink.format, file.prefix )

	return( obj );
}

fg.load.cvf<-function( file.plink.bed, file.plink.bim, file.plink.fam, plink.path, options=list())
{
	obj <- fg_load_cvf( file.plink.bed, file.plink.bim, file.plink.fam, plink.path, options );

	return( obj );
}


fg.load.plink<-function( file.plink.bed, file.plink.bim, file.plink.fam, plink.command=NULL, chr=NULL, options=list())
{
	obj <- fg_load_plink( file.plink.bed, file.plink.bim, file.plink.fam, plink.command, chr, options );

	## dont load plink data, just plink data object
	return( obj );
}

fg.load.simple <- function( file.simple.snp, options=list())
{
	obj <- fg_load_simple( file.simple.snp, options )

	return( obj );
}

fg.load.phenotype<-function( file.phe.long, file.phe.cov, file.phe.time, curve.type=NULL, covariance.type=NULL, file.plot.pdf=NULL, intercept=T, options=list())
{
	obj <- fg_load_phenotype( file.phe.long, file.phe.cov, file.phe.time, curve.type, covariance.type, file.plot.pdf, intercept, options);

	return( obj );
}

fg.snpscan <-function( fgwas.gen.obj, fgwas.phe.obj, method="optim-fgwas", curve.type=NULL, covariance.type=NULL, permutation=NULL, snp.sub=NULL, options=list())
{
	if( ! toupper(method) %in% toupper(c("gls", "mixed", "fast", "fast-norm", "fgwas", "optim-fgwas")))
		stop("the parameter 'method' has 5 optional values: 'gls', 'mixed', 'fast', 'fgwas', 'optim-fgwas'. ")

	ret <- fg_snpscan ( fgwas.gen.obj, fgwas.phe.obj, method=method, curve.type=curve.type, covariance.type=covariance.type,permutation=permutation, snp.sub=snp.sub, options=options );
	return( ret );
}

fg.permutation <-function( fgwas.scan.obj, permutation=NULL, options=list() )
{
	ret <- fg_permutation( fgwas.scan.obj, permutation, options );
	return( ret );
}

fg.report <-function( fgwas.scan.obj, file.pdf, options=list() )
{
	ret <- fg_report( fgwas.scan.obj, file.pdf, options );
	return(ret);
}

fg.select.sigsnp <- function( fgwas.scan.obj, sig.level=0.05, pv.adjust="bonferroni",  options=list() )
{
	ret <- fg_select_sigsnp( fgwas.scan.obj, sig.level, pv.adjust,  options=list() );
	return(ret);
}

summary.fgwas.gen.obj <- function(object)
{
	stopifnot(class(object)=="fgwas.gen.obj");
	summary_fgwas_gen_obj(object)
}


summary.fgwas.phe.obj <- function(object)
{
	stopifnot(class(object)=="fgwas.phe.obj");
	summary_fgwas_phe_obj(obejct)
}

summary.fgwas.scan.obj <- function(object)
{
	stopifnot(class(object)=="fgwas.scan.obj");
	summary_fgwas_scan_obj(obejct)
}

print.fgwas.gen.obj <- function(object, useS4 = FALSE)
{
	stopifnot(class(object)=="fgwas.gen.obj");
	print_fgwas_gen_obj(object);
}

print.fgwas.phe.obj <- function(object,useS4 = FALSE)
{
	stopifnot(class(object)=="fgwas.phe.obj");
	print_fgwas_phe_obj(object);
}

print.fgwas.scan.obj <- function(object,useS4 = FALSE)
{
	stopifnot(class(object)=="fgwas.scan.obj");
	print_fgwas_scan_obj(object);
}

plot.fgwas.scan.obj <- function(x, y, ..., file.pdf=NULL)
{
	stopifnot(class(x)=="fgwas.scan.obj");
	plot_fgwas_scan_obj(x, y, file.pdf, ...)
}

plot.fgwas.phe.obj <- function(x, y, ..., file.pdf=NULL)
{
	stopifnot(class(x)=="fgwas.phe.obj");
	plot_fgwas_phe_obj(x, y, file.pdf, ...)
}


plot.fgwas.curve <- function( object, snp.sub, file.pdf=NULL )
{
	stopifnot(class(object)=="fgwas.scan.obj");

	if(is.null(object$ret.fast) && is.null(object$ret.fgwas))
		stop("No curve information in the result, only for FAST and fGWAS model");

	plot_fgwas_curve( object, snp.sub, file.pdf );
}

profile.fgwas.curve <- function( object, snp.sub )
{
	stopifnot(class(object)=="fgwas.scan.obj");

	if(is.null(object$ret.fast) && is.null(object$ret.fgwas))
		stop("No curve information in the result, only for FAST and fGWAS model");

	r <- profile_fgwas_curve( object, snp.sub );
	return(r);
}
