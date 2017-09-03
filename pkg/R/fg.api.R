cur<-function()
{
	return(FG_ENV);
}

fg.simulate<-function( curve.type, covariance.type, n.obs, n.snp, time.points, par0=NULL, par1=NULL, par2=NULL, par.covar=NULL, par.X=NULL,
		sig.pos=NULL, phe.missing=0.03, snp.missing=0.03, plink.format=FALSE, file.prefix=NULL )
{
	if (is.null( fg.getCurve( curve.type ) ))
		stop("No curve object specified by the parameter 'curve.type'. ")

	if (is.null( fg.getCovariance( covariance.type ) ))
		stop("No covariance object specified by the parameter 'covariance.type'. ")

	obj <- fg_simulate( curve.type, covariance.type, n.obs, n.snp, time.points, par0, par1, par2, par.covar, par.X,
		phe.missing, snp.missing, sig.pos, plink.format, file.prefix )
	return( obj );
}

#optional items: verbose=F
#                forece.split=T

fg.load.plink<-function( file.plink.bed, file.plink.bim, file.plink.fam, plink.command=NULL, chr=NULL, options=list(verbose=F))
{
	## dont load plink data, just plink data object
	obj <- fg_load_plink( file.plink.bed, file.plink.bim, file.plink.fam, plink.command, chr, options );

	return( obj );
}

#optional items: verbose=F
fg.load.simple <- function( file.simple.snp, options=list(verbose=F))
{
	obj <- fg_load_simple( file.simple.snp, options )
	return( obj );
}

#optional items: verbose=F
#                max.optim.failure= 100
#                min.optim.success= 20,
#                R2.loop = 5,
fg.load.phenotype<-function( file.phe.long, file.phe.cov=NULL, file.phe.time=NULL, curve.type=NULL, covariance.type=NULL, file.plot.pdf=NULL, intercept=TRUE, options=list(verbose=F))
{
	if (!is.null(curve.type) && tolower(curve.type)!="auto" && is.null( fg.getCurve( curve.type ) ))
		stop("No curve object specified by the parameter 'curve.type'. ")

	if (!is.null(covariance.type) && tolower(covariance.type)!="auto" && is.null( fg.getCovariance( covariance.type ) ))
		stop("No covariance object specified by the parameter 'covariance.type'. ")

	obj <- fg_load_phenotype( file.phe.long, file.phe.cov, file.phe.time, curve.type, covariance.type, file.plot.pdf, intercept, options);
	return( obj );
}

#optional items: verbose=F
#                min.optim.failure= 100
#                min.optim.success= 20,
#                R2.loop = 5,
fg.data.estimate<-function( obj.phe, curve.type="auto", covariance.type="auto", file.plot.pdf=NULL, options=list(verbose=F) )
{
	if (!is.null(curve.type) && tolower(curve.type)!="auto" && is.null( fg.getCurve( curve.type ) ))
		stop("No curve object specified by the parameter 'curve.type'. ")

	if (!is.null(covariance.type) && tolower(covariance.type)!="auto" && is.null( fg.getCovariance( covariance.type ) ))
		stop("No covariance object specified by the parameter 'covariance.type'. ")

	obj <- fg_dat_est( obj.phe, curve.type, covariance.type, file.plot.pdf, options );
	return( obj );
}

#optional items: verbose=FALSE
#                ncores=1,
#                use.snowfall=TRUE
#                max.optim.failure=20 for fgwas and optim-fgwas
#                min.optim.success=2 for fgwas and optim-fgwas
#                use.gradient=T for fgwas and optim-fgwas
#                degree=3 for mixed
fg.snpscan <-function( fgwas.gen.obj, fgwas.phe.obj, method="optim-fgwas", curve.type=NULL, covariance.type=NULL, snp.sub=NULL, options=list(verbose=F))
{
	if( ! toupper(method) %in% toupper(c("gls", "mixed", "fast", "fgwas", "optim-fgwas")))
		stop("the parameter 'method' has 5 optional values: 'gls', 'mixed', 'fast', 'fgwas', 'optim-fgwas'. ")

	if (!is.null(curve.type) && tolower(curve.type)!="auto" && is.null( fg.getCurve( curve.type ) ))
		stop("No curve object specified by the parameter 'curve.type'. ")

	if (!is.null(covariance.type) && tolower(covariance.type)!="auto" && is.null( fg.getCovariance( covariance.type ) ))
		stop("No covariance object specified by the parameter 'covariance.type'. ")

	ret <- fg_snpscan ( fgwas.gen.obj, fgwas.phe.obj, method=method, curve.type=curve.type, covariance.type=covariance.type,permutation=NULL, snp.sub=snp.sub, options=options );
	return( ret );
}

fg.select.sigsnp <- function( fgwas.scan.obj, sig.level=0.05, pv.adjust="bonferroni",  options=list() )
{
	ret <- fg_select_sigsnp( fgwas.scan.obj, sig.level, pv.adjust,  options=list() );
	return(ret);
}

summary.fgwas.gen.obj <- function(object, ...)
{
	stopifnot(class(object)=="fgwas.gen.obj");
	summary_fgwas_gen_obj(object)
}


summary.fgwas.phe.obj <- function(object, ...)
{
	stopifnot(class(object)=="fgwas.phe.obj");
	summary_fgwas_phe_obj(object)
}

summary.fgwas.scan.obj <- function(object, ...)
{
	stopifnot(class(object)=="fgwas.scan.obj");
	summary_fgwas_scan_obj(object)
}

print.fgwas.gen.obj <- function(x, ..., useS4 = FALSE)
{
	stopifnot(class(x)=="fgwas.gen.obj");
	print_fgwas_gen_obj(x  )
}

print.fgwas.phe.obj <- function(x,..., useS4 = FALSE)
{
	stopifnot(class(x)=="fgwas.phe.obj");
	print_fgwas_phe_obj(x )
}

print.fgwas.scan.obj <- function(x, ..., useS4 = FALSE )
{
	stopifnot(class(x)=="fgwas.scan.obj");
	print_fgwas_scan_obj(x )
}

plot.fgwas.scan.obj <- function(x, y, ..., file.pdf=NULL)
{
	stopifnot(class(x)=="fgwas.scan.obj");
	plot_fgwas_scan_obj(x, y, file.pdf, ...)
}

plot.fgwas.phe.obj <- function(x, y, ..., curve.fitting=T, file.pdf=NULL)
{
	stopifnot(class(x)=="fgwas.phe.obj");
	plot_fgwas_phe_obj(x, file.pdf, curve.fitting, ...)
}

plot.fgwas.curve <- function( object, snp.sub, file.pdf=NULL )
{
	stopifnot(class(object)=="fgwas.scan.obj");

	if(is.null(object$ret.fast) && is.null(object$ret.fgwas))
		stop("No curve information in the result, only for FAST and fGWAS model");

	plot_fgwas_curve( object, snp.sub, file.pdf );
}

#Inner function, not public
profile.fgwas.curve <- function( object, snp.sub )
{
	stopifnot(class(object)=="fgwas.scan.obj");

	if(is.null(object$ret.fast) && is.null(object$ret.fgwas))
		stop("No curve information in the result, only for FAST and fGWAS model");

	r <- profile_fgwas_curve( object, snp.sub );
	return(r);
}
