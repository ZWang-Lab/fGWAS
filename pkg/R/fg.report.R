###############################################################
#
# report utility
#    1). proc_report_dat
#    2). proc_report_snpscan
#    3). proc_report_perm
#    4). proc_report_sig
#    5). proc_report_topn
#
# History:
# 12/15/2011 Version 1.1
#
###############################################################

proc_report_dat<-function( dat )
{
	if ( class(dat) != "FG.DAT" )
	{
		sErrMsg<- "Error: Not a data set for fGWAS.";
		return(list(error=T, err.info=sErrMsg));
	}

	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Curve:", "",       dat$mu.f$name) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Covariance:", "",  dat$cov.f$name ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Pheno. file:", "", dat$phe.csv ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Geno. file:",  "", dat$geno.dat ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Sample size:", "", dat$n.obs) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("SNP size:",    "", dat$n.snp) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Sample times:","", length(dat$times) ) );

	Report.AddParaBreak();

	c1 <- call("fpt.plot_tiled_curves", dat );
	Report.AddFigure( c1, " ", c(5, 5)*254, left.margin=0.1*254);

	c2 <- call("fpt.plot_overlapping_curves", dat );
	Report.AddFigure( c2, " ", c(5, 4)*254, left.margin=0.1*254);

	Report.AddParaBreak();

	return(list(error=F));
}

proc_report_snpscan<-function( dat, r.snpscan, r.perm=NULL )
{
	if ( is.null(r.snpscan) || class(r.snpscan)!="FG.SCAN" )
	{
		sErrMsg<- "Error: Not a result for fGWAS.";
		return(list(error=T, err.info=sErrMsg));
	}

	p.05 <-NA;
	p.01 <-NA;
	if(!is.null(r.perm) || class(r.perm)!="FG.PERM")
	{
		p.05 <-r.perm$p.05;
		p.01 <-r.perm$p.01;
	}

	c1 <- call("fpt.plot_manhattan", r.snpscan$full, p.05, p.01, "Manhattan" );
	Report.AddFigure( c1, " ", c(6, 3)*254, left.margin=0.1*254);

	Report.AddParaBreak();

	return(list(error=F));

}

proc_report_perm<-function( dat, r.perm )
{
	if ( is.null(r.perm) || class(r.perm)!="FG.PERM" )
	{
		sErrMsg<- "Error: Not a result for fGWAS.";
		return(list(error=T, err.info=sErrMsg));
	}

	c1 <- call("fpt.plot_perm_curve", r.perm$p.cut );
	Report.AddFigure( c1, " ", c(6, 3)*254, left.margin=0.1*254);

	Report.AddParaBreak();

	return(list(error=F));

}

proc_report_snps<-function( dat, snp.res)
{
	par_cov = c();
	par_QQ = c();
	par_Qq = c();
	par_qq = c();

	if( dim(snp.res)[1] < 1)
		return("");

	for(i in 1:dim(snp.res)[1] )
	{
		n.start <- 5 + dat$cov.f$pnum +  dat$mu.f$pnum+1
		par <- round( snp.res[i, c(n.start:(n.start+dat$cov.f$pnum-1))], 2);
		par.str <- paste(par, collapse=",");
		par_cov <- c(par_cov, par.str);

		n.start <- n.start + dat$cov.f$pnum;
		par <- round( snp.res[i, c(n.start:(n.start+dat$mu.f$pnum-1))], 2);
		par.str <- paste(par, collapse=",");
		par_QQ <- c( par_QQ, par.str);

		n.start <- n.start + dat$mu.f$pnum;
		par <- round( snp.res[i, c(n.start:(n.start+dat$mu.f$pnum-1))], 2);
		par.str <- paste(par, collapse=",");
		par_Qq <- c( par_Qq, par.str);

		n.start <- n.start + dat$mu.f$pnum;
		par <- round( snp.res[i, c(n.start:(n.start+dat$mu.f$pnum-1))], 2);
		par.str <- paste(par, collapse=",");
		par_qq <- c( par_qq, par.str);
	}

	sigs_list <- data.frame(
				Grp = c(snp.res[,1, drop=F]),
				Pos = sprintf("%.2f", snp.res[,2, drop=F] ),
				LR  = sprintf("%.2f", snp.res[,4, drop=F] ),
				cov = par_cov,
				QQ  = par_QQ,
				Qq  = par_Qq,
				qq  = par_qq);

	Report.AddTable( sigs_list,
				 title = "Significant SNPs",
				 frame = "wang",
				 col.names = c("Grp.", "Pos.", "LR2", "Covariance", "QQ",  "Qq",  "qq" ),
				 col.width = c(100,    100,    200,    320,          400,   400,   400  ),
				 col.align = c("L",    "L",    "R",    "L",          "L",   "L",   "L"),
				 offset.x = 50,
				 max.show = 20)

	Report.AddParaBreak(50);


	for(i in 1:dim(snp.res)[1] )
	{
		n.start <- 5 + dat$cov.f$pnum +  dat$mu.f$pnum+1
		par.cov <- snp.res[i, c(n.start:(n.start+dat$cov.f$pnum-1))];
		n.start <- n.start + dat$cov.f$pnum;
		par.QQ  <- snp.res[i, c(n.start:(n.start+dat$mu.f$pnum-1))];
		n.start <- n.start + dat$mu.f$pnum;
		par.Qq  <- snp.res[i, c(n.start:(n.start+dat$mu.f$pnum-1))];
		n.start <- n.start + dat$mu.f$pnum;
		par.qq  <- snp.res[i, c(n.start:(n.start+dat$mu.f$pnum-1))];

		str <- sprintf("Group=%d, Postion=%.2f", snp.res[i,1], snp.res[i,2]);
		Report.AddHeadline( str, level=2 );

		str <- paste(round( par.cov, 2 ), collapse=",");
		Report.AddHeadline( paste(dat$cov.f$name, str, sep=":"), level=3 );
		str <- paste(round( par.QQ, 2 ), collapse=",");
		Report.AddHeadline( paste("QQ", str, sep=":"), level=3 );
		str <- paste(round( par.Qq, 2 ), collapse=",");
		Report.AddHeadline( paste("Qq", str, sep=":"), level=3 );
		str <- paste(round( par.qq, 2 ), collapse=",");
		Report.AddHeadline( paste("qq", str, sep=":"), level=3 );

		c2 <- call("fpt.plot_com_curve", min(dat$times), max(dat$times), dat$mu.f$func, dat,
				QQ_par=par.QQ, Qq_par=par.Qq, qq_par=par.qq,
				xlab="Time", ylab="Model");

		Report.AddFigure( c2, "", c(3.5, 3 )*254, left.margin=0.1*254);
		Report.AddParaBreak();
	}

	Report.AddParaBreak();

	return(str);
}

proc_report_topn<-function( dat, topn )
{
	proc_report_snps( dat, topn)
}

proc_report_sig<-function( dat, sig.05, sig.01=NA )
{
	if ( is.null(sig.05) )
	{
		sErrMsg<- "Error: Not a result for fGWAS.";
		return(list(error=T, err.info=sErrMsg));
	}

	proc_report_snps( dat, sig.05);
}

fg_report<-function( obj.scan, pdf.file=NULL, options=list())
{
	if (is.null(pdf.file))
		pdf.file <- tempfile(pattern=obj.scan$obj.phe$params$file.phe.long, fileext=".PDF");

	Report.new( pdf.file, options );
	Report.title( "Functional GWAS Report", "fGWAS", "http://ccb.bjfu.edu.cn/" );
	Report.par( "dat.file", obj.scan$obj.phe$params$file.phe.long);

	Report.AddHeadline( "Data Summary", level=1 );
	proc_report_dat(obj.scan$obj.phe);

	if( !is.null(obj.scan$ret.gls) )
	{
		Report.AddHeadline( "GLS method", level=1 );
		proc_report_snpscan(obj.scan$obj.phe, obj.scan$ret.gls$result ) ;
	}

	if( !is.null(obj.scan$ret.fast) )
	{
		Report.AddHeadline( "FAST method", level=1 );
		proc_report_snpscan(obj.scan$obj.phe, obj.scan$ret.fast$result ) ;
	}

	if( !is.null(obj.scan$ret.fgwas) )
	{
		Report.AddHeadline( "fGWAS method", level=1 );
		proc_report_snpscan(obj.scan$obj.phe, obj.scan$ret.fgwas$result ) ;
	}

	#if( !is.null( r.perm) )
	#{
	#	Report.AddHeadline( "SNP Permutation", level=1 );
	#	proc_report_perm(dat, r.perm) ;
	#}

	#if( !is.null( r.perm) && !is.null( r.snpscan) )
	#{
	#	r.sig <- fg_detect_sig(obj.scan$dat, r.snpscan, r.perm)
	#
	#	Report.AddHeadline( "Significant SNP", level=1 );
	#	proc_report_sig( obj.scan$dat, r.sig$sig.05, r.sig$sig.01 ) ;
	#}

	#if( is.null( r.perm) && !is.null( r.snpscan) )
	#{
	#	Report.AddHeadline( "Top SNP", level=1 );
	#	proc_report_topn( obj.scan$obj.phe, r.snpscan) ;
	#}

	Report.Output( pdf.file );
}
