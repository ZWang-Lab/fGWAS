#################################################################
#
# New Systems Mapping Application(SysMap1)
#
# plot utility
#
#    1) plot_qtl_map
#    2) plot_qtl_pos
#    3) plot_tiled_curves
#    4) plot_overlapping_curves
#	 5) plot_com_curve
#    6) plot_permutation
#    7) plot_correlation
#    8) plot_eval_ret
#
# History:
# 12/15/2011 Version 1.1
#
##################################################################

#--------------------------------------------------------------
# ftp.plot_tiled_curves
#
# Draw the growth curve of the individuals in some tiled plots,
# the plot's layout can be automatically calculated.
#
# Input:
#       dat : data object( phenos_table );
#--------------------------------------------------------------
fpt.plot_tiled_curves<-function( dat, rows=NA, cols=NA, max_curves = NULL, selected=NULL )
{
	max_log <- 0;
	min_log <- 0;

	if (is.null(max_curves))
		max_curves <- min(100,dat$n.obs);

	if ( is.null(selected) )
		selected <- c(1:max_curves)

	selected.greater <- which( selected > dat$n.obs);
	if (length (selected.greater) )
		selected <- selected[- selected.greater ];

	if ( is.na (rows) )
	{
		cols <- ceiling(  (length(selected))^0.5 );
		rows <- ceiling( (length(selected))/cols);
	}

	minv <- min(as.matrix(dat$phe[selected,]), na.rm=TRUE)*0.9;
	maxv <- max(as.matrix(dat$phe[selected,]), na.rm=TRUE)*1.1;

	px<-dat$times;
	xmax <- max(px, na.rm=TRUE);
	xmin <- min(px, na.rm=TRUE);

	py <- length ( dat$phe[1, ] );
	rx <- ceiling( py/8 );
	xunit <- c(1, 2, 3, 5, 8, 10, 20, 50, 100, 200, 500, 1000, 2000);
	rx <- xunit[min(which(xunit>rx))];
	if (rx>2000) rx<-as.integer(py/8);

	ry<- ceiling((maxv - minv)/8)
	yunit <- c(1, 2, 3, 5, 8, 10, 20, 50, 100, 200, 500, 1000, 2000);
	ry <- yunit[min(which(yunit>ry))];
	if (ry>2000) ry<-as.integer((maxv - minv)/8);

	p.width <- (xmax-xmin)*cols/9*10;
	p.height <- (maxv-minv)*rows/9*10;

	plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=T, xlab="",ylab="", xlim=c(0, p.width), ylim=c(0, p.height) );
	for (i in 0:(rows-1) )
	{
		for (j in 0:(cols-1) )
		{
			sub_fig<- c( 1/cols*j, 1/cols*(j+1), 1/rows*i,  1/rows*(i+1) )*0.98+0.02;
			sub_rc <- c( p.width * sub_fig[1], p.height* sub_fig[3],
			             p.width * sub_fig[2], p.height* sub_fig[4] );
			rect(sub_rc[1], sub_rc[2], sub_rc[3], sub_rc[4], border="black", lwd=1);
			sub_fc <- c( sub_rc[1]-xmin, sub_rc[2]-minv );

			if ( (i*cols+j+1) <= length(selected) )
			{
				idx <- i*cols+j+1;
				px<-dat$times;
				py<-as.numeric( dat$phe[selected[idx], ] );
				if (length(which(py==-1))>0)
				{
					px<-px[-(which(py==-1))];
					py<-py[-(which(py==-1))];
				}

				if (length(which(is.na(py)))>0)
				{
					px<-px[-(which(is.na(py)))];
					py<-py[-(which(is.na(py)))];
				}

				s01 <- smooth.spline(px, y=py);
				xx1 <- seq(xmin,xmax,by=0.2)

				#plot(c(1:2), c(1:2),xlim=c(xmin, xmax), ylim=c(minv, maxv),
				#	type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i", );

				if ( .RR( "curve_smooth", def=0 ) == 1)
					lines( predict(s01,xx1) + sub_fc, lwd=1)
				else
					lines(px+sub_fc[1], py+sub_fc[2], lwd=1);

				for (k in 1:length(px))
					points(px[k]+sub_fc[1], py[k]+sub_fc[2] ,type="o", pch=19, cex=0.5)

				h<-( maxv - minv )
				if (i==0)
					for (k in as.integer(xmin/rx):as.integer(xmax/rx))
					{
						if ( k*rx+sub_fc[1]+xmin <= sub_rc[1] || k*rx+sub_fc[1]+xmin >= sub_rc[3]) next;
						segments(k*rx+sub_fc[1]+xmin, minv+sub_fc[2], k*rx+sub_fc[1]+xmin, minv+h/40+sub_fc[2]);
					}

				if (j==0)
			  		for (k in as.integer(minv/ry):as.integer(maxv/ry) )
					{
						if ( k*ry+sub_fc[2] <= sub_rc[2] || k*ry+sub_fc[2] >= sub_rc[4]) next;
						segments(xmin+sub_fc[1], k*ry+sub_fc[2], xmin+(px-1)/40+sub_fc[1], k*ry+sub_fc[2]);
					}
			}
			else
			{
				#plot(c(1:2), c(1:2),xlim=c(1, 10), ylim=c(minv, maxv),
		    	#		type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i");
			}

			if (i*cols+j+1 <= dat$n.obs)
				text( 1+(xmax-1)/5 + sub_fc[1], h*8.5/10+minv + sub_fc[2], paste(i*cols+j+1, "", sep=""), font=4);

			if (j==0)
			{
		  		for (k in as.integer(minv/ry):(as.integer(maxv/ry)-1) )
				{
					if ((k*ry)<maxv && (k*ry)>minv && (k*ry-minv)/(maxv-minv)<0.95 && (k*ry-minv)/(maxv-minv)>0.01)
						text(0.5+sub_fc[1], k*ry+ sub_fc[2], k*ry, adj=c(1, 0.5),cex=0.5);
				}
			}

			if (i==0)
			{
				for (k in as.integer(xmin/rx):as.integer(xmax/rx))
				{
					if ((k*rx)/xmax<0.9 && (k*rx)/xmax>0.05)
						text(k*rx+ sub_fc[1]+xmin, 0.5+ sub_fc[2]+minv, (k*rx), adj=c(0.5,1),cex=0.5)
				}
			}
		}
	}
}

#--------------------------------------------------------------
# ftp.plot_overlapping_curves
#
# Draw the growth curve of the individuals in one plots, the curves
# will be overlapped.
#
# Input:
#       dat : data object( phenos_table );
#--------------------------------------------------------------
fpt.plot_overlapping_curves<-function( dat, options=list() )
{
	xmax <- max(dat$times, na.rm=TRUE);
	xmin <- min(dat$times, na.rm=TRUE);
	minv <- floor( min(as.matrix(dat$phe), na.rm=TRUE) )*0.9;
	maxv <- ceiling( max(as.matrix(dat$phe), na.rm=TRUE) )*1.1;

	rx <- xmax/16 ;
	xunit <- c(1, 2, 3, 5, 8, 10, 15, 20, 25, 30, 50, 75, 100, 150, 200, 250, 500, 750, 1000, 1500, 2000);
	rx <- xunit[min(which(xunit>rx))];
	if (rx>2000) rx<-as.integer(xmax/16);

	ry<- (maxv - minv)/25 ;
	yunit <- c(1, 2, 3, 5, 8, 10, 15, 20, 25, 30, 50, 75, 100, 150, 200, 250, 500, 750, 1000, 1500, 2000);
	ry <- yunit[min(which(yunit>ry))];
	if (ry>2000) ry<-as.integer((maxv - minv)/25);

	if (!is.null(options$fixed_yticks))
	{
		minv <-0;
		ry <- ceiling( (maxv - minv)/options$fixed_yticks );
		maxv <- minv + ry * options$fixed_yticks*1.1
	}

	op <- par(mar=c( 3.5, 3, 1, 1), mgp=c(1.75,0.5,0) );

	plot(c(1:2), c(1:2),xlim=c(xmin-abs(xmax-xmin)*1/100, xmax+abs(xmax-xmin)*1/100), ylim=c(minv, maxv ),
		type="n", main="", xlab="Time", ylab="Phenotypic Values", xaxt="n", yaxt="n" ); #xaxs="i", yaxs="i",

	#rect(xmin-abs(xmax-xmin)*1/100, minv, xmax+abs(xmax-xmin)*1/100, (maxv-minv)*0.99+minv, col = c(NA,"black"), lwd=1);

	for (i in 1:length(dat$phe[,1]) )
	{
		px<-dat$times;
		py<-as.numeric( dat$phe[i, ] );
		if (length(which(is.na(py)))>0)
		{
			px<-px[-(which(is.na(py)))];
			py<-py[-(which(is.na(py)))];
		}

		s01<-smooth.spline(px, y=py );
		xx1<-seq(xmin,xmax,by =0.2 )

		if ( .RR("curve_smooth",def=0) == 1)
		{
			pt<-predict(s01,xx1);
			lines(pt$x, pt$y, lwd=1, col="gray")
		}
		else
			lines(px, py, lwd=1, col="gray")


		for (k in 1:length(px))
			points(px[k], py[k], type="o", pch=19, cex=0.5,col="darkgray")
    }

	k <- c( as.integer(minv/ry):as.integer(maxv/ry) )
	axis(side=2, at=k*ry, col="black", col.axis="black", col.lab="black", font.axis=1, font= 1, cex.axis = ifelse(is.null(options$cex_axis), 1, options$cex_axis)  )

	k <- c( as.integer(xmin/rx):as.integer(xmax/rx) )
	axis(side=1, at=k*rx, col="black", col.axis="black", col.lab="black", font.axis=1, font= 1, cex.axis = ifelse(is.null(options$cex_axis), 1, options$cex_axis), padj=0.5 )

	points(dat$times, colMeans(dat$phe), type="o", pch=19, cex=0.5, col="red" );

	x10 <- seq(xmin,xmax,0.1);
	if (!is.null(dat$mu.f$par))
	{
		y10 <- dat$mu.f$func(dat$mu.f$par, x10);
		lines(x10, y10, lwd=2, col="black" , lty="22");
	}

	if (!any(is.na(dat$mu.f$lower)))
	{
		y10  <- dat$mu.f$func(dat$mu.f$lower, x10);
		lines(x10, y10, lwd=2, col="red" , lty="31");
	}

	if (!any(is.na(dat$mu.f$upper)))
	{
		y10  <- dat$mu.f$func(dat$mu.f$upper, x10);
		lines(x10, y10, lwd=2, col="red" , lty="31");
	}

	par(op);

}

#--------------------------------------------------------------
# ftp.plot_com_curve
#
# Draw the np curve for QQ, Qq, qq genotypes. the paramter of
# the pharmacology curve is E0, E50, Emax.
#
# Input:
#     nMesa : Mesuared times
#     nLong : total times
#     data  : data
#     QQ_Par: E0, E50 and Emax for QQ type
#     Qq_Par: E0, E50 and Emax for Qq type
#     qq_Par: E0, E50 and Emax for qq type
#--------------------------------------------------------------

fpt.plot_com_curve<-function( nMesa, nLong, f_curve_mu, dat = NULL, QQ_par=NULL, Qq_par=NULL, qq_par=NULL, simu_QQ=NULL, simu_Qq=NULL, simu_qq=NULL, xlab="Time", ylab="Phenotype" )
{
	if ( nMesa > nLong)
	 	nLong <- nMesa + 4;

	limit1 <- 10;
	if (!is.null(dat)) limit1 <- max( dat$phe, na.rm=TRUE );

	if (!is.null(QQ_par) && !any(is.na(QQ_par)) )
		limit1 <- max( limit1, f_curve_mu( QQ_par, 1:nMesa, options=list(tmin=QQ_par[1], tmax=nMesa ) ) );
	if (!is.null(Qq_par) && !any(is.na(Qq_par)) )
		limit1 <- max( limit1, f_curve_mu( Qq_par, 1:nMesa, options=list(tmin=Qq_par[1], tmax=nMesa ) ) );
	if (!is.null(qq_par) && !any(is.na(qq_par)) )
		limit1 <- max( limit1, f_curve_mu( qq_par, 1:nMesa, options=list(tmin=qq_par[1], tmax=nMesa ) ) );

	op <- par( mar=c(3,3,1,1), mgp=c(2,0.5,0) );
	plot(c(1:2), c(1:2),xlim=c( 1-1, nLong+1 ), ylim=c(0, limit1*1.1),
		type="n", main="",xlab=xlab, ylab=ylab, xaxt="s", yaxt="s"  );

	if (!is.null(dat))
	{
		tLen <- length( dat$phe[1,] );
		for (i in 1:length(dat$phe[,1]) )
		{
			#lines(c(1:tLen), dat$phe[i,], col=rgb(0.75,0.75,0.75), lwd=1);

			py <- as.numeric( dat$phe[i,] );
			px <- dat$times ;

			if (length(which(py==-1))>0)
			{
				px<-px[-(which(py==-1))];
				py<-py[-(which(py==-1))];
			}

			lines(px, py, col=rgb(0.75,0.75,0.75), lwd=1);
		}
	}

	x10 <- seq(0,nLong,0.1)	;
	p0  <- which(x10>=nMesa)[1];
	p1  <- length(x10);

	if (!is.null(QQ_par)  && !any(is.na(QQ_par)) )
	{
		y20 <- f_curve_mu( QQ_par, x10, options=list(tmin=QQ_par[1], tmax=nMesa) );
		lines(x10[1:p0],  y20[1:p0], lwd=2, col="darkgreen", lty="solid");
		lines(x10[p0:p1], y20[p0:p1], lwd=2, col="darkgreen", lty="solid");
	}

	if (!is.null(simu_QQ)  && !any(is.na(simu_QQ)) )
	{
		y20 <- f_curve_mu( simu_QQ, x10, options=list(tmin=simu_QQ[1], tmax=nMesa)) ;
		lines(x10[1:p0],  y20[1:p0], lwd=2, col="black", lty="dotted");
		lines(x10[p0:p1], y20[p0:p1], lwd=2, col="black", lty="dotted");
	}


	if (!is.null(Qq_par)  && !any(is.na(Qq_par)) )
	{
		y10 <- f_curve_mu( Qq_par, x10, options=list(tmin=Qq_par[1], tmax=nMesa) );
		lines(x10[1:p0],  y10[1:p0], lwd=2, col="blue" , lty="solid");
		lines(x10[p0:p1], y10[p0:p1], lwd=2, col="blue" , lty="solid");

		#t1<-as.integer(Qq_par[1]*nMesa);
		#segments(t1, 0, t1, limit1*1.1, lty="solid");
	}

	if (!is.null(simu_Qq) && !any(is.na(simu_Qq)))
	{
		y10 <- f_curve_mu( simu_Qq, x10, options=list(tmin=simu_Qq[1], tmax=nMesa) );
		lines(x10[1:p0],  y10[1:p0], lwd=2, col="black", lty="dotted");
		lines(x10[p0:p1], y10[p0:p1], lwd=2, col="black", lty="dotted");
	}

	if (!is.null(qq_par) && !any(is.na(qq_par)) )
	{
		y00 <- f_curve_mu( qq_par, x10, options=list(tmin=qq_par[1], tmax=nMesa)) ;
		lines(x10[1:p0],  y00[1:p0], lwd=2, col="red" , lty="solid");
		lines(x10[p0:p1], y00[p0:p1], lwd=2, col="red" , lty="solid");

		#t2<-as.integer(qq_par[1]*nMesa);
		#segments(t2, 0, t2, limit1*1.1, lty="solid");
	}


	if (!is.null(simu_qq) && !any(is.na(simu_qq)) )
	{
		y00 <- f_curve_mu( simu_qq, x10, options=list(tmin=simu_qq[1], tmax=nMesa));
		lines(x10[1:p0],  y00[1:p0], lwd=2, col="black" , lty="dotted");
		lines(x10[p0:p1], y00[p0:p1], lwd=2, col="black" , lty="dotted");
	}

	segments( nMesa, 0, nMesa, limit1*1.1, lty="dotted", col="black");

	plotchar <-c(20,21);
	sLegend<- c();
	colors <- c();

	if (!is.null(QQ_par) && !any(is.na(QQ_par)) )
	{
		sLegend<- c( sLegend, sprintf("(QQ%i):%3.2f,%3.2f,...",2, QQ_par[1], QQ_par[2]) );
		colors <- c( colors, "darkgreen");
	}
	if (!is.null(Qq_par) && !any(is.na(Qq_par)) )
	{
		sLegend<- c( sLegend, sprintf("(Qq%i):%3.2f,%3.2f,...",1, Qq_par[1], Qq_par[2]) );
		colors <- c( colors, "blue");
	}
	if (!is.null(qq_par) && !any(is.na(qq_par)) )
	{
		sLegend<- c( sLegend, sprintf("(qq%i):%3.2f,%3.2f,...",0, qq_par[1], qq_par[2]) );
		colors <- c( colors, "red");
	}

	if (length(sLegend)>0)
		legend(x="topright", y=NULL, sLegend, cex=0.6, col=colors,  pch=plotchar, lty=1, title="Legend")

	par(op);

}

#--------------------------------------------------------------
# ftp.plot_perm_curve
#
# Draw the permutation graph.
#
# Input:
#     pmdat :
#--------------------------------------------------------------
fpt.plot_perm_curve<-function( pmdat )
{
	op <- par(mar=c( 5, 5, 3, 1), mgp=c(3,1,0) );

	plot( 1, 1,xlim=c(-6,1), ylim=c(0, max(pmdat[,2]) ),
		type="n", xlab="p-value", ylab="Cutoff", xaxt="n", yaxt="s", xaxs="i", yaxs="i", main="Permutation result");

	lines( log10(pmdat[,1]), pmdat[,2], lty=1, col="blue")

	legend_left <- c();
	idx05 <- which(pmdat[,1]==0.05);
	str_legend <- c();
	legend_left <- c();

	if (length(idx05)>0)
	{
		segments(log10(0.05), 0, log10(0.05), pmdat[idx05[1], 2], lty=2 );
		str_legend <- sprintf("%-7.6f: %6.2f", 0.05, pmdat[idx05[1], 2]);
		legend_left <- c( str_legend );
	}

	sig_lines<-c( 0.01, 0.001, 0.0001, 0.00001, 0.000001);
	for (i in 1:length(sig_lines) )
	{
		idx01 <- which(pmdat[,1]==sig_lines[i]);
		if (length(idx01)>0)
		{
			segments(log10(sig_lines[i]), 0, log10(sig_lines[i]), pmdat[idx01[1], 2] , lty=2);
			str_legend <- sprintf("%-7.6f: %6.2f", sig_lines[i], pmdat[idx01[1], 2]);
			legend_left <- c( legend_left, str_legend );
		}
	}

	sig_lables<-c( 1, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001);
	axis(1, at=log10(sig_lables),labels=sig_lables, las=2);

	if (length(legend_left)>0)
	{
		str_max <- max(strwidth(str_legend));
		legend("topright", legend = legend_left,text.width =str_max , title = "Cutoffs");
	}

	par(op)
}


#--------------------------------------------------------------
# ftp.plot_correlation
#
# Draw the corelation graph.
#
# Input:
#     dat :
#     sTitle :
#     labels:
#     log10:
#     significance:
#--------------------------------------------------------------

fpt.plot_correlation<-function( dat, sTitle, labels=NULL, log10=TRUE, significance=0.05)
{
	nLenArm <- max( max(dat[,1]), max(dat[,2]) ) ;

	len <- length(dat[,1]);
	rm <- dat[,3 ];
	if (log10)
	{
		rm[ rm>1 ] <- 1;
		rm <- -log10 ( dat[,3 ] );
		rm[ rm<0 ] <- 0;
	}

	maxv <- max( rm );
	if (maxv <= -log10(0.05))
		maxv <- -log10(0.04);
	minv <- min( rm );
	minv <- 0;

	ramp <- colorRamp(c( "green", "blue" ));
	cols <- rgb ( ramp(seq(0, 1, length = 1000) ), maxColorValue=255 );
	mc <- array("#FFFFFF", dim = c( nLenArm, nLenArm));
	cls <- round( (rm-minv)/(maxv-minv)*1000 ) +1;
	cls[which(cls>1000)]<-1000;

	for (n in 1:len)
		mc[ dat[n,1], dat[n,2] ] <- cols[  cls[n] ];
	mc[is.na(mc)]<-"#FFFFFF";

	par(mar=c(2,2,2,2)+0.1);
	plot(c(0, nLenArm*1.45), c(0, nLenArm*1.45), type= "n", xlab="", ylab="", xaxt="n", yaxt="n");
	for (x in 1:nLenArm)
	for (y in 1:nLenArm)
	{
		ox <- (nLenArm-x+1);
		oy <- y;
		if (oy+ox<(nLenArm+1)) next;

		x0 <- -sqrt(2)/2*(ox-1) + sqrt(2)/2*(oy-1) + nLenArm*sqrt(2)/2;
		y0 <- 2*(nLenArm-1) - sqrt(2)/2*(oy-1+ox-1) - nLenArm*sqrt(2)/2*0.4;
		xs <- c(x0, x0+sqrt(2)/2, x0, x0-sqrt(2)/2, x0);
		ys <- c(y0, y0+sqrt(2)/2, y0+sqrt(2), y0+sqrt(2)/2, y0);

		if ( mc[x,y] != "#FFFFFF" )
			polygon(xs, ys, col=mc[x,y], border="gray", angle=-45)
		else
			polygon(xs, ys, col=mc[x,y], border="white", angle=-45);
	}


	l <- nLenArm * sqrt(2)/2;
	x0 <- nLenArm*sqrt(2)/4;
	for (i in 0:100)
	{
		rect(x0+i*(l/100), nLenArm*sqrt(2)/2*0.1,
		     x0+(i+1)*l/100, nLenArm*sqrt(2)/2*0.2,
		     col=cols[i*10+1], border=cols[i*10+1])
	}

	text( x0+0* (l/100), nLenArm*sqrt(2)/2*0.1-strheight("1"), "1",cex=0.5 );
	x50 <- round( (-log10(0.5) - minv)/(maxv-minv)*100 ) +1;
	if (x50<100)
	{
		segments(x0+x50*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x50*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x50*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".5",cex=0.5 );
	}
	x005 <- round( (-log10(0.05) - minv)/(maxv-minv)*100 ) +1;
	if (x005<100)
	{
		segments(x0+x005*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x005*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x005*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".05",cex=0.5 );
	}
	x001 <- round( (-log10(0.01) - minv)/(maxv-minv)*100 ) +1;
	if (x001)
	{
		segments(x0+x001*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x001*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x001*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".01",cex=0.5 );
	}
	x0001 <- round( (-log10(0.001) - minv)/(maxv-minv)*100 ) +1;
	if (x0001<100)
	{
		segments(x0+x0001*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x0001*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x0001*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".001",cex=0.5 );
	}
	x00001 <- round( (-log10(0.0001) - minv)/(maxv-minv)*100 ) +1;
	if (x00001<100)
	{
		segments(x0+x00001*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x00001*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x00001*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".0001",cex=0.5 );
	}

	for (n in 1:len )
	{
		if ( dat[n,3] <= significance )
		{
			ox <- (nLenArm - dat[n,1] + 1);
			oy <- dat[n,2];
			if (oy+ox<(nLenArm+1)) next;

			x0 <- -sqrt(2)/2*(ox-1) + sqrt(2)/2*(oy-1) + nLenArm*sqrt(2)/2;
			y0 <- 2*(nLenArm-1) - sqrt(2)/2*(oy-1+ox-1) - nLenArm*sqrt(2)/2*0.4;
			if (nLenArm<10)
			{
				str <- sprintf("%5.3f", dat[n,3]);
				text( x0, y0+sqrt(2)/2, str, col="black", srt=-45);
			}
			else
			{
				points( x0, y0+sqrt(2)/2, pch=3,cex=0.8, col="yellow");
			}
		}
	}

	if (!is.null( labels ))
	{
		for (n in 1:length(labels) )
		{
			x <- (n-1)*sqrt(2)+sqrt(2)/2;
			y <- sqrt(2)/2*nLenArm * (1.4) ;
			text(x, y-0.5, labels[n], adj=c(1, 0), srt=-90);
		}
	}


	title( sTitle );
}

fpt.plot_manhattan<-function( res, p.05=NA, p.01=NA, sig.level=NA, p.sig=NA, map.title="" )
{
	op <- par( mgp = c( 1.5, 0.5,0 ), mar=c( 4, 4, 1.5, 1.5) );

	y.max <- round(max(-log10(res[,3]), na.rm=T))+1;

	#ylab=expression(-log[10](italic(p))),
	plot( 1,1, type="n", xlab="SNP", ylab="-log10(pvalue)",  cex.axis=0.7, xlim=c(1, dim(res)[1]), ylim=c(0, y.max ) );

	if(!is.na(p.05))
	{
		abline( h=-log10(p.05), col="gray", lwd=1, lty="dashed");
		text( x=0, -log10(p.05) + 0.1, "p=0.05", cex=0.6, srt=90, adj=c(0.5, -1));
	}

	if(!is.na(p.01))
	{	abline( h=-log10(p.01), col="gray", lwd=1, lty="dashed");
		text( x=0, -log10(p.01) + 0.1, "p=0.01", cex=0.6, srt=90, adj=c(0.5, -1));
	}

	if(!is.na(p.sig))
	{	abline( h= -log10(p.sig), col="gray", lwd=1, lty="dashed");
		text( x=0, -log10(p.sig) + 0.1, paste("p=",sig.level, sep=""), cex=0.6, srt=90, adj=c(0.5, -1));
	}

	cols <- c( "green","black",  "orange",  "red", "yellow", "blue", "purple");
	points( 1:NROW(res), -log10(res[,3]), pch=20, col=cols[ (res[,1]%%7+1)], cex=0.5 );

	if (map.title !="" )
		title( map.title );

	par(op);
}

fpt.plot_curve <- function( obj.phe, par_h0=NULL, par_h1=NULL, snp.vec=NULL, include.rawdata=TRUE, include.meanvector=TRUE, extra=list(), ...)
{
	xlines<-function(x, y, col="black", lwd=1.0, lty="solid")
	{
		y.na <- which(is.na(y));
		if(length(y.na)>0)
		{
			x <- x[-y.na];
			y <- y[-y.na];
		}
		lines(x, y, col=col, lwd=lwd, lty=lty);
	}

	ylines<-function(x, y0, y1, col="black", lwd=1.0, lty="solid")
	{
		y.na <- which(is.na(y0));
		if(length(y.na)>0)
		{
			x <- x[-y.na];
			y0 <- y0[-y.na];
			y1<- y1[-y.na]
		}

		for(i in 1:length(x))
		{
			segments(x, y0 + 0.5*y1, x, y0 - 0.5*y1, col=col, lwd=lwd, lty=lty);
		}
	}

	get_mean_vector<-function(pheY0, pheT0)
	{
		pheT_vec <- c(pheT0)
		pheY_vec <- c(pheY0)
		ti <- unique(pheT_vec);
		ti <- sort( ti [ !is.na(ti) ] );
		y.mu <- rep(NA, length(ti));
		y.count <- rep(NA, length(ti));
		for(i in 1:length(ti) )
		{
			y.mu[i] <-  mean(pheY_vec[pheT_vec==ti[i]], na.rm=T);
			y.count[i] <- sum(pheT_vec==ti[i], na.rm=T)
		}
		return(list(time=ti, value=y.mu, count=y.count))
	}

	get_est_vector<-function( pheY0, pheX0, pheT0, par_X, par_curve, ti)
	{
		if( !is.null(par_X) && length(par_X)>0 && !is.null(pheX0) )
			X.cov <- matrix( pheX0 %*% par_X, ncol=1, byrow=T) %*% matrix(rep(1,NROW(ti)), nrow=1, byrow=T)
		else
			X.cov <- matrix(rep(0,NROW(ti)), nrow=1, byrow=T);

		y.est <- get_curve( obj.phe$obj.curve, par_curve, ti, options=options )  + colMeans( X.cov, na.rm=T);

		return(list(time=ti, value=y.est))
	}

	pheY <- obj.phe$pheY;
	pheX <- obj.phe$pheX;
	pheT <- obj.phe$pheT;
    options <- list(min.time=min(pheT, na.rm=T), max.time=max(pheT, na.rm=T));
	h0.mean <- h0.est <- h1A.mean <- h1A.est <- h1B.mean<- h1B.est <- h1C.mean<- h1C.est <- NULL

	if( obj.phe$intercept )
		if (is.null(pheX))
			pheX <- matrix(1, nrow=NROW(pheY), ncol=1 )
		else
			pheX <- cbind(1, pheX);

	par_x_num <- ifelse( is.null( pheX ), 0 , NCOL( pheX ) )
	par_curve_num <- get_param_info(obj.phe$obj.curve, pheT )$count;

	ti <- seq( min( pheT, na.rm=T ), max( pheT, na.rm=T ), length=100 )
	if(!is.null(par_h0))
	{
		h0.mean<- get_mean_vector( pheY, pheT ) ;
		h0.est <- get_est_vector( pheY, pheX, pheT, if(is.null(pheX)) NULL else par_h0[1:par_x_num], par_h0[(par_x_num+1):(par_x_num+par_curve_num)], ti);
	}

	snp0.idx <- snp1.idx <- snp2.idx <- NULL;
	if(!is.null(snp.vec))
	{
		snp0.idx <- which(snp.vec==0);
		snp1.idx <- which(snp.vec==1);
		snp2.idx <- which(snp.vec==2);
	}

	if(!is.null(par_h1))
	{
		if(length(snp0.idx)>0)
		{
			h1A.mean<- get_mean_vector( pheY[snp0.idx,,drop=F], pheT[snp0.idx,,drop=F] ) ;
			h1A.est <- get_est_vector( pheY[snp0.idx,,drop=F], if(!is.null(pheX)) pheX[snp0.idx,,drop=F] else NULL, pheT[snp0.idx,,drop=F], par_h1[1:par_x_num], par_h1[(par_x_num+1):(par_x_num+par_curve_num)], ti);
		}

		if(length(snp1.idx)>0)
		{
			h1B.mean<- get_mean_vector( pheY[snp1.idx,,drop=F], pheT[snp1.idx,,drop=F] ) ;
			h1B.est <- get_est_vector( pheY[snp1.idx,,drop=F], if(!is.null(pheX)) pheX[snp1.idx,,drop=F] else NULL, pheT[snp1.idx,,drop=F], par_h1[1:par_x_num], par_h1[(par_x_num+1+par_curve_num):(par_x_num+par_curve_num*2)], ti);
		}

		if(length(snp2.idx)>0)
		{
			h1C.mean<- get_mean_vector( pheY[snp2.idx,,drop=F], pheT[snp2.idx,,drop=F] ) ;
			h1C.est <- get_est_vector( pheY[snp2.idx,,drop=F], if(!is.null(pheX)) pheX[snp2.idx,,drop=F] else NULL, pheT[snp2.idx,,drop=F], par_h1[1:par_x_num], par_h1[(par_x_num+1+par_curve_num*2):(par_x_num+par_curve_num*3)], ti);
		}
	}

	gen.col <- c("pink", "khaki", "skyblue");
	dots <- list(...)    ;
	if (is.null(dots$xlim)) dots$xlim = c(min(pheT, na.rm=T), max(pheT, na.rm=T));
	if (is.null(dots$ylim)) dots$ylim = c(min(pheY, na.rm=T), max(pheY, na.rm=T));
	if (is.null(dots$xlab)) dots$xlab = "Time";
	if (is.null(dots$ylab)) dots$ylab = "Phenotype";

	plot(1,1, type="n", xlim=dots$xlim, ylim=dots$ylim, xlab=dots$xlab, ylab=dots$ylab, cex=0.75);
	if(include.rawdata)
		for(i in 1:NROW(pheY))
			lines(pheT[i,], pheY[i,], col=ifelse( is.null(snp.vec), "gray", gen.col[snp.vec[i]+1]), lwd=0.2);

	if(!is.null(par_h0))
	{
		if(include.meanvector)
		{
			xlines(h0.mean$time, h0.mean$value, col="black", lwd=1, lty=11)
			points(h0.mean$time,  h0.mean$value,  col="black", pch=21, cex=h0.mean$count/NROW(pheY)*1.5);
		}
		xlines(h0.est$time,  h0.est$value,  col="black", lwd=1.5)
	}

	if(!is.null(par_h1))
	{
		if (length(snp0.idx)>0 && include.meanvector) xlines(h1A.mean$time, h1A.mean$value, col="red", lwd=0.8, lty=11)
		if (length(snp0.idx)>0) xlines(h1A.est$time, h1A.est$value, col="red", lwd=0.8)

		if (length(snp1.idx)>0 && include.meanvector) xlines(h1B.mean$time, h1B.mean$value, col="darkgreen", lwd=0.8, lty=11)
		if (length(snp1.idx)>0) xlines(h1B.est$time, h1B.est$value, col="darkgreen", lwd=0.8)

		if (length(snp2.idx)>0 && include.meanvector) xlines(h1C.mean$time, h1C.mean$value, col="blue", lwd=1, lty=11)
		if (length(snp2.idx)>0) xlines(h1C.est$time, h1C.est$value, col="blue", lwd=0.8)
	}

	if(!is.null(snp.vec))
	{
		legend("topleft", legend=c(paste("G0=",length(snp0.idx)), paste("G1=",length(snp1.idx)), paste("G2=",length(snp2.idx)) ),
			col=c("red", "darkgreen", "blue" ),lty=c("solid","solid","solid"), cex=0.75  )
		title(main = paste( extra$METHOD, "SNP=", extra$NAME, "LR2=", round(extra$LR2, 1), sep=" "), cex=0.75)

		text(median(c(pheT), na.rm=T), min(pheY, na.rm=T), paste( "MAF=", round(extra$MAF, 3), "NMISS=", extra$NMISS,sep=" "), cex=0.75, adj=c(0.5, 0.5) );
	}

	return
}