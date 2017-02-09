\name{fg.report}
\alias{fg.report}
\title{
Making report for fGWAS data analysis
}
\description{
}
\usage{
fg.report(  fgwas.scan.obj, file.pdf, options=list()  )
}
\arguments{
  \item{fgwas.scan.obj}{}
  \item{file.pdf}{}
  \item{options}{}
}
\details{
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fgwas.scan.obj, file.pdf, options = list()) 
{
    ret <- fg_report(fgwas.scan.obj, file.pdf, options)
    return
  }
}
\keyword{ fGWAS }
\keyword{ Report }
