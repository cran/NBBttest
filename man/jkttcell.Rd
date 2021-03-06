\name{jkttcell}
\alias{jkttcell}
\docType{data}
\title{
Jurkat T-cell transcritomic data
}
\description{
The data are transcriptomic count data of 
\emph{RNA} reads generated by next 
generation sequencing from Jurkat T-cells.
}
\usage{data("jkttcell")}
\format{
  A data frame with 13409 observations on the 
  following 13 variables.
  \describe{
    \item{\code{tagid}}{: a numeric vector}
    \item{\code{geneid}}{: a numeric vector}
    \item{\code{name}}{: a string vector}
    \item{\code{chr}}{: a string vector} 
    \item{\code{strand}}{: a charactor vector}
    \item{\code{pos}}{: a numeric vector}
    \item{\code{anno}}{: a string vector}
    \item{\code{Jurk.NS.A}}{: a numeric vector}
    \item{\code{Jurk.NS.B}}{: a numeric vector}
    \item{\code{Jurk.NS.C}}{: a numeric vector}
    \item{\code{Jurk.48h.A}}{: a numeric vector}
    \item{\code{Jurk.48h.B}}{: a numeric vector}
    \item{\code{Jurk.48h.C}}{: a numeric vector}
  }
}
\details{
The data are count data generated by next generation 
sequencing from Jurkat T-cells. The T-cells were 
treated by resting and stimulating with \emph{CD3/CD28} 
for 48 hours. The data have 7 columns for the information 
of \emph{poly(A)} sites: tagid, geneid, gene name, 
chromosome, strand, \emph{poly(A)} site position, 
\emph{poly(A)} site annotation and 6 columns 
for count data: \code{Jurk.NS.A},	\code{Jurk.NS.B},	
\code{Jurk.NS.C},	\code{Jurk.48h.A},	\code{Jurk.48h.B},	
\code{Jurk.48h.C} where NS means normal state or no 
sitmulation and 48h means 48 hours after \emph{CD3/CD28} 
stimulatuin of T-cells. 13409 \emph{RNA} isoforms 
were detected to have alternative 
\emph{poly(A)} sites.
}
\source{
Real transcriptomic count data
}
\value{ID, information, count data of RNA reads}
\references{
Yuan-De Tan Anita M. Chandler, Arindam Chaudhury, 
and Joel R. Neilson(2015) A Powerful Statistical 
Approach for Large-scale Differential Transcription 
Analysis. \emph{Plos One}. 
DOI: 10.1371/journal.pone.0123658. 
}
\examples{
data(jkttcell)
}
\keyword{datasets}
