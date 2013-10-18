# Makefile to use knitr for package vignettes

# put all PDF targets here, separated by spaces
PDFS= CNVrd2.pdf

all: $(PDFS) 

clean:
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.toc
	rm -f figures/CNVrd2-*	CNVrd2-*
	rm -f CNVrd2.Rnw
	cp CNVrd2.Rnw.stub CNVrd2.Rnw
	
%.pdf: %.Rnw
	cp CNVrd2.Rnw.knitr CNVrd2.Rnw
	$(R_HOME)/bin/Rscript -e "library(knitr); knit2pdf('$*.Rnw')"

