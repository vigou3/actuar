TEXI2DVI = texi2dvi --pdf
RM = rm -rf

all: actuar.pdf credibility.pdf coverage.pdf lossdist.pdf risk.pdf simulation.pdf clean

actuar.pdf : actuar.tex
	$(TEXI2DVI) actuar.tex

credibility.pdf : credibility.tex
	$(TEXI2DVI) credibility.tex

coverage.pdf : coverage.tex
	$(TEXI2DVI) coverage.tex

lossdist.pdf : lossdist.tex
	$(TEXI2DVI) lossdist.tex

risk.pdf : risk.tex
	$(TEXI2DVI) risk.tex

simulation.pdf : simulation.tex
	$(TEXI2DVI) simulation.tex

clean:
	$(RM) *-[0-9][0-9][0-9].eps *-[0-9][0-9][0-9].pdf *.aux *.log  *.blg \
	*.bbl *.out *.tex *.rel *~ auto Rplots.ps Rplots.pdf
