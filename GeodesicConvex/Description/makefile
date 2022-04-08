# makefile for IPOL articles
# to be used with the IPOL LaTeX class verson 0.4

DOC	= article
TEX	= $(wildcard *.tex)
BIB	= $(wildcard *.bib)
LATEX	= pdflatex -interaction nonstopmode -halt-on-error
BIBTEX	= bibtex
DPILR	= 96
PDFLR	= ps2pdf -dCompatibilityLevel=1.4 \
		-dEmbedAllFonts=true \
		-dAutoRotatePages=/None \
		-dPDFSETTINGS=/screen \
		-dColorConversionStrategy=/LeaveColorUnchanged \
		-dColorImageResolution=$(DPILR) \
		-dGrayImageResolution=$(DPILR) \
		-dMonoImageResolution=$(DPILR)

default :	pdf
pdf:	$(DOC).pdf
all:	$(DOC).pdf $(DOC)_lr.pdf

# LaTeX->PDF compilation
%.pdf	: %.tex $(TEX) $(BIB)
	-$(LATEX) -draftmode $*.tex >/dev/null 2>&1
	-$(BIBTEX) $*.aux
	-$(LATEX) -draftmode $*.tex >/dev/null 2>&1
	$(LATEX) $*.tex

# PDF filtering
%_lr.pdf	: %.pdf %.info
	$(PDFLR) $< $@
	pdftk $@ multibackground lr_header.pdf output $@.tmp \
		&& mv -f $@.tmp $@
	pdftk $@ update_info $*.info output $@.tmp \
		&& mv -f $@.tmp $@
	qpdf --linearize --stream-data=compress $@ $@.tmp \
		&& mv -f $@.tmp $@

# metadata
%.info	: %.pdf
	pdftk $< dump_data output $@

# cleanup
clean	:
	$(RM) $(DOC).aux $(DOC).log $(DOC).out
	$(RM) $(DOC).blg $(DOC).bbl $(DOC).toc
	$(RM) $(DOC).info
distclean	: clean
	$(RM) $(DOC).pdf $(DOC)_lr.pdf

.PHONY	: pdf all clean distclean
