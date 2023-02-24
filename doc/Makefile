
all: mpeg-pcc-ges-tm-sw-manual.pdf

MD_SRC=\
	README.about.md \
	README.build.md \
	README.options.md \
	README.usage.md \
	$(NULL)

TEX_SRC=\
	mpegdoc.cls \
	mpeg-pcc-ges-tm-sw-manual.tex \
	$(patsubst %.md,%.tex,$(MD_SRC)) \
	$(NULL)

mpeg-pcc-ges-tm-sw-manual.pdf: $(TEX_SRC)

LATEX:=$(shell which xelatex || which pdflatex || which latex)

# for bibtex
%.aux: %.tex
	 echo "x" | $(LATEX) $(LATEXFLAGS) $< || echo "Hack to skip error..."

# normally this has a rule to use bibtool to pull in mpeg bib entries
%.bib: %.tex
	touch $@

%.bbl: %.tex %.aux %.bib
	-bibtex $(*F)

ifneq ($(notdir $(LATEX)),latex)
# for direct pdf generation
%.pdf: %.tex %.bbl
	$(LATEX) $(LATEXFLAGS) $<
	$(LATEX) $(LATEXFLAGS) $<
	$(LATEX) $(LATEXFLAGS) $<
endif

%.pdf: %.dvi
	dvipdfm $*.dvi

%.dvi: %.tex %.bbl
	$(LATEX) $(LATEXFLAGS) $<
	$(LATEX) $(LATEXFLAGS) $<

# convert markdown text to latex for document generation
%.tex: %.md
	pandoc -f markdown -t latex -o $@ $<

clean:
	rm -f mpeg-pcc-ges-tm-sw-manual.log mpeg-pcc-ges-tm-sw-manual.out mpeg-pcc-ges-tm-sw-manual.pdf mpeg-pcc-ges-tm-sw-manual.toc mpeg-pcc-ges-tm-sw-manual.aux mpeg-pcc-ges-tm-sw-manual.bbl mpeg-pcc-ges-tm-sw-manual.bib README.about.tex README.build.tex README.options.tex README.usage.tex