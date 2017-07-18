SVG2PNG = inkscape --export-dpi=300 --export-area-drawing --export-png=
SVG2PDF = inkscape --export-area-drawing --export-pdf=
SVG_SOURCES = $(abspath $(wildcard *.svg))
PNG_IMAGES = $(SVG_SOURCES:%.svg=%.png)
PDF_IMAGES = $(SVG_SOURCES:%.svg=%.pdf)

all: pdf png

clean: pdf-clean png-clean

distclean: clean
	rm -f *~

png: $(PNG_IMAGES)

png-clean:
	rm -f *.png

pdf: $(PDF_IMAGES)

pdf-clean:
	rm -f *.pdf

.PHONY: all clean distclean pdf pdf-clean png png-clean

%.png: %.svg
	$(SVG2PNG)$@ $<

%.pdf: %.svg
	$(SVG2PDF)$@ $<
