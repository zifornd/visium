WORKSPACE = $(realpath analysis)
OUTPUT = $(realpath output)

all: index.html \
     about.html \
     $(OUTPUT)/01-data-loading.html \
     $(OUTPUT)/02-quality-control.html \
     $(OUTPUT)/03-normalisation.html \
     $(OUTPUT)/04-reduced-dimensions.html \
     $(OUTPUT)/05-clustering.html \
     $(OUTPUT)/06-merge-samples.html \
     $(OUTPUT)/07-integrate-samples.html \
     $(OUTPUT)/08-marker-detection.html \
	 $(OUTPUT)/09-fea-annotation.html \
     $(OUTPUT)/10-spatial-features.html \
	 $(OUTPUT)/11-cell-atlas.html


index.html: index.qmd
	quarto render $<

about.html: about.qmd
	quarto render $<

$(OUTPUT)/01-data-loading.html: $(WORKSPACE)/01-data-loading.qmd
	quarto render $<

$(OUTPUT)/02-quality-control.html: $(WORKSPACE)/02-quality-control.qmd $(OUTPUT)/01-data-loading.html
	quarto render $<

$(OUTPUT)/03-normalisation.html: $(WORKSPACE)/03-normalisation.qmd $(OUTPUT)/02-quality-control.html
	quarto render $<

$(OUTPUT)/04-reduced-dimensions.html: $(WORKSPACE)/04-reduced-dimensions.qmd $(OUTPUT)/03-normalisation.html
	quarto render $<

$(OUTPUT)/05-clustering.html: $(WORKSPACE)/05-clustering.qmd $(OUTPUT)/04-reduced-dimensions.html
	quarto render $<

$(OUTPUT)/06-merge-samples.html: $(WORKSPACE)/06-merge-samples.qmd $(OUTPUT)/05-clustering.html
	quarto render $<

$(OUTPUT)/07-integrate-samples.html: $(WORKSPACE)/07-integrate-samples.qmd $(OUTPUT)/06-merge-samples.html
	quarto render $<

$(OUTPUT)/08-marker-detection.html: $(WORKSPACE)/08-marker-detection.qmd $(OUTPUT)/07-integrate-samples.html
	quarto render $<

$(OUTPUT)/09-fea-annotation.html: $(WORKSPACE)/09-fea-annotation.qmd $(OUTPUT)/08-marker-detection.html
	quarto render $<

$(OUTPUT)/10-spatial-features.html: $(WORKSPACE)/10-spatial-features.qmd $(OUTPUT)/09-fea-annotation.html
	quarto render $<

$(OUTPUT)/11-cell-atlas.html: $(WORKSPACE)/11-cell-atlas.qmd $(OUTPUT)/10-spatial-features.html
	quarto render $<
	
	