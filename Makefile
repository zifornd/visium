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
$(OUTPUT)/07-integrate-samples.html

#  $(OUTPUT)/06-cell-annotation.html \
#  $(OUTPUT)/07-merge-samples.html \
#  $(OUTPUT)/08-integrate-samples.html \
#  $(OUTPUT)/09-marker-detection.html \

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


# $(OUTPUT)/01-data-loading.html: $(WORKSPACE)/01-data-loading.qmd
# 	Rscript -e 'quarto::quarto_render("$<", output_file = "$@")'

# $(OUTPUT)/01-data-loading.html: $(WORKSPACE)/01-data-loading.Rmd
# 	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
# $(OUTPUT)/02-quality-control.html: $(WORKSPACE)/02-quality-control.Rmd $(OUTPUT)/01-data-loading.html
# 	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
# $(OUTPUT)/03-normalisation.html: $(WORKSPACE)/03-normalisation.Rmd $(OUTPUT)/02-quality-control.html
# 	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
# $(OUTPUT)/04-reduced-dimensions.html: $(WORKSPACE)/04-reduced-dimensions.Rmd $(OUTPUT)/03-normalisation.html
# 	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'

# $(OUTPUT)/05-clustering.html: $(WORKSPACE)/05-clustering.Rmd $(OUTPUT)/04-reduced-dimensions.html
# 	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
# $(OUTPUT)/06-cell-annotation.html: $(WORKSPACE)/06-cell-annotation.Rmd $(OUTPUT)/05-clustering.html
# 	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
# $(OUTPUT)/07-merge-samples.html: $(WORKSPACE)/07-merge-samples.Rmd $(OUTPUT)/06-cell-annotation.html
# 	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
# $(OUTPUT)/08-integrate-samples.html: $(WORKSPACE)/08-integrate-samples.Rmd $(OUTPUT)/06-cell-annotation.html
# 	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
# $(OUTPUT)/09-marker-detection.html: $(WORKSPACE)/09-marker-detection.Rmd $(OUTPUT)/08-integrate-samples.html
# 	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
	