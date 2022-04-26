WORKSPACE = /sbgenomics/workspace
OUTPUT = /sbgenomics/output-files

all: $(OUTPUT)/01-data-loading.html \
		 $(OUTPUT)/02-quality-control.html \
		 $(OUTPUT)/03-normalisation.html \
		 $(OUTPUT)/04-reduced-dimensions.html \
		 $(OUTPUT)/05-clustering.html \
		 $(OUTPUT)/06-cell-annotation.html \
		 $(OUTPUT)/07-merge-samples.html \
		 $(OUTPUT)/08-integrate-samples.html \
		 $(OUTPUT)/09-marker-detection.html \


$(OUTPUT)/01-data-loading.html: $(WORKSPACE)/01-data-loading.Rmd
	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
$(OUTPUT)/02-quality-control.html: $(WORKSPACE)/02-quality-control.Rmd $(OUTPUT)/01-data-loading.html
	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
$(OUTPUT)/03-normalisation.html: $(WORKSPACE)/03-normalisation.Rmd $(OUTPUT)/02-quality-control.html
	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
$(OUTPUT)/04-reduced-dimensions.html: $(WORKSPACE)/04-reduced-dimensions.Rmd $(OUTPUT)/03-normalisation.html
	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'

$(OUTPUT)/05-clustering.html: $(WORKSPACE)/05-clustering.Rmd $(OUTPUT)/04-reduced-dimensions.html
	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
$(OUTPUT)/06-cell-annotation.html: $(WORKSPACE)/06-cell-annotation.Rmd $(OUTPUT)/05-clustering.html
	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
$(OUTPUT)/07-merge-samples.html: $(WORKSPACE)/07-merge-samples.Rmd $(OUTPUT)/06-cell-annotation.html
	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
$(OUTPUT)/08-integrate-samples.html: $(WORKSPACE)/08-integrate-samples.Rmd $(OUTPUT)/06-cell-annotation.html
	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
$(OUTPUT)/09-marker-detection.html: $(WORKSPACE)/09-marker-detection.Rmd $(OUTPUT)/08-integrate-samples.html
	Rscript -e 'rmarkdown::render("$<", output_file = "$@")'
	
	