FROM bioconductor/bioconductor_docker

# Install required packages
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("GenomicFeatures")'
RUN R -e 'BiocManager::install("ChIPseeker")'
RUN R -e 'BiocManager::install("rtracklayer")'

# Copy the R script to the container
COPY src/gene_feature_annotator.R /usr/local/bin/gene_feature_annotator.R
RUN chmod +x /usr/local/bin/gene_feature_annotator.R

# Set the entry point
ENTRYPOINT ["/usr/local/bin/gene_feature_annotator.R"]
