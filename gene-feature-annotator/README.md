# Gene Feature Annotator

This is a Docker container for a gene feature annotator tool. The container will contain all the necessary dependencies and configurations, making it easy to deploy and run the tool on any machine.

## Installation

```sh
git clone https://github.com/theopavlove/gene-feature-annotator.git

cd gene-feature-annotator

docker build -t gene-feature-annotator .
```

## How to run

### From shell

```sh
docker run -v /path/to/data/:/data \
    gene-feature-annotator \
    /data/input/input.bed \
    /data/input/input.gtf \
    /data/output/output.tsv
```

### From Python

TODO

<!-- ```py
import docker

client = docker.from_env()
client.containers.run("gene-feature-annotator", "sleep infinity", detach=True)
``` -->

## Output format

The resulting table contains the following columns:

- `geneID` -- Gene ID.
- `seqnames` -- Peak chromosome / sequence name.
- `start` -- Peak Start coordinate.
- `end` -- Peak End coordinate.
- `width` -- Peak length.
- `strand` -- Peak strand.
- `V4`, `V5`... -- Remaining columns from the initial bed-file.
- `annotation` -- Annotation type for the given peak.
- `geneChr` -- Gene chromosome.
- `geneStart` -- Gene Start coordinate.
- `geneEnd` -- Gene End coordinate.
- `geneLength` -- Gene length.
- `geneStrand` -- Gene strand.
- `transcriptId` -- Transcript ID.
- `distanceToTSS` -- Distance from the peak to the TSS of the associated gene.
- `gene_name` -- Gene name column from the `.gtf` file. 
- `gene_type` -- Gene type column from the `.gtf` file.

## Example

```sh
docker run -v $(pwd)/data:/data \
    gene-feature-annotator \
    /data/input/mm10.Z-DNA.bed \
    /data/input/mm10.gencode.basic.annotation.gtf.gz \
    /data/output/mm10.Z-DNA.annotation.tsv
```

## Future plans

- [ ] Add support for multiple bed files input
- [ ] Speed up the program (if possible)