# pombase-gocam-tool

[PomBase](https://www.pombase.org) tool for processing GO-CAM models.

## Compilation

`cargo build --release`

## Sub-commands

 - `all-genes`
 - `check-allowed-relations`
 - `connected-genes`
 - `cytoscape`
 - `cytoscape-model-connections`
 - `cytoscape-model-connections-with-rel-nodes`
 - `cytoscape-simple`
 - `cytoscape-simple-merged`
 - `detached-chemicals`
 - `detached-genes`
 - `find-holes`
 - `find-missing-evidence`
 - `find-obsolete-terms`
 - `genes-enabling-activities`
 - `gocam-py-parse-test`
 - `graph-viz-dot`
 - `joining-chemicals`
 - `make-chado-data`
 - `overlapping-nodes`
 - `per-model-stats`
 - `print-edges`
 - `print-individuals`
 - `print-nodes`
 - `print-tuples`
 - `print-unconnected-individuals`
 - `serialize`
 - `total-stats` - Find MF, BP or CC with missing evidence
 - `find-missing` - Find missing BP or CC
 - `write-annotation`

## Examples commands

```
curl https://raw.githubusercontent.com/pombase/pombase-gocam/refs/heads/main/tests/data/gomodel%3A66187e4700001744.json > gomodel:66187e4700001744.json
pombase-gocam-tool find-holes gomodel:66187e4700001744.json | less -S
pombase-gocam-tool stats gomodel:66187e4700001744.json | less -S
pombase-gocam-tool print-nodes gomodel:66187e4700001744.json | less -S

```

