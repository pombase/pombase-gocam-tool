# pombase-gocam-tool

[PomBase](https://www.pombase.org) tool for processing GO-CAM models.

## Compilation

`cargo build --release`

## Sub-commands

 - `find-holes`
 - `stats`
 
## Examples commands

```
curl https://raw.githubusercontent.com/pombase/pombase-gocam/refs/heads/main/tests/data/gomodel%3A66187e4700001744.json > gomodel:66187e4700001744.json
pombase-gocam-tool find-holes gomodel:66187e4700001744.json | less -S
pombase-gocam-tool stats gomodel:66187e4700001744.json | less -S
```

