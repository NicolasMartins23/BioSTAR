# BioSTAR

BioSTAR is an object-oriented Python library for biological sequence analysis.

## Architecture

The **BioSTAR package is the engine**. It contains the biological domain model and reusable analysis functionality and does not depend on a web framework or API layer.

- `BioStar/domain/` — DNA, RNA, Protein and other biological objects.
- `BioStar/data/` — biochemical constants and reference tables.
- `BioStar/analysis/` — sequence-analysis algorithms.
- `BioStar/io/` — input parsing such as FASTA.
- `BioStar/utils/` — small framework-independent helpers.

The package is intentionally independent of HTTP, FastAPI, databases and authentication.

## Future API

A future v2 API can be a separate FastAPI package/project:

```
FastAPI API -> BioSTAR engine -> biological domain
```

The API should depend on BioSTAR, never the other way around.

## Public API

```python
from BioStar import DNA, RNA, Protein, NucleicAcid, OpenReadFrame
```

Legacy module imports are retained as compatibility shims while new code should use the organized package structure.
