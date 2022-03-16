# Use case diagrams

## High-level

```{mermaid}
flowchart LR
subgraph uc[pyHiM]
	f1([Preprocessing])
	f12([Identification])
	f7([Matching])
	f17([Postprocessing])
	f9([Visualization])
end
Image --- f1 & f12 & f7 & f17 & f9
uc -.- a2[Parallelization servers]
f12 --- a1[IA models]
```
## Preprocessing

```{mermaid}
flowchart LR
subgraph uc[Preprocessing]
	f1([Projection]) -->
    f2([Compute shift]) -->
	f12([Apply registration])
end
a0[Image 3D] --- f1
f12 --- a1[Identification]
```

## Identification

```{mermaid}
flowchart LR
subgraph uc[Identification]
	f1([Segment mask])
    f2([Detect spot])
end
a0[Masks 2D] --- f1
a1[Barcodes 2D or 3D] --- f2
f1 & f2 --- a2[Matching]
```

## Matching

```{mermaid}
flowchart LR
subgraph uc[Matching]
	f1([Match masks within masks])
    f2([Match barcodes within masks])
end
a0[Segmented masks] --- f1
a1[Segmented barcodes] --- f2
f1 & f2 --- a2[Postprocessing]
```

## Postprocessing

```{mermaid}
flowchart LR
subgraph uc[Postprocessing]
	f1([Filter data]) -->
    f2([Build PWD matrix])
end
a0[Localization matrix] --- f1
f2 --- a2[PWD matrix]
```

## Visualization

```{mermaid}
flowchart LR
subgraph uc[Visualization]
	f1([Plot])
    f2([Logging])
end
a0[Data] --- f1
a1[Runtime] --- f2
f1 --- a2[PNG image]
f2 --- a3[LOG or Markdown file]
```

