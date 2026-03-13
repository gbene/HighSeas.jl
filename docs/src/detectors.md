# Detectors

Here are all the detectors that can be used to detect if an event is happening or not. When a detector is used, the simulation state is saved in the experiments' output folder.

## Simple detectors

These are simple detectors that either do nothing (necessary for code homogeneity) or very simple things (i.e. detect the event but do nothing with the data) 

```@docs
EmptyDetector
```

```@docs
SimpleDetector
```


## Complex detectors

These are more advanced detectors that detect events and do something (i.e. define a catalog entry) 

```@docs
CatalogDetector
```

