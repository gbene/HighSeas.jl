# Saving and reading

Here are objects and functions related to saving and reading simulation results.

## Saving

Saving data is managed automatically when solving. The user can choose what to save using the appropriate AbstractSavers. The data will be saved as in JLD2 format
for all savers except for the SnapshotSaver which will be saved as png or svg depending on the backend (GPU will save as .png, CPU .svg)

```@docs
StepSaver
```

```@docs
CatalogSaver
```

```@docs
SamplerSaver
```

```@docs
SnaptshotSaver
```




## Reading

```@docs
loadData
```

```@docs
loadSSH
```

```@docs
readSheet
```

```@docs
HighSeas.LoadedStep
```

```@docs
HighSeas.LoadedSamplers
```

