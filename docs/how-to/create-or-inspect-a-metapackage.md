# Create or inspect a metapackage

Create a metapackage from SingleM packages.

```bash
singlem metapackage --singlem-packages marker-a.spkg marker-b.spkg \
  --metapackage markers.smpkg --threads 4
```

Inspect it:

```bash
singlem metapackage --metapackage markers.smpkg --describe
```
