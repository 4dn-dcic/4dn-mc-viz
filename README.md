# 4dn-mc-viz
Experimenting with multi-contact genomics file format

General h5 schema:
```console
/
 ├── chroms
 │   ├── length int32
 │   └── name |S64
 ├── bins
 │   ├── chrom int32
 │   ├── start int32
 │   └── end int32
 ├── clusters (attrs: [max])
 │   ├── bin int32
 │   ├── cluster_name int32
 │   └── cluster_length int32
 └── index
     └── bin_offset int32
```
