
## UMI tools in rust

**UMI** stands for Unique Molecular Identifier.

The goal is reproduce some functionalities from [`umi-tools`](https://umi-tools.readthedocs.io/en/latest/QUICK_START.html)
especially the `extract` command.

For example, convert the following read

```
@VH00666:90:AAAWVCCHV:1:1101:24026:1000
GTCAGTTATAGCGGGCGCGCAAAAAAAAAAAAAAAAAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTCCC
[...]
```

into:

```
@VH00666:90:AAAWVCCHV:1:1101:24026:1000_GTCAGT
GCGGGCGCGCAAAAAAAAAAAAAAAAAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTCCC
```

The UMI was **GTCAGT** and append to the read **name**, while being removed from the sequence along with the **TATA** linker.

## TODO

- Record where the TATA is found with stats.pos_linker (`HashMap<u32, u32>`)
- Read / write compressed FASTQ (https://github.com/wckdouglas/fq-filter-reads/)
- Use multithreading