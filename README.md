# vcf_extractor
Extracts VCF fields to a tab-delimited text file.

## Running extract_vcf

```extract_vcf``` is written in C, and uses HTSLib functions. HTSLib is required to compile and execute.

If an HTSLib module is already available on your system, it should be loaded before compiling or running the executable:

```
module spider htslib
module load bioinf/htslib/1.9
```

If HTSLib is unavailable on your system, you can [download HTSLib here](https://www.htslib.org/download/) and [install from source](https://github.com/samtools/htslib/blob/develop/INSTALL). If ```extract_vcf``` can't find HTSLib at runtime, execution will fail with ```./extract_vcf: error while loading shared libraries: libhts.so.2: cannot open shared object file: No such file or directory```. I compile and run with HTSLib 1.9; other versions of HTSLib have not been tested and may cause unexpected behaviour.

```extract_vcf``` requires an input VCF file and output tab-delimited text file as command line arguments or redirection from STDIN/to STDOUT using -. VCF files can be compressed or uncompressed, indexed or unindexed.

By default, ```extract_vcf``` will output the SAMPLE, CHROM, POS, REF, and ALT fields, where SAMPLE is the sample name in the VCF header. However, it accepts several additional options to extract other VCF fields:
* --id: Extracts the ID field with column name ID
* --info: Extracts the specified INFO fields with column names set as INFO tags
* --format: Extracts the specified FORMAT fields with column names set as FORMAT tags
* --sample-names: Add sample names from command line input rather than using the sample names present in the VCF header; the number of names specified must match the number of samples within the VCF file

Additionally, ```extract_vcf``` accepts options designed to transform VCF fields parsed by ```bcftools +split-vep``` or natively comma-separated by transcript into long format:
* --split-fields (and --delimiter): INFO fields specified in --split-fields are split on --delimiter (by default ','), with each element appearing on a new line

```
./extract_vcf --info vep_csq in.vcf out.tsv
chr1 654823 A T missense_variant,missense_variant,synonymous_variant

./extract_vcf --info vep_csq --split-fields vep_csq in.vcf out.tsv

chr1 654823 A T missense_variant
chr1 654823 A T missense_variant
chr1 654823 A T synonymous_variant
```

To test whether ```extract_vcf``` is working correctly, you can run the following command:

```
./extract_vcf --id --info AC,AF,SOR --format GT,DP test/test.vcf test/test.tsv
```

## Compiling extract_vcf

If you need to compile ```extract_vcf``` from source, first install or load HTSLib. You can then compile the executable using the following command:

```gcc -O3 -o extract_vcf extract_vcf.c -lhts```

Compilation was done using gcc 10.2.0 on Linux.
