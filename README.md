# fastqc-rs
A fast quality control tool for FASTQ files written in rust.

Available statistics are:
- Read length
- Sequence quality score
- Sequence quality per base
- Sequence content per base
- k-mer content
- GC content

## Usage

```
fastqc-rs -q path/to/my_sequence.fastq > report.html
```

Arguments: 

| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| -q --fastq 	       |	-           |The path to the FASTQ file to use
| -k          | 5           |The length k of k-mers for k-mer counting
