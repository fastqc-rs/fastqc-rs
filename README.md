# fastqc-rs
A fast quality control tool for FASTQ files written in rust. 

## Usage

```
fastqc-rs -q path/to/my_sequence.fastq > report.html
```

Arguments: 

| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| -q --fastq 	       |	-           |The path to the FASTQ file to use
| -k          | 5           |The length k of k-mers for k-mer counting
