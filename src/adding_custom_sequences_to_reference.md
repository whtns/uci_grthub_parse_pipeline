


[Skip to main content](#main-content)
[![Parse Biosciences Home](Adding%20Custom%20Sequences%20and%20Gene%20Annotation%20File%20Formatting%20%E2%80%93%20Support%20Suite%20-%20Parse%20Biosciences_files/logo-white.png)](https://www.parsebiosciences.com/)

[Products](https://www.parsebiosciences.com/products/)
[Technology](https://www.parsebiosciences.com/technology/)
[Company](https://www.parsebiosciences.com/company/)
[Resources](https://www.parsebiosciences.com/resources/)

[Support Suite Home](https://support.parsebiosciences.com/hc/en-us)

![](Adding%20Custom%20Sequences%20and%20Gene%20Annotation%20File%20Formatting%20%E2%80%93%20Support%20Suite%20-%20Parse%20Biosciences_files/default_avatar.png)
Kevin Stachelek

[My activities](https://support.parsebiosciences.com/hc/en-us/requests)
[My profile](https://support.parsebiosciences.com/hc/en-us/profiles/41564339448596)
[Sign out](https://support.parsebiosciences.com/access/logout?return_to=https%3A%2F%2Fsupport.parsebiosciences.com%2Fhc%2Fen-us)

1. [Support Suite - Parse Biosciences](https://support.parsebiosciences.com/hc/en-us)
2. [Data Analysis](https://support.parsebiosciences.com/hc/en-us/categories/360004765711-Data-Analysis)
3. [FASTQ File Processing](https://support.parsebiosciences.com/hc/en-us/sections/27082017166484-FASTQ-File-Processing)
4. [Running the command line pipeline (advanced)](https://support.parsebiosciences.com/hc/en-us/sections/25114901087380-Running-the-command-line-pipeline-advanced)

Articles in this section

* [Adding Custom Sequences and Gene Annotation File Formatting](https://support.parsebiosciences.com/hc/en-us/articles/28706669023636-Adding-Custom-Sequences-and-Gene-Annotation-File-Formatting)
* [Workflow for Processing Long Read Data with the Parse Biosciences Pipeline](https://support.parsebiosciences.com/hc/en-us/articles/17921932763028-Workflow-for-Processing-Long-Read-Data-with-the-Parse-Biosciences-Pipeline)
* [Read and Transcript Filtering with the CRISPR Detect Pipeline](https://support.parsebiosciences.com/hc/en-us/articles/17195970123028-Read-and-Transcript-Filtering-with-the-CRISPR-Detect-Pipeline)
* [Pipeline Options (Current Version)](https://support.parsebiosciences.com/hc/en-us/articles/17173670279828-Pipeline-Options-Current-Version)

# Adding Custom Sequences and Gene Annotation File Formatting

* First Published: July 22, 2024
* Updated: 1 month ago

* [Genome and Gene Annotation Formatting requirements](#h_01K3VNB3B0XC9NF417HDA8GS4Y)
  + [GTF Formatting checklist](#h_01J60CWT8QHJHVPNPCTTSJ6HC2)
* [Adding custom sequences](#h_01J3E2455FS22ZZJTDCM5EGHMX)

Prior to using a genome for processing sequencing data from a Parse
assay, the analysis pipeline requires the preparation of sequence and
gene model data provided in various genome files. Here we'll discuss the
format and features required for the `--genes` (gene model), and the `--fasta` argument in the `split-pipe --mode mkref` command used to create a reference genome. The full `mkref` command is described in the "[Running the Pipeline](https://support.parsebiosciences.com/hc/en-us/articles/23060220494228)" article. Here is an example of a `mkref` command:

```
split-pipe \
--mode mkref \
--genome_name hg38 \
--fasta /newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--genes /newvolume/genomes/Homo_sapiens.GRCh38.108.gtf.gz \
--output_dir /newvolume/genomes/hg38
```
## Genome and Gene Annotation Formatting requirements

### Genome formatting requirements

Genomes are typically formatted as FASTA file. In FASTA format the
line before the nucleotide sequence, called the FASTA definition line,
must begin with a carat (">"), followed by a unique, but arbitrary,
SeqID (sequence identifier). An example is shown below:

```
>sequence_name1
TGACCG..
>sequence_name2
TGACCG..
...
```

**Note:** Some SeqIDs include a description which is
separated by a single space. The description is ignored by the pipeline
and only the SeqID is parsed from the FASTA file.

### Gene annotation formatting requirements

Currently, `split-pipe --mode mkref` accepts and processes GTF and GFF3 formatted files as an argument to the parameter `--genes`. Please note that GTF formatting is equivalent to GFF2.

**Note:** UCSC formatted GTF files are not supported. UCSC gene annotations are typically sourced through [Enesmbl](https://ensembl.org/) or [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/), and we recommend using either source as an alternative.

A GTF file is a tab-delimited file that consists of 9 separate
columns. The contents of each column are described below; Definitions
were taken directly from the [Ensembl website.](http://ensembl.org/)

### **GTF columns:**

1. Seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. **Important**:
   each seqname in the annotation must be identical to the FASTA SeqID
   (described above), without any additional content such as species or
   assembly. See the example GTF output below.
2. Source - name of the program that generated this feature, or the data source (database or project name)
3. Feature - feature type name, e.g. Gene, Variation, Similarity
4. Start - Start position\* of the feature, with sequence numbering starting at 1.
5. End - End position\* of the feature, with sequence numbering starting at 1.
6. Score - A floating point value.
7. Strand - defined as + (forward) or - (reverse).
8. Frame - One of '0', '1' or '2'. '0' indicates that the first base of
   the feature is the first base of a codon, '1' that the second base is
   the first base of a codon, and so on..
9. Attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

The Parse pipeline `mkref` command requires the standard
9-column GTF file. In addition, there are minimal requirements for data
elements in columns 3 and 9. For *each* annotated "feature" (gene
or non-coding region), please go through each step of the formatting
checklist below to ensure your annotation is parsed correctly.

### GTF Formatting checklist

* The file must be **tab delimited**, otherwise it can't
  be parsed. A delimiter is the character used to separate columns on a
  file formatted as a table. A common mistake is to use four spaces
  instead of tabs when adding a custom sequence. The pipeline's  `mkref` command will return "Problem: Problem parsing gtf attributes" in the log file if this is an issue.
* Column 1 - The seqname must have a corresponding and identical FASTA
  identifier. For example, if a FASTA identifier is labeled
  ">CP049705.1", then the seqname in the annotation must be
  "CP049705.1"
* Column 3 - a "**gene**" *and* "**exon**" label. This means that you will have at least **two entries** for every feature (see the example below).
* Column 7 - A positive strand label **"+"**. Reads that map to the reverse strand are not used.
* Column 9 - a "**gene\_id**" and '**gene\_biotype "protein coding**"' *or* '**gene\_type protein coding**'. Note that "biotype" is typically used in Ensembl GTFs, while "type" is more commonly used in Gencode GTFs.

Here is an example of a single feature that has the minimal required inputs. Note that there is one line each with **"gene" and "exon"** as features in column 3. Both lines include a plus sign in the 7th column, and "gene\_biotype" keys in the in column 9.

```
1 havana  gene  11869  14409  .  +  .  gene_id "ENSG00000223972"; gene_biotype "protein_coding";
1 havana  exon  11859  14389  .  +  .  gene_id "ENSG00000223972"; gene_biotype "protein_coding";
```

How can you quickly check to see if your GTF file meets these requirements?

* For column 3 features, run the command `cut -f3 gene_model.gtf | sort | uniq` and you should see both "gene" and "exon" show up in the results.
* For column 9 attributes, simply run `grep -c "gene_type" gene_model.gtf` and `grep -c "gene_biotype" gene_model.gtf`.
  If a match for either entry shows up (i.e. the number reported is
  greater than zero), then you have the required input for this field.

## Adding custom sequences

There may be a scenario where one wishes to detect the expression of a
particular sequence that does not exist in their reference genome or
gene annotation. For example, a researcher may insert a vector
containing a GFP sequence that is only expressed in a specific
population of cells and want to view the resulting expression in their
pipeline output. In the following tutorial, we’ll go over the steps
required to add a custom sequence to a genome fasta and corresponding
GTF file.

### Example: adding a GFP sequence to the human genome

Let's start off by creating a copy of our genome to work with. It's
generally a good idea to create backups when modifying files in case you
encounter any problems along the way.

```
cp Homo_sapiens.GRCh38.dna.primary_assembly.fa hg38_gfp.fa
```

Next, using an example GFP sequence from [Ensembl](http://ensembl.org/),
we'll append the GFP sequence to our genome copy. For simplicity, we'll
change the fasta header containing the chromosome name to "EGFP".

```
# download sequence
wget -O gfp_append.fa "http://ensembl.org/Saccharomyces_cerevisiae/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;g=YMR122W-A;output=fasta;r=XIII:511315-511569;strand=feature;t=YMR122W-A_mRNA;peptide=yes;intron=yes;utr3=yes;exon=yes;genomic=unmasked;cdna=yes;utr5=yes;coding=yes;_format=Text"

# change header
sed -i "1s/.*/>EGFP/g" gfp_append.fa

# append sequence to new copy
cat gfp_append.fa >> hg38_gfp.fa
```

Now we can append the annotation for our new gene to the GTF file. A full description of GTF formatting can be found in the "**[Gene annotation formatting requirements](#01J3E23V5PYZYGQ71G0RNZT226)**" section below, and on **[Ensembl's website.](https://ensembl.org/)** **IMPORTANT:**
Most errors related to custom sequences not showing up are related to
improper GTF formatting. If your feature of interest does not show up in
your results, please carefully review the **[Gene annotation formatting requirements](#01J3E23V5PYZYGQ71G0RNZT226)** section below to ensure the appropriate fields are filled out.

You can use the following `echo` command below as a template for appending your own custom sequence, but there are a few important items to note:

* The `\t` between each column should not be changed. These are tab delimiters required for current GTF standards.
* It's critical that the chromosome name in the first column of the GTF file is an **exact** match with chromosome name in the fasta header you added in the previous step.
* You'll want to change the start and end point of your gene of
  interest in the 4th and 5th columns. "1" (4th column) is the first base
  in our example gene and "4813" (5th column) is the last base. You can
  use the bash command `tail -n +2 gfp_append.fa | wc | awk '{print $3-$1}'` to find the total number of characters (bases) in the sequence. IMPORTANT:
  If the read mapping position does not agree with the start and end
  positions in the GTF file, then you will not see your region of interest
  show up in your expression matrix. This can potentially generate empty
  transcript assignment files and crash the pipeline if you are working
  with targeted library (i.e. a library comprised of all custom
  sequences).
* You can change the values between each `\"` to suit the description for your sequence of interest, but the backslash and double quote should be left untouched.
* The characters ".", "+", and "." (columns 6, 7, and 8) do need any modification.
* Lastly, you'll append the transgene to the full GTF with two nearly
  identical lines, except one contains "gene" in the feature column and
  the other contains "exon"; both of these are required by STAR. It's a
  good idea to test your command to ensure the output is correct before
  appending the new text to your GTF file.  We've included this in
  the second command below.

```
# creating a copy
cp Homo_sapiens.GRCh38.93.gtf hg38_ens93_gfp.gtf

# testing our echo command
echo -e "EGFP\tcustom\tgene\t1\t4813\t.\t+\t.\tgene_id \"id_EGFP\"; gene_name \"EGFP\"; gene_biotype \"protein_coding\";"
echo -e "EGFP\tcustom\texon\t1\t4813\t.\t+\t.\tgene_id \"id_EGFP\"; gene_name \"EGFP\"; gene_biotype \"protein_coding\";"

# using the echo command to append both "gene" and "exon" to our GTF; both are required.
echo -e "EGFP\tcustom\tgene\t1\t4813\t.\t+\t.\tgene_id \"id_EGFP\"; gene_name \"EGFP\"; gene_biotype \"protein_coding\";"  >> hg38_ens93_gfp.gtf
echo -e "EGFP\tcustom\texon\t1\t4813\t.\t+\t.\tgene_id \"id_EGFP\"; gene_name \"EGFP\"; gene_biotype \"protein_coding\";"  >> hg38_ens93_gfp.gtf
```

That's it! In order to apply the new changes you'll need to generate a new alignment index with  `split-pipe mode mkref` and follow up with `split-pipe mode all`.

Was this article helpful?
Yes
No

0 out of 1 found this helpful

Have more questions? Please contact support.

[Return to top](#article-container)

## Recently viewed articles

* [Pipeline Options (Current Version)](https://support.parsebiosciences.com/hc/en-us/articles/17173670279828-Pipeline-Options-Current-Version)
* [Pipeline Change Log](https://support.parsebiosciences.com/hc/en-us/articles/17194120578708-Pipeline-Change-Log)
* [Running the Pipeline (Current Version)](https://support.parsebiosciences.com/hc/en-us/articles/23060220494228-Running-the-Pipeline-Current-Version)
* [Pipeline Download (Current Version)](https://support.parsebiosciences.com/hc/en-us/articles/17200056667924-Pipeline-Download-Current-Version)
* [Pipeline Installation (Current Version)](https://support.parsebiosciences.com/hc/en-us/articles/23060102930580-Pipeline-Installation-Current-Version)

## Related articles

* [Running the Pipeline (Current Version)](https://support.parsebiosciences.com/hc/en-us/related/click?data=BAh7CjobZGVzdGluYXRpb25fYXJ0aWNsZV9pZGwrCJRBQCD5FDoYcmVmZXJyZXJfYXJ0aWNsZV9pZGwrCJRRz8obGjoLbG9jYWxlSSIKZW4tdXMGOgZFVDoIdXJsSSJLL2hjL2VuLXVzL2FydGljbGVzLzIzMDYwMjIwNDk0MjI4LVJ1bm5pbmctdGhlLVBpcGVsaW5lLUN1cnJlbnQtVmVyc2lvbgY7CFQ6CXJhbmtpBg%3D%3D--6b2329227a014de88d02bfc57a71466b08276996)
* [Trailmaker User Guide](https://support.parsebiosciences.com/hc/en-us/related/click?data=BAh7CjobZGVzdGluYXRpb25fYXJ0aWNsZV9pZGwrCJQ%2BBUigGDoYcmVmZXJyZXJfYXJ0aWNsZV9pZGwrCJRRz8obGjoLbG9jYWxlSSIKZW4tdXMGOgZFVDoIdXJsSSI8L2hjL2VuLXVzL2FydGljbGVzLzI3MDc2NjgyMTM3MjM2LVRyYWlsbWFrZXItVXNlci1HdWlkZQY7CFQ6CXJhbmtpBw%3D%3D--565dc35810fff05bcab6d113857b85446d9fc38a)
* [Seurat Tutorial - 65k PBMCs](https://support.parsebiosciences.com/hc/en-us/related/click?data=BAh7CjobZGVzdGluYXRpb25fYXJ0aWNsZV9pZGwrCEz41dRTADoYcmVmZXJyZXJfYXJ0aWNsZV9pZGwrCJRRz8obGjoLbG9jYWxlSSIKZW4tdXMGOgZFVDoIdXJsSSI%2BL2hjL2VuLXVzL2FydGljbGVzLzM2MDA1MzA3ODA5Mi1TZXVyYXQtVHV0b3JpYWwtNjVrLVBCTUNzBjsIVDoJcmFua2kI--907ec1ef421dfad661e460b9785bf72995877753)
* [Content of fastq and bam files (barcode and other annotations)](https://support.parsebiosciences.com/hc/en-us/related/click?data=BAh7CjobZGVzdGluYXRpb25fYXJ0aWNsZV9pZGwrCBRoa6YDBDoYcmVmZXJyZXJfYXJ0aWNsZV9pZGwrCJRRz8obGjoLbG9jYWxlSSIKZW4tdXMGOgZFVDoIdXJsSSJiL2hjL2VuLXVzL2FydGljbGVzLzQ0MTM3MjM0Njk4NDQtQ29udGVudC1vZi1mYXN0cS1hbmQtYmFtLWZpbGVzLWJhcmNvZGUtYW5kLW90aGVyLWFubm90YXRpb25zBjsIVDoJcmFua2kJ--24d252ef8371845cec86d4be0c807d3e454bb814)
* [How Do I Analyze my Parse Biosciences Data?](https://support.parsebiosciences.com/hc/en-us/related/click?data=BAh7CjobZGVzdGluYXRpb25fYXJ0aWNsZV9pZGwrCJRx6uKdGDoYcmVmZXJyZXJfYXJ0aWNsZV9pZGwrCJRRz8obGjoLbG9jYWxlSSIKZW4tdXMGOgZFVDoIdXJsSSJRL2hjL2VuLXVzL2FydGljbGVzLzI3MDY2Mzk1OTQ3NDEyLUhvdy1Eby1JLUFuYWx5emUtbXktUGFyc2UtQmlvc2NpZW5jZXMtRGF0YQY7CFQ6CXJhbmtpCg%3D%3D--76f957da80d3c888fa63336b73ad36e9727f3e11)

[![Parse Biosciences](Adding%20Custom%20Sequences%20and%20Gene%20Annotation%20File%20Formatting%20%E2%80%93%20Support%20Suite%20-%20Parse%20Biosciences_files/logo-white.svg)](https://parsebiosciences.com/ "Home page")

700 Dexter Ave

Suite 600

Seattle, WA 98109

info@parsebiosciences.com

#### Products

* [Evercode Whole Transcriptome](https://parsebiosciences.com/products/evercode-whole-transcriptome/)
* [Evercode TCR](https://parsebiosciences.com/products/evercode-tcr/)
* [Gene Capture](https://parsebiosciences.com/products/gene-select/)
* [CRISPR Detect](https://parsebiosciences.com/products/crispr-detect/)
* [Data Analysis](https://parsebiosciences.com/products/data-analysis/)

#### Company

* [Technology](https://parsebiosciences.com/technology/)
* [Company](https://parsebiosciences.com/company/)
* [News](https://parsebiosciences.com/news/)
* [Upcoming Events](https://parsebiosciences.com/events/)
* [Careers](https://parsebiosciences.com/careers/)
* [Contact](https://parsebiosciences.com/contact/)

#### Resources

* [Technical Brochure](https://www.parsebiosciences.com/datasets/tech-brochure/evercode-whole-transcriptome-technical-brochure/)
* [Datasets & Product Literature](https://www.parsebiosciences.com/resources/datasets/)
* [Publications](https://www.parsebiosciences.com/resources/publications/)
* [Recorded Webinars](https://www.parsebiosciences.com/resources/webinars/)
* [Support Suite](https://support.parsebiosciences.com/hc/en-us)
* [Privacy Policy](https://parsebiosciences.com/privacypolicy)
* [Terms of Service](https://parsebiosciences.com/terms-of-service)

Stay up to date on the latest news and announcements.

Email\*Generate new mask

© 2025 Parse Biosciences. All rights reserved.


