# NAME

**intron-prospector** — Extract putative intron junctions from RNA-Seq alignments

# SYNOPSIS

`intron-prospector [options] [readaligns ...]`

# DESCRIPTION

Find putative intron junctions in a RNA-Seq alignment. The *readaligns* files maybe in SAM, BAM, or CRAM format and does not need to be sorted or indexed and maybe streamed. If omitted, stdin is used.

This program allows for integrating splice junction calling into an alignment pipeline.  Pass-through mode copys the alignment file on `stdin` to `stdout`.  It can then sit between an aligner outputting SAM and `samtools` that converts to BAM/CRAM.

Both short-read (Illumina) and long-read RNA-Seq can be process.  For
short-reads, it is recommended to use `--min-confidence-score=1.0`, which is the default. 
For long-read experiments this should be disable with`--min-confidence-score=0.0`.
Low coverage experiment may wish to disable or use a lower confidence filter. 
Introns are as determined by the aligner and indicated in the BAM by `N` operations.

If a genome is provided, only canonical introns are kept (AT/AG, GC/AG, and AT/AC).  This can be overridden with `--sj-filter=all`

Process multiple BAMs at once will use more more memory, as intron counts are
not collected per chromosome.  This mode does not require a sorted BAM and is
incompatible with `--pass-through`.  It may also improve the confidence
scoring in some cases.

The `intron-prospector-merge` program can be used to convert to from the
`--intron-calls` format to other formats as well as merge the output from
multiple `intron-prospector` runs.

TSV and BED files will be automatically compressed with `gzip` if they end in `.gz`.

# OPTIONS

`-h, --help`

> Prints brief usage information.

`-v, --version`

> Prints the current version number.

`-g fasta, --genome-fasta=fasta`

> Genome FASTA file, must be indexed by `samtools faidx`.  Donor and acceptor dinuclotides are identified if provided. THIS SHOULD NORMALLY BE PROVIDE!

`-S, --skip-missing-targets`

> If a target sequence is not in the genome, skip the alignment rather than generate an error.

`-f spec, --sj-filter=spec`

> Filter based on splice junctions motif. A spec of `canon` keeps only canonical splice junctions (GT/AG, GC/AG, and AT/AC).  A spec of `all` passes all through.  If `canon` is specified, a genome must be supplied, and if a genome is supplied, `canon` becomes the default.

`-a INT, --min-anchor-length=INT`

> Minimum anchor length. Junctions which satisfy a minimum anchor length on both sides are reported.  Mismatch bases don't count towards the meeting this threshold.  The default is 8 bases.

`-i INT, --min-intron-length=INT`

> Minimum intron length. The default is 70 bases.

`-I INT,  --max-intron-length=INT`

> Maximum intron length. The default is 500000 bases.

`-d,  --allow-anchor-indels`

> Are indels allowed in anchors?  Useful for ONT reads.  Indels bases don't count towards the meeting the anchor size threshold.  The default is don't allow indels. .

`-m INT,  --max-anchor-indel-size=INT`

> Maximum size of any contiguous indel in an anchor if indels are allowed in anchors.  Exceeding this length discards the intron.  The default is 36.

`-C FLOAT, --min-confidence-score=FLOAT`

> Calculate the Shannon-Wiener Diversity Index to use as a confidence score and discard intron calls below this value. A value of 1.0 is a good threshold for filtering. A value of 0.0 disables filtering.  The default is 1.0. This methodology it taken from JuncBASE DOI: 10.1101/gr.108662.110.

`-r INT, --min-read-count=INT`

> Minimum number of reads supporting a junction for it to be reported. Junctions with fewer reads than this threshold are discarded. The default is 0 (no filtering).

`-s STRING, --strandness=STRING`

> Strand specificity of RNA library preparation.  Use `UN` for unstranded, `RF` for first-strand, `FR` for second-strand (case-insensitive).  The default is `UN`.  This is used to set the strand in the junction-format BED file.

`-u, --unsorted`

> SAM/BAM files are not sorted.  This will require more memory, as intron counts are
> not collected per chromosome.

`-X categories, --exclude=categories`

> Skip alignments in these categories, given as a comma separated list or by
> repeating the option.  Four of the categories name a SAM flag:
> *qcfail* - did not pass the sequencer's filters, flag 0x200;
> *dup* - PCR or optical duplicate, flag 0x400;
> *secondary* - another alignment of a read that aligns in more than one place,
> flag 0x100;
> *supplementary* - part of a chimeric alignment, flag 0x800.

> The fifth, *multimapping*, is a read that aligns in more than one place.  It
> is not a SAM flag: it is the NH tag, which counts the places the read aligns.
> It also removes something different from *secondary*: skipping *secondary*
> drops the extra alignment records and still uses the read's primary
> alignment, while skipping *multimapping* drops every alignment of that read.
> A read with no NH tag is not *multimapping*, since nothing says how many
> places it aligns; its reads are counted in `unsure_mapped_count`.

> By default *qcfail*, *dup* and *secondary* are skipped, and *supplementary*
> and *multimapping* are used.  An unmapped alignment contains no introns, so
> it is always skipped and is not a category.

`--include=categories`

> Use alignments in these categories, undoing the default for them.  As with
> `--exclude`, give a comma separated list or repeat the option.
> `--include=qcfail` counts the reads the sequencer's chastity filter flagged,
> which is what STAR counts in `SJ.out.tab`.  Naming a category in both
> `--include` and `--exclude` is an error.

`--set-XS-strand-tag`

> Set the XS:A tag to the strand if the strand for a read can be determined from the
> introns based on a count of recognized vs unrecognized splice sites.
> If the strand can't be determined, the tags are not modified, possible leaving
> an existing XS:A tag in place. Requires `--pass-through` and `--genome-fasta`.

`--set-TS-strand-tag`

> Set the TS:A tag in the same manner as `--set-XS-strand-tag`.

`-c FILE, --intron-calls=FILE`

> Write a tab-separated values (TSV) file, with the junctions calls.  It will contain the following columns :
* chrom - chromosome name
* intron_start - zero based starting coordinates of intron.
* intron_end - zero based ending coordinates of intron.
* strand - `+`, `-`, or `.` if not known.  This is set only if experiment is stranded, it is not set from genomic sequence.
* uniq_mapped_count - number of uniquely mapped reads supporting the junction.
* multi_mapped_count - number of multi-mapped reads supporting the junction.
* unsure_mapped_count - number of reads that are either in discordant or partial mapped pair or do not have the BAM `NH` tag.
* max_left_overhang - maximum number of bases overlapping the exon upstream of the intron.
* max_right_overhang - maximum number of bases overlapping the exon downstream of the intron.
* confidence - confidence score.
* splice_sites - Donor-acceptor pair, if genome provided. In upper-case for known splicing patterns, lower-case for unknown.  If strand is not known, this field is empty.

`-j FILE, --junction-bed=FILE`

> Write the junction calls and support anchors to this BED 12 file.  This is in the same format as TopHat `junctions.bed` and `regtools junction extract` output.  It is a UCSC BED track, with each junction consists of two connected BED blocks, where each block is as long as the maximal overhang of any read spanning the junction. The score is the number of alignments spanning the junction, with a maximum score of 1000 for UCSC browser compatibility.  If genome is supplied, BED is colored green for U2 junctions, blue for U12, or red for unknown.

`-n FILE, --intron-bed=FILE`

> Write the intron BED 9 with the bounds of the introns. The score is the number of alignments spanning the junction, with a maximum score of 1000 for UCSC browser compatibility.  If genome is supplied, BED is colored green for U2 junctions, blue for U12, or red for unknown.

`-b FILE, --intron-bed6=FILE`

> Write the intron BED 6 with the bounds of the introns. The score is the number of alignments spanning the junction.   This is for software not wanting to create a browser track.

`-p FILE, --pass-through=FILE`

> Pass through input BAM/SAM records to this file, used for constructing pipelines with `/dev/stdout` is specified.  CRAM output is not support, these will be written as BAM.

`-D FILE, --debug-trace=FILE`

> Output records, in TSV format, for reach read intron indicating the information going into classifying it, including read name.  First few columns are BED-like for easy conversion.

# NOTES
The computation of strand is problematic.  If the strandness of the experiment is specified, then that is used to determine stand.  If the alignment provides the XS attribute, that is used.  Otherwise, the strand can't be determined from the BAM.  If the genome is provided and a known splice sites are detected, this is then used if the stand is not identified by other methods.

Secondary alignments are not used to support introns unless
`--include=secondary` is given.

# COMPARISON WITH STAR SJ.out.tab

STAR writes its own junction table while aligning.  Its counts and
intron-prospector's are not expected to be identical, for reasons that are not
defects in either program.

*Which alignments are counted*

- STAR counts the alignments the sequencer's chastity filter flagged, SAM flag
  0x200; these are skipped unless `--include=qcfail` is given.  This is
  normally the largest difference.  On an Illumina WTC11 BAM where 5.3% of
  reads are flagged, with the same anchor and read count thresholds applied
  to both sides, `SJ.out.tab` held 6,376 junctions the calls did not;
  `--include=qcfail` left 270 of them.
- STAR counts secondary alignments in column 8; these are skipped unless
  `--include=secondary` is given.  `--exclude=multimapping` goes further and
  drops every alignment of a read that aligns in more than one place, which
  STAR has no option to do when writing the tab.
- `uniq_mapped_count` and `multi_mapped_count` do not partition the reads the
  way STAR's columns 7 and 8 do.  An alignment with no NH tag, a discordant
  pair, or a mate that did not map is counted in `unsure_mapped_count`, which
  `SJ.out.tab` has no counterpart for.
- A read supporting a junction is counted once, as STAR counts it.  A pair
  whose mates both cross one junction is one read, and so is a multi-mapping
  read whose alignment records repeat a junction.  The overhangs take every
  alignment, counted or not, so a mate with a longer anchor still sets
  `max_left_overhang` or `max_right_overhang`.

*Which junctions are written*

- STAR applies `--outSJfilterOverhangMin` (30 12 12 12, by motif),
  `--outSJfilterCountUniqueMin` and `--outSJfilterCountTotalMin` (3 1 1 1),
  `--outSJfilterDistToOtherSJmin` and `--outSJfilterIntronMaxVsReadN` to
  `SJ.out.tab` only, not to the alignments it writes.  A junction can be in the
  BAM and not in the tab.
- Junctions in an annotation given to STAR bypass those filters, so
  `SJ.out.tab` holds rows with an overhang of 1.
- `--min-anchor-length=12` is the closest match to STAR's overhang threshold
  for canonical junctions.

*Columns that do not line up*

- Column 9 of `SJ.out.tab` is the maximum over reads of the smaller of the two
  overhangs.  `max_left_overhang` and `max_right_overhang` are independent
  maxima that can come from different reads, so the smaller of the two is an
  upper bound on STAR's value, not the same number.
- Columns 2 and 3 of `SJ.out.tab` are the first and last intron base, 1-based
  inclusive.  The intron calls TSV and the BEDs are 0-based half open, so the
  start differs by one and the end is the same number.
- STAR classifies the intron motif itself; `splice_sites` is reported only with
  `--genome-fasta` and is empty otherwise.
- STAR marks junctions that are in the annotation it was given, column 6.
  intron-prospector reads no annotation.
- `confidence` has no counterpart in `SJ.out.tab`.

# EXAMPLES

Call junctions from a BAM file, also creating BEDs of junctions and introns:
```
intron-prospector --genome-fasta=thegenome.fa.gz --intron-calls=introns.tsv --junction-bed=juncs.bed --intron-bed=introns.bed reads.bam
```

Pipeline to call introns and create a CRAM file:
```
cat reads.sam \
    | samtools sort -O sam  \
    | ./intron-prospector -c introns.tsv --genome-fasta=thegenome.fa.gz -p /dev/stdout \
    | samtools view -O CRAM -T grch38.fa >reads.cram
```
Note that the `cat` command could be an aligner outputting a SAM file and that the genome FASTA file must be index by `samtools faidx`.

 
# BUGS

See GitHub Issues: <https://github.com/diekhans/intron-prospector/issues>

# AUTHOR

Mark Diekhans <markd@ucsc.edu>

Source available from <https://github.com/diekhans/intron-prospector>

Base on code from RegTools <https://github.com/griffithlab/regtools>
by Avinash Ramu <aramu@genome.wustl.edu>.

