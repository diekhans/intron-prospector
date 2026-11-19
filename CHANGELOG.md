# Major and user-visible changes

## version 1.6.0 2026-11-19
* A read supporting a junction is counted once, as `SJ.out.tab` counts it.  A
  pair whose mates both cross one junction was counted twice, and a
  multi-mapping read was counted once for each of its alignment records at the
  junction.  The overhangs still take every alignment, so an uncounted mate
  with a longer anchor still sets `max_left_overhang` or
  `max_right_overhang`.  The confidence is computed over the offsets of the
  reads counted, so a junction supported by one pair whose mates both cross it
  now has one offset rather than two and a confidence of 0, which the default
  `--min-confidence-score=1.0` drops.  On the WTC11 benchmark the counts fall
  by a median of 24% and then equal the `SJ.out.tab` counts on 99.1% of the
  junctions the two share. [incompatibility]
* BED color for GC/AG introns is now U2 green rather than unknown red.
* Help text and man page regenerated; --strandness description was attached to
  --unsorted and --set-TS-strand-tag had a stray copy of the --skip-missing-targets
  text.
* Fixed --unsorted --pass-through writing a header and no alignments.
* Fixed junctions being reported on the wrong chromosome when several input files
  order their headers differently.
* Fixed crash on an input file with no readable SAM header; it is now an error
  naming the file.
* Fixed strand of single-end reads under --strandness, which assumed the paired-end
  flags.
* Anchor length no longer counts indels seen before the anchor started, which
  dropped junctions whose reported overhang passed --min-anchor-length.
* Fixed -d and -m short options being rejected.
* Fixed the TSV column count error reporting the expected count as the actual one.
* Initialized the unsorted flag, which was read uninitialized.
* The splice sites now decide a junction's orientation: sites that are
  canonical only reverse complemented set the strand, rather than being
  reported as non-canonical and dropped by --sj-filter=canon.  A wrong
  --strandness or XS tag no longer loses junctions.
* --sj-filter=canon records its decision per junction, so the reads of a
  rejected junction are no longer counted into a later one and the counts no
  longer depend on read order.
* A failed write is now an error rather than a silent truncation, in both
  programs.
* Intron coordinates are 64 bit, so positions past 2^32 are no longer wrapped.
* Read offsets in the confidence score are 32 bit, so offsets past 65535 no
  longer collide and zero the score of a long read junction.
* A padding CIGAR operation is ignored rather than ending the run.
* --debug-trace names its junc_strand column, which the header left out while
  the rows wrote it.
* The strand of a read whose mate is unmapped or on another target comes from
  its own orientation, as for an unpaired read, rather than being undetermined.
* The error for several input BAMs names --pass-through, the option that
  exists.
* --min-anchor-length is checked rather than parsed with atoi, where a bad
  value silently became 0 and turned the anchor test off.
* A TSV read by intron-prospector-merge may have CRLF line endings; the
  carriage return is no longer left on the last column.
* The build runs with -Wall -Wextra and is warning free.
* intron-prospector-merge no longer discards the reads of a record whose
  strand disagrees with the others, issue 15.  It warns, sums the counts, and
  takes the strand with the most reads, counted over every merged record so
  the result does not depend on the order the files are given in.
* A failure while running, such as an unsorted BAM or a missing genome
  sequence, is no longer reported as a command line error with a pointer to
  --help.
* --exclude names categories of alignment by SAM flag: qcfail, dup, secondary
  and supplementary.  --include uses a category the default skips, so
  --include=qcfail counts the reads a sequencer chastity filter flagged, which
  is what STAR counts in SJ.out.tab.  Either option takes a comma separated
  list and may be repeated.  The default is what it was: qcfail, dup and
  secondary skipped, supplementary used, unmapped always skipped.
* The category of multi-mapped reads is named multimapping.  It is the one
  category that is not a SAM flag: it reads the NH tag, and it drops every
  alignment of a read that aligns in more than one place, where secondary
  drops only the extra alignment records. [incompatibility]
* --excludes is deprecated.  It warns and is otherwise --exclude.
* A build from a git tree reports the last release version with +main
  appended, such as 1.5.1+main, so it is not mistaken for a release.
* The release process is documented in DEVELOPMENT.md.

## version 1.5.1 2026-04-22
* Fixed bug where -S option consumed the next argument

## version 1.5.0 2026-04-04
* Add --min-read-count option to filter introns based on number of supporting reads.
* If htslib is not found by configure, error message makes suggestion when using Conda.
* The intronProspector and intronProspectorMerge commands have been renamed to 
  intron-prospector and intron-prospector-merge to follow more common conventions and
  decrease wear on the shift key.  Wrapper scripts are provide for compatibility, which
  will be removed in the next major release.

## version 1.4.0 2025-08-12
* Add --sj-filter=canon option `intronProspector` to keep only canonical splice junctions.
  This is the default if a genome is provide.  Filtering can be disabled with
  --sj-filter=all. [incompatibility]
* Add --sj-filter=canon option `intronProspectorMerge` to keep only canonical splice
  junctions.
* Make --min-confidence-score=1.0 the default. May want to lower or disable for low-coverage
  or long-reads. [incompatibility]


## version 1.3.0 2025-04-10
* Added ability to read and write gzip compressed TSV and BED files based on an name
  extension of '.gz'


## version 1.2.0 2025-03-15
* Added ability to process multiple BAMs at once.  This will use more memory
  and possibly improve the confidence scoring in some cases
* Added option to allow input BAM to be unsorted.  This will use more memory.
