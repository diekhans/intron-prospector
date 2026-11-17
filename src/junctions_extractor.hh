/*  junctions_extractor.hh -- splice junction extraction.

    Copyright (c) 2018-2025, Mark Diekhans, University of California, Santa Cruz

This code is derived from the regtools package available at:

  https://github.com/griffithlab/regtools

We thank them for providing this software under a flexible license:

    Copyright (c) 2015, The Griffith Lab

    Author: Avinash Ramu <aramu@genome.wustl.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef JUNCTIONS_EXTRACTOR_H
#define JUNCTIONS_EXTRACTOR_H

#include <assert.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <map>
#include <set>
#include <vector>
#include "htslib/sam.h"
#include "junctions.hh"
#include "splice_juncs.hh"

class Genome;

using namespace std;

// Strandness of data; 0 = unstranded, 1 = RF, 2 = FR
// DO NOT CHANGE VALUES, code depends on it.
typedef enum {
    UNSTRANDED = 0,
    RF_STRANDED = 1,
    FR_STRANDED = 2
} Strandness;


/* Categories of alignment that --include and --exclude name, as a mask.
 *
 * Four of them are SAM flags, so those bits of the mask are the flag bits and
 * are tested against the alignment's flags.  The fifth, multimapping, is read from
 * the NH tag rather than from a flag and uses a bit of its own, above any SAM
 * flag.  It also removes a different thing: secondary drops one alignment
 * record of a read that aligns in several places, leaving its primary
 * alignment, while multimapping drops every alignment of such a read.
 *
 * An unmapped alignment is not a category.  It has no junctions and is always
 * skipped.
 */
static const unsigned ALIGN_CATS_NONE = 0;

// not a SAM flag; the NH tag says how many places the read aligns
static const unsigned ALIGN_CAT_MULTIMAPPING = 0x10000;

// the categories that are SAM flags, so can be tested against aln->core.flag
static const unsigned ALIGN_CATS_FLAGS = BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY;

// excluded unless --include names them; supplementary alignments and
// multi-mapped reads are used
static const unsigned ALIGN_CATS_EXCLUDED_DEFAULT = BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY;

// parse a comma separated category list into a mask of SAM flags
unsigned align_categories_parse(const string& spec);

// category names of a mask, comma separated, for messages
string align_categories_to_string(unsigned cats);

// the category names, for help and error messages
string align_categories_help();


/* Which templates have already counted a junction.
 *
 * A read pair whose mates both cross one junction is one read of support, as
 * SJ.out.tab counts it, and so is a multi-mapping read whose records repeat a
 * junction.  The first alignment to count a junction makes an entry here,
 * keyed by the read name and the junction, and any further alignment of that
 * template finds it and does not count the junction again.
 *
 * An alignment supporting a junction has matched bases before the intron, so
 * it starts before intron_start.  Once a coordinate sorted scan has passed
 * intron_start no further alignment can support that junction, which is when
 * retire_before drops the entry.  Only the junctions within a read length of
 * the scan are held.  The unsorted path has no such point, so its entries
 * accumulate until the file has been read.
 *
 * The read name is stored as a 64 bit hash.  Two names of one junction that
 * hash alike would be counted as one read rather than two.
 */
class CountedTemplates {
    private:
    typedef pair<uint64_t, const Junction*> Entry;

    // templates that have counted a junction
    set<Entry> entries_;

    // entries by the junction start that retires them
    multimap<hts_pos_t, Entry> by_intron_start_;

    static uint64_t name_hash(bam1_t *aln) {
        return hash<string>()(string(bam_get_qname(aln)));
    }

    /* Can this template have another alignment at this junction?  Its mate
     * can, when both are on this target, and so can another record of a read
     * that aligns in more than one place.  A single-end read whose NH tag says
     * it aligns once cannot, so nothing is recorded for it; an NH that is
     * missing leaves it open, so an entry is made.
     */
    static bool template_can_repeat(bam1_t *aln, int num_aligns) {
        return (((aln->core.flag & BAM_FPAIRED)
                 and not (aln->core.flag & BAM_FMUNMAP)
                 and (aln->core.mtid == aln->core.tid))
                or (num_aligns != 1));
    }

    public:
    // has another alignment of this template counted this junction?
    bool already_counted(bam1_t *aln, const Junction *junc) const {
        // the common case is nothing recorded, which skips hashing the name
        return (not entries_.empty())
            and (entries_.count(Entry(name_hash(aln), junc)) > 0);
    }

    // record that this alignment counted this junction
    void record(bam1_t *aln, int num_aligns, const Junction *junc) {
        if (template_can_repeat(aln, num_aligns)) {
            Entry entry(name_hash(aln), junc);
            entries_.insert(entry);
            by_intron_start_.insert(pair<hts_pos_t, Entry>(junc->intron_start, entry));
        }
    }

    // drop the entries for junctions the scan has passed
    void retire_before(hts_pos_t pos) {
        multimap<hts_pos_t, Entry>::iterator it = by_intron_start_.begin();
        while ((it != by_intron_start_.end()) and (it->first < pos)) {
            entries_.erase(it->second);
            by_intron_start_.erase(it++);
        }
    }

    void clear() {
        entries_.clear();
        by_intron_start_.clear();
    }
};

// The class that deals with creating the junctions.
// Coordinate-sorted BAM are require to process one target
// at a time, reducing memory usage.
class JunctionsExtractor {
private:
    // Minimum anchor length for junctions
    // Junctions need atleast this many bp overlap
    // on both ends.  Mismatch bases are not included in the count.
    uint32_t min_anchor_length_;

    // maximum size of a non-intron indel to allow in an anchor.  Normally zero for short reads.
    uint32_t max_anchor_indel_size_;
    
    // Minimum length of an intron, i.e min junction width
    uint32_t min_intron_length_;

    // Maximum length of an intron, i.e max junction width
    uint32_t max_intron_length_;

    // filter for types of junctions to keep
    JunctionFilter junction_filter_;

    
    //strandness of data; 0 = unstranded, 1 = RF, 2 = FR
    Strandness strandness_;

    // SAM flags of the alignment categories to skip
    unsigned exclude_flags_;

    // skip reads the NH tag says align in more than one place
    bool exclude_multimapping_;
    
    // Alignment file
    string bam_;

    // Object to get genome sequences, which maybe NULL
    Genome *genome_;

    // missing genomic targets
    bool skip_missing_targets_;
    
    // target index to target (chrom) name
    vector<string> targets_;
    
    // Map to store the junctions by intron coordinates
    JunctionTable junctions_;

    // junctions the splice site filter rejected, so the decision is made once
    set<JunctionKey> rejected_;

    // junctions already counted for a template, so each is counted once
    CountedTemplates counted_templates_;

    // missing genomic sequence that have been warned about
    set<string> missing_genomic_warned_;

    // update tags from observed orientation
    bool set_XS_strand_tag_;
    bool set_TS_strand_tag_;

    // checking for not sorted
    int previous_target_;

    // alignment buffer, used to handle one level of read-ahead
    bam1_t *aln_buf_;
    bool aln_pending_;

    // used for generate BED names
    uint32_t ijunc_;

    // BAM files
    samFile *in_sam_;
    bam_hdr_t *in_header_;
    samFile *out_sam_;
    
    // debugging trace output if not NULL
    ostream *trace_fh_;
    
    // internal functions
    void save_targets(bam_hdr_t *header);
    bool junction_qc(hts_pos_t anchor_start, hts_pos_t intron_start,
                     hts_pos_t intron_end, hts_pos_t anchor_end,
                     uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt,
                     uint32_t left_indel_cnt, uint32_t right_indel_cnt);
    void process_alignment(bam1_t *aln);
    void write_pass_through(bam1_t *aln,
                            bam_hdr_t *in_header,
                            samFile* out_sam);
    
    Junction *create_junction(bam1_t *aln, const JunctionKey &key,
                              const string& chrom, char strand,
                              hts_pos_t anchor_start, hts_pos_t intron_start,
                              hts_pos_t intron_end, hts_pos_t anchor_end,
                              const string& splice_sites);
    Junction *update_junction(bam1_t *aln, const JunctionKey &key,
                              const string& chrom, char strand,
                              hts_pos_t anchor_start, hts_pos_t intron_start,
                              hts_pos_t intron_end, hts_pos_t anchor_end);
    Junction *add_junction(bam1_t *aln, const string& chrom, char strand,
                           hts_pos_t anchor_start, hts_pos_t intron_start, 
                           hts_pos_t intron_end, hts_pos_t anchor_end);
    char get_junction_strand_XS(bam1_t *aln);
    char get_junction_strand_flag(bam1_t *aln);
    string get_splice_sites(const string& chrom, char strand,
                            hts_pos_t intron_start, hts_pos_t intron_end);
    char correct_strand(char strand,
                        string& splice_sites);
    bam1_t* read_align();
    bool process_target_alignment(int target_index,
                                  bam1_t *aln);

    JunctionsExtractor() {
        assert(false); // Default constructor not allowed
    }
public:
    JunctionsExtractor(uint32_t min_anchor_length,
                       uint32_t max_anchor_indel_size,
                       uint32_t min_intron_length,
                       uint32_t max_intron_length,
                       JunctionFilter junction_filter,
                       Strandness strandness,
                       unsigned exclude_cats,
                       Genome *genome,
                       bool skip_missing_targets,
                       bool set_XS_strand_tag,
                       bool set_TS_strand_tag,
                       ostream *trace_fh):
        min_anchor_length_(min_anchor_length),
        max_anchor_indel_size_(max_anchor_indel_size),
        min_intron_length_(min_intron_length),
        max_intron_length_(max_intron_length),
        junction_filter_(junction_filter),
        strandness_(strandness),
        exclude_flags_(exclude_cats & ALIGN_CATS_FLAGS),
        exclude_multimapping_((exclude_cats & ALIGN_CAT_MULTIMAPPING) != 0),
        genome_(genome),
        skip_missing_targets_(skip_missing_targets),
        set_XS_strand_tag_(set_XS_strand_tag),
        set_TS_strand_tag_(set_TS_strand_tag),
        previous_target_(-1),
        aln_buf_(NULL),
        aln_pending_(false),
        ijunc_(0),
        in_sam_(NULL),
        in_header_(NULL),
        out_sam_(NULL),
        trace_fh_(trace_fh) {
    }

    ~JunctionsExtractor();

    // Get the strand
    char get_junction_strand(bam1_t *aln);

    // Validate a junction and save if it passes.
    void record_junction(bam1_t *aln, const string& chrom, char strand,
                         hts_pos_t anchor_start, hts_pos_t intron_start,
                         hts_pos_t intron_end, hts_pos_t anchor_end,
                         uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt,
                         uint32_t left_indel_cnt, uint32_t right_indel_cnt,
                         int *orientCnt);

    // free junctions found so far
    void clear() {
        junctions_.clear();
        rejected_.clear();
        counted_templates_.clear();
    }
    
    // Open a BAM file for process
    void open(const string& bam);

    // open the pass through, must be done after open.
    void open_pass_through(const string& bam_pass_through);
    
    // close current BAM file, also closes pass-through if open
    void close();

    // get the number of targets
    int get_num_targets() const {
        return targets_.size();
    }
    
    // Identify exon-exon junctions. This requires BAM to be sorted and
    // works with pass-through.  It also optimizes memory usage.
    void identify_junctions_for_target(int target_index);

    // Copy remaining reads (unaligned) to pass-through, if it is open
    void copy_unaligned_reads();

    // Identifies junctions for an entire BAM. This does not require
    // the BAM to be sorted, however it uses more memory
    void identify_junctions_for_bam();

    // Get a vector of all the junctions
    JunctionVector get_junctions() {
        return junctions_.get_junctions();
    }
};
#endif
