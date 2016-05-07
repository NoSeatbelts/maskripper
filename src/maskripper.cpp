#include <getopt.h>
#include "dlib/bam_util.h"

int usage(char **argv, int retcode=EXIT_FAILURE) {
    fprintf(stderr, "maskripper <opts> in.bam out.bam\n"
                    "Flags:\n-m: Minimum trimmed read length. Default: 0.\n"
                    "-l: output compression level. Default: 6.\n"
                    "Use - for stdin or stdout.\n");
    return retcode;
}
struct opts_t {
    uint8_t *data;
    uint32_t l_data:16;
    uint32_t m_data:16;
    uint32_t n_cigar:8;
    uint32_t min_trimmed_len:8;
    ~opts_t() {if(data) free(data);}
    void resize(uint32_t new_min) {
        if(new_min > m_data) {
            m_data = new_min;
            kroundup32(m_data);
            data = (uint8_t *)realloc(data, m_data * sizeof(uint8_t));
        }
    }
};

void trim_array_tags(opts_t *op, bam1_t *b) {
    return;
}

static int trim_ns(bam1_t *b, void *data) {
    opts_t *op((opts_t *)data);
    int tmp;
    uint8_t *const seq(bam_get_seq(b));
    uint32_t *const cigar(bam_get_cigar(b));
    op->n_cigar = b->core.n_cigar;
    op->resize(b->l_data); // Make sure it's big enough to hold everything.
    memcpy(op->data, b->data, b->core.l_qname);

    // Get #Ns at the beginning 
    for(tmp = 0; bam_seqi(seq, tmp) == dlib::htseq::HTS_N; ++tmp);
    const int n_start(tmp);

    if(tmp == b->core.l_qseq - 1) // all bases are N -- garbage read
         return 1;

    // Get #Ns at the end
    for(tmp = b->core.l_qseq - 1; bam_seqi(seq, tmp) == dlib::htseq::HTS_N; --tmp);
    const int n_end(b->core.l_qseq - 1 - tmp);

    // Get new length for read
    const int final_len(b->core.l_qseq - n_end - n_start);
    if(final_len < op->min_trimmed_len) // Too short.
        return 1;

    if((tmp = bam_cigar_oplen(cigar[b->core.n_cigar - 1]) - n_end) == 0) {
        --op->n_cigar;
    } else {
        cigar[b->core.n_cigar - 1] = bam_cigar_gen(tmp, BAM_CSOFT_CLIP);
    }

    // Get new n_cigar. 
    if((tmp = bam_cigar_oplen(*cigar) - n_start) == 0) {
        --op->n_cigar;
        memcpy(op->data + b->core.l_qname, cigar + 1, op->n_cigar << 2); // << 2 for 4 bit per cigar op
    } else {
        *cigar = bam_cigar_gen(tmp, BAM_CSOFT_CLIP);
        memcpy(op->data + b->core.l_qname, cigar, op->n_cigar << 2);
    }
    uint8_t *opseq(op->data + b->core.l_qname + (op->n_cigar << 2));
    for(tmp = 0; tmp < final_len >> 1; ++tmp)
        opseq[tmp] = ((bam_seqi(seq, tmp * 2 + n_start) << 2) | (bam_seqi(seq, tmp * 2 + n_start + 1)));

    // Copy in qual and all of aux.
    tmp = bam_get_l_aux(b);
    memcpy(opseq + ((final_len + 1) >> 1), bam_get_qual(b) + n_start, final_len + tmp);
    // Switch data strings
    std::swap(op->data, b->data);
    b->core.n_cigar = op->n_cigar;
    b->l_data = b->core.l_qname + (op->n_cigar << 2) + ((final_len + 1) >> 1) + final_len + tmp;
    trim_array_tags(op, b);
    b->core.l_qseq = final_len;
    return 0;
}

int main(int argc, char *argv[]) {
    if(argc < 3)
        return usage(argv);

    if(strcmp(argv[1], "--help") == 0)
        return usage(argv, EXIT_SUCCESS);

    int c;
    char out_mode[4] = "wb";
    opts_t opts{0};
    while((c = getopt(argc, argv, "m:l:h?")) > -1) {
        switch(c) {
        case 'm': opts.min_trimmed_len = (uint32_t)atoi(optarg); break;
        case 'l': out_mode[2] = atoi(optarg) % 10 + '0'; break;
        case 'h': case '?': return usage(argv, EXIT_SUCCESS);
        }
    }
    if(argc - 2 != optind)
        LOG_EXIT("Required: precisely two positional arguments (in bam, out bam).\n");

    // Actually this function. You can't really apply a null function....
    dlib::BamHandle inHandle(argv[optind]);
    dlib::BamHandle outHandle(argv[optind + 1], inHandle.header, out_mode);
    dlib::abstract_single_iter(inHandle.fp, inHandle.header, outHandle.fp,
                               &trim_ns, (void *)&opts);
    return EXIT_SUCCESS;
}
