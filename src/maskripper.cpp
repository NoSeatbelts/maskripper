#include <getopt.h>
#include "dlib/bam_util.h"

int usage(char **argv, int retcode=EXIT_FAILURE) {
    fprintf(stderr, "MASKRIPPER version %s.\n"
                    "Usage: maskripper <opts> in.bam out.bam\n"
                    "Flags:\n-m: Minimum trimmed read length. Default: 0.\n"
                    "-l: output compression level. Default: 6.\n"
                    "-S: output in sam format.\n"
                    "Use - for stdin or stdout.\n", MASKRIPPER_VERSION);
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

static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}


static inline uint8_t *skip_aux(uint8_t *s)
{
    int size = aux_type2size(*s); ++s; // skip type
    uint32_t n;
    switch (size) {
    case 'Z':
    case 'H':
        while (*s) ++s;
        return s + 1;
    case 'B':
        size = aux_type2size(*s); ++s;
        memcpy(&n, s, 4); s += 4;
        return s + size * n;
    case 0:
        abort();
        break;
    default:
        return s + size;
    }
}

uint8_t *aux_aux_get(uint8_t *s, const char tag[2], uint8_t *end)
{
    LOG_DEBUG("Trying\n");
    int y = tag[0]<<8 | tag[1];
    while (s < end) {
        int x = (int)s[0]<<8 | s[1];
        s += 2;
        if (x == y) return s;
        s = skip_aux(s);
    }
    return 0;
}


template<typename T>
inline void trim_array_tag(bam1_t *b, const char *tag, const int n_start, const int n_end, const int final_len) {
    T tmp[300];
    if((n_start + n_end) == 0) return;
    LOG_DEBUG("Trying for tag %s.\n", tag);
    uint8_t *data = bam_aux_get(b, tag);
    uint8_t *td = data + 4;
    LOG_DEBUG("Len: %i.\n", *(int *)td);
    LOG_DEBUG("Type: %c, %i.\n", (char)*(data + 3), *(data + 3));
    LOG_DEBUG("n_start, final_len: %i, %i.\n", n_start, final_len);
    for(int i = 0; i < n_start + n_end + final_len; ++i) fprintf(stderr, ",%u", ((int *)(data + 8))[i]);
    memcpy(tmp, data + 8 + n_start * sizeof(T), final_len * sizeof(T)); // 8 for the tag, len, type.
    LOG_DEBUG("Old len: %i. New len: %i.\n", final_len + n_start + n_end, final_len);
    bam_aux_del(b, data);
    assert((data = bam_aux_get(b, tag)) == nullptr);
    for(int i = 0; i < final_len; ++i) fprintf(stderr, ",%u", tmp[i]);
    dlib::bam_aux_array_append(b, tag, 'I', sizeof(uint32_t), final_len, (uint8_t *)tmp);
}

inline void trim_array_tags(bam1_t *b, const int n_start, const int n_end, const int final_len) {
    LOG_DEBUG("Getting PV.\n");
    trim_array_tag<uint32_t>(b, "PV", n_start, n_end, final_len);
    LOG_DEBUG("Getting FA.\n");
    trim_array_tag<uint32_t>(b, "FA", n_start, n_end, final_len);
}


static int trim_ns(bam1_t *b, void *data) {
    // Currently passes all reads to keep balanced pairs. TODO: rewrite for pairs of reads and filter if both fail.
    opts_t *op((opts_t *)data);
    assert(bam_aux_get(b, "PV"));
    std::vector<uint8_t> aux(bam_get_aux(b), bam_get_aux(b) + bam_get_l_aux(b));
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
         return 0; // Currently outputting to avoid 

    // Get #Ns at the end
    for(tmp = b->core.l_qseq - 1; bam_seqi(seq, tmp) == dlib::htseq::HTS_N; --tmp);
    const int n_end(b->core.l_qseq - 1 - tmp);

    // Get new length for read
    const int final_len(b->core.l_qseq - n_end - n_start);
    if(final_len < op->min_trimmed_len) // Too short.
        return 0;

    if(n_end) {
        if((tmp = bam_cigar_oplen(cigar[b->core.n_cigar - 1]) - n_end) == 0) {
            LOG_DEBUG("Entire cigar operation is the softclip. Decrease the number of new cigar operations.\n");
            --op->n_cigar;
        } else {
            LOG_DEBUG("Updating second cigar operation in-place.\n");
            cigar[b->core.n_cigar - 1] = bam_cigar_gen(tmp, BAM_CSOFT_CLIP);
        }
    }

    // Get new n_cigar.
    if((tmp = bam_cigar_oplen(*cigar) - n_start) == 0) {
        --op->n_cigar;
        memcpy(op->data + b->core.l_qname, cigar + 1, op->n_cigar << 2); // << 2 for 4 bit per cigar op
    } else {
        if(n_start) *cigar = bam_cigar_gen(tmp, BAM_CSOFT_CLIP);
        memcpy(op->data + b->core.l_qname, cigar, op->n_cigar << 2);
    }
    uint8_t *opseq(op->data + b->core.l_qname + (op->n_cigar << 2)); // Pointer to the seq region of new data field.
    for(tmp = 0; tmp < final_len >> 1; ++tmp) {
        opseq[tmp] = (bam_seqi(seq, ((tmp << 1) + n_start)) << 4) | (bam_seqi(seq, (tmp << 1) + n_start + 1));
        assert(bam_seqi(opseq, tmp * 2) == bam_seqi(seq, (2 * tmp + n_start)));
        assert(bam_seqi(opseq, tmp * 2 + 1) == bam_seqi(seq, (2 * tmp + n_start + 1)));
    }
#if 0
    char tmpbuf[300];
    for(tmp = 0; tmp < final_len; ++tmp) {
        tmpbuf[tmp] = seq_nt16_str[bam_seqi(opseq, tmp)];
    }
    tmpbuf[tmp] = '\0';
    assert(strcmp(tmpbuf, "ANNCNNNGNNNNTNNNNCNNGGNNNNNNNNNCNNNNCNNNNNNNNAAAANNTNNNAAAAAAAAAAGAGAGAGGGAGAGAGACTATACACAGGCACCACCACATTTGGCTAATTTTT") == 0);
#endif

    // Copy in qual and all of aux.
    tmp = bam_get_l_aux(b);
    memcpy(opseq + ((final_len + 1) >> 1), bam_get_qual(b) + n_start, final_len + tmp);
    // Switch data strings
    std::swap(op->data, b->data);
    b->core.n_cigar = op->n_cigar;
    //b->l_data = b->core.l_qname + (op->n_cigar << 2) + ((final_len + 1) >> 1) + final_len + tmp;
    b->core.l_qseq = final_len;
    memcpy(bam_get_aux(b), aux.data(), aux.size());
    b->l_data = (bam_get_aux(b) - b->data) + aux.size();
    //trim_array_tags(b, n_start, n_end, final_len);
    //bam_aux_append(b, "NE", 'i', sizeof(int), (uint8_t *)&n_end);
    //bam_aux_append(b, "NS", 'i', sizeof(int), (uint8_t *)&n_start);
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
    while((c = getopt(argc, argv, "m:l:h?S")) > -1) {
        switch(c) {
        case 'm': opts.min_trimmed_len = (uint32_t)atoi(optarg); break;
        case 'l': out_mode[2] = atoi(optarg) % 10 + '0'; break;
        case 'S': sprintf(out_mode, "w"); break;
        case 'h': case '?': return usage(argv, EXIT_SUCCESS);
        }
    }
    if(argc - 2 != optind)
        LOG_EXIT("Required: precisely two positional arguments (in bam, out bam).\n");

    LOG_DEBUG("ZOUNDS\n");
    // Actually this function. You can't really apply a null function....
    dlib::BamHandle inHandle(argv[optind]);
    dlib::BamHandle outHandle(argv[optind + 1], inHandle.header, out_mode);
    dlib::abstract_single_iter(inHandle.fp, inHandle.header, outHandle.fp,
                               &trim_ns, (void *)&opts);
    return EXIT_SUCCESS;
}
