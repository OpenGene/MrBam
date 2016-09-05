import pysam

def get_reads(o, p):
    "get all reads covers chr:pos"

    samfile = pysam.AlignmentFile(o.reads, "rb")

    ch, pos = p.split(':')
    pos = int(pos) - 1 # 0-based leftmost coordinate

    for read in samfile.fetch(ch, pos, pos+1):
        aligned_pairs = read.get_aligned_pairs(matches_only=True)

        try:
            query_pos = next(qpos for qpos, rpos in aligned_pairs if rpos == pos)
        except:
            continue

        yield (
            read.query_name,
            read.query_sequence[query_pos],
            read.query_qualities[query_pos],
            read.reference_start,
            read.next_reference_start,
            read.template_length
        )

    samfile.close()
