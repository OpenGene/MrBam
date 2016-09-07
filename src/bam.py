from pysam import AlignmentFile

def get_reads(o, file, chr, pos):
    "get all reads covers chr:pos"

    samfile = AlignmentFile(file, "rb")

    pos = int(pos) - 1 # 0-based leftmost coordinate

    for read in samfile.fetch(chr, pos, pos+1):
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
            read.reference_length,
            read.next_reference_start,
            abs(read.template_length),
            read.is_reverse,
            read.is_paired and not read.mate_is_unmapped
        )

    samfile.close()
