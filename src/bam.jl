using PyCall

@pyimport pysam

"""
get all reads covers chr:pos
"""
function get_reads_by_pileup(p)
    samfile = pysam.AlignmentFile(o.reads, "rb")

    chr, pos = split(p, ':')
    pos = parse(Int, pos) - 1 # 0-based leftmost coordinate

    @task for pileup in samfile[:pileup](chr, pos, pos+1)
        o.verbos && println("coverage at base $(pileup[:reference_pos]+1) = $(pileup[:n])")

        pileup[:reference_pos] != pos && continue

        for read in pileup[:pileups]
            if read[:is_del] == 1 || read[:is_refskip] == 1
                o.verbos && println(STDERR, "isdel or isrefskip")
                break
            end

            query_pos = read[:query_position] + 1 # julia use 1-based indexing
            alignment = read[:alignment]

            produce(
                alignment[:query_name],
                alignment[:query_sequence][query_pos],
                alignment[:query_qualities][query_pos],
                alignment[:reference_start],
                alignment[:next_reference_start],
                alignment[:template_length]
            )
        end

        return samfile[:close]()
    end
end

"""
get all reads covers chr:pos
"""
function get_reads_by_fetch(p)
    samfile = pysam.AlignmentFile(o.reads, "rb")

    chr, pos = split(p, ':')
    pos = parse(Int, pos) - 1 # 0-based leftmost coordinate

    @task begin
        for read in samfile[:fetch](chr, pos, pos+1)
            query_pos = let
                aligned_pairs = read[:get_aligned_pairs](matches_only=true)
                i = findfirst(x->cadr(x)==pos, aligned_pairs)
                i == 0 && continue
                car(aligned_pairs[i]) + 1 # julia use 1-based indexing
            end

            produce(
                read[:query_name],
                read[:query_sequence][query_pos],
                read[:query_qualities][query_pos],
                read[:reference_start],
                read[:next_reference_start],
                read[:template_length]
            )
        end

        return samfile[:close]()
    end
end

const get_reads = get_reads_by_fetch
