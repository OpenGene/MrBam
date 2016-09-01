using PyCall

@pyimport pysam

"""
get all reads of chr:pos
"""
function get_reads(p)
    samfile = pysam.AlignmentFile(o.reads, "rb")

    chr, pos = split(p, ':')
    pos = parse(Int, pos)

    @task for pileup in samfile[:pileup](chr, pos, pos+1)
        o.verbos && println("coverage at base $(pileup[:pos]) = $(pileup[:n])")

        pileup[:pos] != pos && continue

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
                alignment[:reference_end]
            )
        end

        return samfile[:close]()
    end
end
