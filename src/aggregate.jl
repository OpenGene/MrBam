"aggregate reads by startpos, endpos and base"
function aggregate_reads(reads)
    pair_dict     = Dict{AbstractString, Vector}() # name -> reads
    unique_pairs  = Dict{Tuple{Int, Int, Char}, Bytes}() # start, end, base -> quality
    unique_single = Dict{Tuple{Int, Int, Char}, Bytes}() # start, end, base -> quality

    nsum,  nerror    = 0, 0 # depth, errors (such as 3 reads share the same name)
    nlowq, ninconsis = 0, 0 # low quality bases, inconsistent pairs (count on reads)

    # 1. find pairs
    for read in reads
        name, info = car(read), cdr(read)

        if name in keys(pair_dict)
            if length(pair_dict[name]) != 1 # 0 or 2
                println(STDERR, "more than 2 reads share the same name; all droped.")
                nerror += length(pair_dict[name]) + 1
                pair_dict[name] = []
            else
                push!(pair_dict[name], info)
            end
        else
            pair_dict[name] = [info]
        end

        nsum += 1
    end

    # 2. find unique pair/single s
    for pair in values(pair_dict)
        if length(pair) == 1 # single
            base, qual, startpos, endpos = car(pair)

            if qual <= o.quality
                nlowq += 1
                continue
            end

            if (startpos, endpos, base) in keys(unique_single)
                push!(unique_single[(startpos, endpos, base)], Byte(qual))
            else
                unique_single[(startpos, endpos, base)] = Byte[qual]
            end
        else # pair
            r1, r2 = pair
            base1, qual1, startpos1, endpos1 = r1
            base2, qual2, startpos2, endpos2 = r2

            if qual1 <= o.quality || qual2 <= o.quality
                nlowq += 2
                continue
            end

            if base1 != base2
                ninconsis += 2
                continue
            end

            startpos = min(startpos1, startpos2)
            endpos   = max(endpos1,   endpos2)
            qual     = max(qual1,     qual2)
            base     = base1

            if (startpos, endpos, base) in keys(unique_pairs)
                push!(unique_pairs[(startpos, endpos, base)], Byte(qual))
            else
                unique_pairs[(startpos, endpos, base)] = Byte[qual]
            end
        end
    end

    unique_pairs, unique_single, nsum, nerror, nlowq, ninconsis
end
