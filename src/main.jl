#!/usr/bin/env julia

using ArgParse
using OhMyJulia
using ProgressMeter

include(rel"bam.jl")
include(rel"aggregate.jl")

const o = let
    s = ArgParseSettings()

    @add_arg_table s begin
        "query"
            help = "vcf file contains mutations to query"
            required = true
        "reads"
            help = "bam file contains origin reads info. There must be a corresponding .bai file in the same directory"
            required = true
        "--quality", "-q"
            help = "drop bases whose qulity is less than this"
            arg_type = Int
            default = 20
        "--ref", "-r"
            help = "reference file. If not provided, consider the mode of reads as the ref"
        "--verbos", "-v"
            help = "output debug info"
            nargs = 0
    end

    parse_args(s) |> to_module
end

open(o.query * ".readsinfo.txt", "w") do fout
    fout << "#Chr" << '\t' << "Position" << '\t' << "Unique_Pairs_Support_Alt" << '\t' << "Unique_Single_Support_Alt" << '\t'
    fout << "N_Reads" << '\t' << "N_Error" << '\t' << "N_Low_Quality" << '\t' << "N_Inconsistence_Pairs" << '\n'

    @showprogress 1 "wait..." for line in open(readlines, o.query)
        startswith(line, '#') && continue

        p = join(split(line, '\t')[[1,2]], ':')

        unique_pairs, unique_single, nsum, nerror, nlowq, ninconsis = p |> get_reads |> aggregate_reads

        refbase = let # TODO: query from ref file if provided, or parse and get from vcf file
            d = Dict('A'=>0, 'T'=>0, 'C'=>0, 'G'=>0)
            for x in (unique_single, unique_pairs), (startpos, endpos, base) in keys(x)
                d[base] += 1
            end
            max(d..., by=x->x.second).first
        end

        unique_pairs_support_alt  = count(x->x[3]!=refbase, keys(unique_pairs))
        unique_single_support_alt = count(x->x[3]!=refbase, keys(unique_single))

        fout << replace(p, ':', '\t', 1) << '\t' << unique_pairs_support_alt << '\t' << unique_single_support_alt << '\t'
        fout << nsum << '\t' << nerror << '\t' << nlowq << '\t' << ninconsis << '\n'
    end
end
