#!/usr/bin/env julia

using ArgParse
using OhMyJulia

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
        "--ref", "-r"
            help = "reference file. If not provided, consider the mode of reads as the ref"
        "--verbos", '-v'
            help = "output debug info"
            nargs = 0
        "--quality", "-q"
            help = "drop bases whose qulity is less than this"
            arg_type = Int
            default = 20
    end

    parse_args(s) |> to_module
end
