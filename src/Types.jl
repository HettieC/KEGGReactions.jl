module Types

using DocStringExtensions

"""
$(TYPEDEF)
A struct for storing KEGG enzyme information.
$(FIELDS)
"""
mutable struct KEGGEnzyme
    id::String
    name::String
    reactions::Vector{String} #one enzyme can have multiple reactions
end


"""
$(TYPEDEF)
A struct for storing KEGG reaction information. Does not store the metabolite
information. 
$(FIELDS)
"""
@kwdef mutable struct KEGGReaction
    id::String
    name::String
    stoichiometry::Dict{String,Int64}
    ec::Vector{String} # multiple ECs can be assigned to a single reaction
    pathway::Vector{String} #multiple pathways possible
    dblinks::Dict{String,Vector{String}}
end

"""
$(TYPEDEF)
A struct for storing KEGG metabolite information.
$(FIELDS)
"""
struct KEGGMetabolite
    id::Int64
    name::String
    charge::Int64
    formula::String
    dblinks::Vector{String}
end

end