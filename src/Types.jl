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
    ec::Union{Vector{String},Nothing} # multiple ECs can be assigned to a single reaction
    pathway::Union{Vector{String},Nothing} #multiple pathways possible
    dblinks::Union{Dict{String,Vector{String}},Nothing}
end

"""
$(TYPEDEF)
A struct for storing KEGG metabolite information.
$(FIELDS)
"""
@kwdef mutable struct KEGGMetabolite
    id::String
    name::String
    formula::Union{String,Nothing}
    dblinks::Union{Dict{String,Vector{String}},Nothing}
end

end