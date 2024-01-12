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
mutable struct KEGGReaction
    id::Int64
    equation::String
    name::String
    ec::Vector{String} # multiple ECs can be assigned to a single reaction
    pathway::Vector{String} #multiple pathways possible
    rhea::Int64
    istransport::Bool
    isbalanced::Bool
end

KEGGReaction() = KEGGReaction(0, "", "", "", nothing, nothing, false, false)

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