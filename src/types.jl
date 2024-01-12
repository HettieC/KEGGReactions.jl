"""
$(TYPEDEF)
A struct for storing KEGG enzyme information.
$(FIELDS)
"""
mutable struct KEGGEnzyme
    id::String
    name::Maybe{String}
    reactions::Maybe{Vector{String}} #one enzyme can have multiple reactions
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
    name::Maybe{String}
    ec::Maybe{Vector{String}} # multiple ECs can be assigned to a single reaction
    pathway::Maybe{Vector{String}} #multiple pathways possible
    rhea::Maybe{Int64} 
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
    name::Maybe{String}
    charge::Maybe{Int64}
    formula::Maybe{String}
    dblinks::Maybe{Vector{String}}
end