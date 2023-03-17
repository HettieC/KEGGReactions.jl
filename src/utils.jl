"""
$(TYPEDSIGNATURES)
Get the whole enzyme list from KEGG 
"""
function get_ec_req(;should_cache = true)
    _is_cached("ec_list", rid) &&
        return _get_cache("reaction_metabolites", rid)

    req = nothing
    # make sure the website is still valid 
    try 
        req = HTTP.request("GET","https://rest.kegg.jp/list/enzyme")
    catch
        req = nothing 
        println("KEGG API has been changed. Check the KEGG website and update the function ```update_ec```.")
    end

    return [String(x) for x in split(String(req.body),"\n")]
end

"""
$(TYPEDSIGNATURES)
Get the updated enzyme numbers, with the ```ec_req_lines``` input for the function being the output of get_ec_req().
"""
function update_ec(ec::String,ec_req_lines)
    for ln in ec_req_lines[1:end-1]
        ec_info = split(ln;limit=2)
        entry = split(ec_info[1],"ec:")[2]
        info = String(strip(ec_info[2]))
        if entry == ec 
            # get the new ec numbers
            if startswith(info,"Transferred to ") 
                x = split(info,"Transferred to ")[2]
                new_ecs = [String(y) for y in split(x," and ")]
                return vcat([update_ec(y,ec_req_lines) for y in new_ecs]...)
            # deleted entries    
            elseif info == "Deleted entry"
                return nothing      
            else # standard entries
                return [ec]
            end
        end
    end
end

"""
$(TYPEDSIGNATURES)
Get the name of and reactions associated with an EC number from the KEGG API.
"""
function get_kegg_ec(ec::String)
    req = nothing 
    try 
        req = HTTP.request("GET","https://rest.kegg.jp/get/$ec")
    catch
        req = nothing 
        print("No entry matching this id: $ec")
        return nothing
    end
    out = Dict{String,Any}()
    lines = split(String(req.body),"\n")
    if split(lines[1])[4] != "Enzyme"
        throw(error("Entry $ec not an enzyme"))
    else
        previous = lines[1]
        for ln in lines
            if startswith(ln,"SYSNAME")
                out["ec_name"] = strip(split(ln;limit=2)[2])
            elseif startswith(ln,"ALL_REAC")
                out["rxns"] = [String(split(x,";")[1]) for x in split(ln) if startswith(x,"R")]
            elseif startswith(ln," ") && startswith(previous,"ALL_REAC")
                append!(out["rxns"],[String(split(x,";")[1]) for x in split(ln) if startswith(x,"R")])
            end
            previous = ln
        end
    end
    return out
end

"""
$(TYPEDSIGNATURES)
Get the reaction name, stoichiometry, and database cross references
"""
function get_kegg_rxn(rxn_id::String)
    if contains(rxn_id,"(G)")
        return nothing
    else
        req = nothing 
        try 
            req = HTTP.request("GET","https://rest.kegg.jp/get/$rxn_id")
        catch
            req = nothing 
            print("No entry matching this id: $rxn_id")
        end
        out = Dict{String,Any}()
        lines = split(String(req.body),"\n")
        stoich = Dict{String,Int64}()
        if split(lines[1])[3] != "Reaction"
            throw(error("Entry $rxn_id not a reaction"))
        else
            for ln in lines 
                if startswith(ln,"NAME")
                    out["rxn_name"] = strip(String(split(ln;limit=2)[2]))
                elseif startswith(ln,"EQUATION")
                    subs_prods = split( strip(split(ln,"EQUATION")[2]) , "<=>")
                    subs = [strip(x) for x in split(subs_prods[1]," + ")]
                    prods = [strip(x) for x in split(subs_prods[2]," + ")]
                    for s in subs 
                        if startswith(s,"C")
                            stoich[s] = -1
                        else
                            x = split(s)
                            stoich[String(x[2])] = isnothing(tryparse(Int64,x[1])) ? 0 : -tryparse(Int64,x[1])
                        end
                    end
                    for p in prods
                        if startswith(p,"C")
                            stoich[p] = 1
                        else
                            x = split(p)
                            stoich[String(x[2])] = isnothing(tryparse(Int64,x[1])) ? 0 : tryparse(Int64,x[1])
                        end
                    end
                    out["stoich"] = stoich
                elseif startswith(ln,"PATHWAY")
                    out["pathway"] = [String(strip(split(ln;limit=2)[2]))]
                elseif startswith(strip(ln),"rn")
                    push!(out["pathway"],String(strip(ln)))
                elseif contains(ln,"RHEA:") && !haskey(out,"RHEA")
                    out["RHEA"] = [parse(Int64,split(ln)[3])]
                elseif contains(ln,"RHEA:")
                    append!(out["RHEA"],parse(Int64,split(ln)[3]))
                end
            end
        end
        return out
    end
end

"""
$(TYPEDSIGNATURES)
Get the name and formula of a compound.
"""
function get_kegg_cmpd(cmpd_id::String)
    req = nothing 
    try 
        req = HTTP.request("GET","https://rest.kegg.jp/get/$cmpd_id")
    catch
        req = nothing 
        print("No entry matching this id: $cmpd_id")
    end
    out = Dict{String,Any}()
    out["dblinks"]=String[]
    lines = split(String(req.body),"\n")
    if split(lines[1])[3] != "Compound"
        throw(error("Entry $cmpd_id not a compound"))
    else
        for ln in lines
            if startswith(ln,"NAME")
                out["name"] = String(strip(split(ln;limit=2)[2]))
            elseif startswith(ln,"FORMULA")
                out["formula"] = String(strip(split(ln;limit=2)[2]))
            elseif contains(ln,"PubChem:")
                push!(out["dblinks"],String(strip(split(ln,"DBLINKS")[end])))
            elseif contains(ln,"ChEBI:")
                push!(out["dblinks"],String(strip(split(ln,"DBLINKS")[end])))
            end
        end
    end
    return out
end

"""
$(TYPEDSIGNATURES)
Make a dictionary of ec against its kegg reactions with all necessary info.
"""
function get_dataset_kegg()
    if !ispath("data/model_building/kegg_dataset.json")
        # update the ec numbers from uniprot
        ec_req_lines = get_ec_req()
        ec_list = String[]
        redundant_ecs = String[]
        open("data/model_building/ec_list.txt") do io 
            for ln in eachline(io) 
                if !isnothing(update_ec(ln,ec_req_lines))
                    append!(ec_list,update_ec(ln,ec_req_lines))
                else        
                    push!(redundant_ecs,ln)
                end
            end
        end
        unique!(ec_list)
        unique!(redundant_ecs)

        # get the kegg entries for each ec
        ec_kegg = Dict{String,Dict{String,Any}}()
        for ec in ec_list 
            ec_kegg[ec] = get_kegg_ec(ec)
        end

        # get kegg entries for each reaction 
        ec_rxns = Dict{String,Dict{String,Dict{String,Any}}}()
        for (k,v) in ec_kegg
            println(k)
            if haskey(v,"rxns")
                for rxn in v["rxns"]
                    if isnothing(get_kegg_rxn(rxn))
                        continue
                    elseif !haskey(ec_rxns,k)
                        ec_rxns[k] = Dict(rxn => get_kegg_rxn(rxn))
                    else
                        ec_rxns[k][rxn] = get_kegg_rxn(rxn)
                    end
                end
            end 
        end
        return ec_rxns
    else
        return JSON.parse(open("data/model_building/kegg_dataset.json"))
    end
end



"""
$(TYPEDSIGNATURES)
Make list of lists of the reaction directions, using the rhea-directions.tsv file
# import the rhea-directions.tsv file so that I only use master reactions
"""
function parse_rhea_directions(path_to_rhea_directions::String)
    directions = Vector{Vector{Int64}}()
    open(path_to_rhea_directions) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline=false;continue) 
            dirs = split(ln, "\t")
            push!(directions,[parse(Int64,x) for x in dirs])
        end
    end 
    return directions
end

"""
$(TYPEDSIGNATURES)
"""
function add_reaction_from_kegg!(
    model;
    kegg_id::String,
    name = "",
    isozymes = nothing,
    subsystem = nothing,
    ec = nothing,
    iso_confidence = nothing
)
    # add gene to model 
    if !isnothing(isozymes)
        for isozyme in isozymes 
            for (_, gid) in isozyme
                gid ∉ genes(model) && add_gene!(model, Gene(gid))
            end
        end    
    end

    # get reaction components
    rxn = get_kegg_rxn(kegg_id)
    metabolite_dict = haskey(rxn,"stoich") ? rxn["stoich"] : nothing
    name = haskey(rxn,"rxn_name") ? rxn["rxn_name"] : nothing
    rhea_id = haskey(rxn,"RHEA") ? ["$x" for x in rxn["RHEA"]] : nothing
    annotations = Dict(
        "EC" => [ec],
        "Iso_Confidence" => ["$iso_confidence"],
    )
    if !isnothing(rhea_id)
        annotations["RHEA"] = rhea_id
    end
    # add reaction to model 
    isos = isnothing(isozymes) ?  nothing : [Isozyme(gene_product_stoichiometry = Dict(gid => num for (num, gid) in iso)) for iso in isozymes]
    add_reaction!(
        model,
        Reaction(
            kegg_id;
            name,
            metabolites = metabolite_dict,
            gene_associations = isos,
            subsystem,
            annotations,
            )
        )

    # add metabolites to model 
    for (met,_) in metabolite_dict
        if met ∉ metabolites(model)
            kegg_met = get_kegg_cmpd(met)
            name = haskey(kegg_met,"name") ? kegg_met["name"] : nothing 
            formula = haskey(kegg_met,"formula") ? kegg_met["formula"] : nothing 
            dblinks = haskey(kegg_met,"dblinks") ? kegg_met["dblinks"] : nothing
            add_metabolite!(
                model,
                Metabolite(
                    met;
                    name,
                    formula,
                    annotations = Dict(
                        "dblinks" => dblinks
                    )
                )
            )
        end
    end
end
