





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
        req = HTTP.request("GET", "https://rest.kegg.jp/get/$ec")
    catch
        req = nothing
        print("No entry matching this id: $ec")
        return nothing
    end
    out = Dict{String,Any}()
    lines = split(String(req.body), "\n")
    if split(lines[1])[4] != "Enzyme"
        throw(error("Entry $ec not an enzyme"))
    else
        previous = lines[1]
        for ln in lines
            if startswith(ln, "SYSNAME")
                out["ec_name"] = strip(split(ln; limit = 2)[2])
            elseif startswith(ln, "ALL_REAC")
                out["rxns"] =
                    [String(split(x, ";")[1]) for x in split(ln) if startswith(x, "R")]
            elseif startswith(ln, " ") && startswith(previous, "ALL_REAC")
                append!(
                    out["rxns"],
                    [String(split(x, ";")[1]) for x in split(ln) if startswith(x, "R")],
                )
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
function get_kegg_rxn(rxn_id::String;cache=true, max_retries = 5)
    rxn = nothing
    if contains(rxn_id, "(G)")
        return nothing
    else
        if _is_cached("reaction",rxn_id)
            return _get_cache("reaction",rxn_id)
        else 
            retry_counter = 0 
            req = nothing
            while retry_counter <= max_retries
                retry_counter += 1
                try
                    req = HTTP.request("GET", "https://rest.kegg.jp/get/$rxn_id")
                catch
                    req = nothing
                    print("No entry matching this id: $rxn_id")
                end
                out = Dict{String,Any}()
                lines = split(String(req.body), "\n")
                stoich = Dict{String,Int64}()
                if split(lines[1])[3] != "Reaction"
                    throw(error("Entry $rxn_id not a reaction"))
                else
                    for ln in lines
                        if startswith(ln, "NAME")
                            out["name"] = strip(String(split(ln; limit = 2)[2]))
                        elseif startswith(ln, "EQUATION")
                            subs_prods = split(strip(split(ln, "EQUATION")[2]), "<=>")
                            subs = [strip(x) for x in split(subs_prods[1], " + ")]
                            prods = [strip(x) for x in split(subs_prods[2], " + ")]
                            for s in subs
                                if startswith(s, "C")
                                    stoich[s] = -1
                                else
                                    x = split(s)
                                    stoich[String(x[2])] =
                                        isnothing(tryparse(Int64, x[1])) ? 0 :
                                        -tryparse(Int64, x[1])
                                end
                            end
                            for p in prods
                                if startswith(p, "C")
                                    stoich[p] = 1
                                else
                                    x = split(p)
                                    stoich[String(x[2])] =
                                        isnothing(tryparse(Int64, x[1])) ? 0 : tryparse(Int64, x[1])
                                end
                            end
                            out["stoich"] = stoich
                        elseif startswith(ln, "ENZYME")
                            out["ec"] = [String(strip(split(ln; limit = 2)[2]))]
                        elseif startswith(ln, "PATHWAY")
                            out["pathway"] = [String(strip(split(ln; limit = 2)[2]))]
                        elseif startswith(strip(ln), "rn")
                            push!(out["pathway"], String(strip(ln)))
                        elseif contains(ln, "RHEA:") && !haskey(out, "RHEA")
                            out["RHEA"] = [String(split(ln)[3])]
                        elseif contains(ln, "RHEA:")
                            append!(out["RHEA"], String(split(ln)[3]))
                        end
                    end
                end
                rxn = KEGGReaction(;
                        id=rxn_id,
                        name = haskey(out, "name") ? out["name"] : nothing,
                        stoichiometry = haskey(out, "stoich") ? out["stoich"] : nothing,
                        ec = haskey(out, "ec") ? out["ec"] : nothing,
                        pathway = haskey(out, "pathway") ? out["pathway"] : nothing,
                        dblinks = haskey(out, "RHEA") ? Dict("RHEA" => out["RHEA"]) : nothing,
                    )
                if cache
                    _cache("reaction",rxn_id,rxn)
                end
            end
            return rxn
        end
    end
end

export get_kegg_rxn

"""
$(TYPEDSIGNATURES)
Get the name and formula of a compound.
"""
function get_kegg_met(met_id::String; cache=true, max_retries = 5)
    met = nothing
    if _is_cached("reaction_metabolites",met_id)
        return _get_cache("reaction_metabolites",met_id)
    else
        retry_counter = 0 
        req = nothing
        while retry_counter <= max_retries
            retry_counter += 1
            try
                req = HTTP.request("GET", "https://rest.kegg.jp/get/$met_id")
            catch
                req = nothing
                print("No entry matching this id: $met_id")
            end
            out = Dict{String,Any}()
            out["dblinks"] = Dict{String,Vector{String}}()
            lines = split(String(req.body), "\n")
            if split(lines[1])[3] != "Compound"
                throw(error("Entry $met_id not a compound"))
            else
                for ln in lines
                    if startswith(ln, "NAME")
                        out["name"] = String(strip(split(ln; limit = 2)[2]))
                    elseif startswith(ln, "FORMULA")
                        out["formula"] = String(strip(split(ln; limit = 2)[2]))

                    elseif contains(ln, "PubChem:")
                        out["dblinks"]["PubChem"] = [String(strip(split(ln, "PubChem: ")[2]))]
                    elseif contains(ln, "ChEBI:")
                        out["dblinks"]["ChEBI"] = [String(strip(split(ln, "ChEBI: ")[2]))]
                    end
                end
            end
            met = KEGGMetabolite(
                id = met_id,
                name = out["name"],
                formula = haskey(out, "formula") ? out["formula"] : nothing,
                dblinks = haskey(out, "dblinks") ? out["dblinks"] : nothing,
            )
            if cache 
                _cache("reaction_metabolites",met_id,met)
            end
        end
        return met
    end
end

export get_kegg_met

