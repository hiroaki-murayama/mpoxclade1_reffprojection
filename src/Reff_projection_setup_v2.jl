sorted_keys = sort(collect(keys(all_susceptibilities)), by = x -> x[1], rev = true)
years = [2010, 2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060] 

# all contact
first_key = sorted_keys[1]
susceptibilities = all_susceptibilities[first_key] 
# obtain R₀ with respect to each year 
for year in years
    susceptibility = susceptibilities[year]
    zmb2015_24_fit.all[2].susceptibility[1] = susceptibility
    eigval = dominanteigval(zmb2015_24_fit.all[2])
    push!(dominant_eigvals_all, eigval)
end

# %%
# physical contact
second_key = sorted_keys[2]
susceptibilities = all_susceptibilities[second_key]
# obtain R₀ with respect to each year
for year in years
    susceptibility = susceptibilities[year]
    zmb2015_24_fit.phys[2].susceptibility[1] = susceptibility
    eigval = dominanteigval(zmb2015_24_fit.phys[2])
    push!(dominant_eigvals_phys, eigval)
end

# %%
# home contact
third_key = sorted_keys[3]
susceptibilities = all_susceptibilities[third_key]
# obtain R₀ with respect to each year
for year in years
    susceptibility = susceptibilities[year]
    zmb2015_24_fit.home[2].susceptibility[1] = susceptibility
    eigval = dominanteigval(zmb2015_24_fit.home[2])
    push!(dominant_eigvals_home, eigval)
end

# %%
# physical home contact
forth_key = sorted_keys[4]
susceptibilities = all_susceptibilities[forth_key]
# obtain R₀ with respect to each year
for year in years
    susceptibility = susceptibilities[year]
    zmb2015_24_fit.physhome[2].susceptibility[1] = susceptibility
    eigval = dominanteigval(zmb2015_24_fit.physhome[2])
    push!(dominant_eigvals_physhome, eigval)
end
