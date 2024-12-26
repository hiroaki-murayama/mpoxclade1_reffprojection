function create_susceptibility_by_year(year::Int, value1::Float64, value2::Float64, value3::Float64, ageprop::Vector{Float64})
    # fixed the first four values
    base_susceptibility = [
        fill(value1), #0-4
        fill(1), #5-9
        fill(1), #10-14
        fill(1), #15-19
        fill(1) #20-29
    ]
    # cumputed the suscaptibility after 2040 with demographic projection used (proportion of 1980 cohort within 50+yo population in each year)
    value4 = value3*ageprop[1] + (1-ageprop[1]) # 2035
    value5 = value3*ageprop[2] + (1-ageprop[2]) # 2040
    value6 = value3*ageprop[3] + (1-ageprop[3]) # 2045
    value7 = value3*ageprop[4] + (1-ageprop[4]) # 2050
    value8 = value3*ageprop[5] + (1-ageprop[5]) # 2055
    value9 = value3*ageprop[6] + (1-ageprop[6]) # 2060
    # susceptibility values change depending on the year
    if year == 2010
        append!(base_susceptibility, [fill(value3), fill(value3), fill(value3)])
    elseif year == 2015
        append!(base_susceptibility, [fill(value2), fill(value3), fill(value3)])
    elseif year == 2020
        append!(base_susceptibility, [fill(1), fill(value3), fill(value3)])
    elseif year == 2025
        append!(base_susceptibility, [fill(1), fill(value2), fill(value3)])
    elseif year == 2030
        append!(base_susceptibility, [fill(1), fill(1), fill(value3)])
    elseif year == 2035
        append!(base_susceptibility, [fill(1), fill(1), fill(value4)])
    elseif year == 2040
        append!(base_susceptibility, [fill(1), fill(1), fill(value5)])
    elseif year == 2045
        append!(base_susceptibility, [fill(1), fill(1), fill(value6)])
    elseif year == 2050
        append!(base_susceptibility, [fill(1), fill(1), fill(value7)])
    elseif year == 2055
        append!(base_susceptibility, [fill(1), fill(1), fill(value8)])
    elseif year == 2060
        append!(base_susceptibility, [fill(1), fill(1), fill(value9)])
    else
        error("Invalid year. Provide a year between 2010 and 2060.")
    end
    return base_susceptibility
end

function create_susceptibilities_for_years(years::Vector{Int}, value1::Float64, value2::Float64, value3::Float64, ageprop::Vector{Float64})
    susceptibilities = Dict{Int, Vector}()
    for year in years
        susceptibilities[year] = create_susceptibility_by_year(year, value1, value2, value3, ageprop)
    end
    return susceptibilities
end

# function for computing dominant eigenvalue of ngm where contact matrix is exchangable
function compute_dominant_eigenvalues(contact_matrix::Matrix{Float64}, susceptibilities::Dict{Int64, Vector}, years::Vector{Int}, β::Float64)
    dominant_eigvals = []  
    for year in years
        σ = susceptibilities[year] 
        K = [β * σ[a] * contact_matrix[a, b] for a in 1:size(contact_matrix, 1), b in 1:size(contact_matrix, 2)]
        eigval = maximum(abs.(eigvals(K)))
        push!(dominant_eigvals, eigval)
    end
    return dominant_eigvals
end