# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

using Pkg
Pkg.activate("../")

# ### Preparation

# load functions
include("../src/eigen.jl");
include("../src/eigen_setup.jl");
include("../src/eigen_output.jl")

gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=8,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300))

# new
tshuapaplot = plot(tshuapa_h2hag,color=:black)
tshuapa2015_fit = output_fit(
    tshuapa_h2hag,
    zmb_skeleton = zmb2015,
    drc_skeleton = drc2015,
    dataplots = tshuapaplot
    );
zmb2015_fit = tshuapa2015_fit.zmb_fit;
endemicplot = plot(plot(drc_endemic_ag),ylim=(0,0.35),xtickfontsize=9)
endemic2015_24_fit = output_fit(
    [tshuapa_h2hag,drc_endemic_ag];
    zmb_skeleton = [zmb2015,zmb2024],
    drc_skeleton = [drc2015,drc2024],
    dataplots = [tshuapaplot,endemicplot]
    );
zmb2015_24_fit = endemic2015_24_fit.zmb_fit
drc2015_24_fit = endemic2015_24_fit.drc_fit;

# ### Demography

function calculate_age_ratios(countries::Vector{String}, years::AbstractVector{Int}, base_year::Int, directory::String)
    country_ratios = Dict{String, Vector{Float64}}()
    for country in countries
        for year in years
            # read csv automatically
            filepath = joinpath(directory, "$(country)-$(year).csv")
            df = CSV.read(filepath, DataFrame)
            # modified age tag
            age_ranges = [split(replace(String(age), "+" => ""), "-") for age in df.Age]
            age_min = parse.(Int, first.(age_ranges))
            age_max = [x == "" ? typemax(Int) : parse(Int, x) for x in last.(age_ranges)]
            # compute the proportion of 1980 cohort within 50+yo population in each year
            total_population = df.M .+ df.F
            x_age = 55 + (year - base_year)
            total_50_plus = sum(total_population[findall(x -> x >= 50, age_min)])
            total_x_age_plus = sum(total_population[findall(x -> x >= x_age, age_min)])
            ratio = total_x_age_plus / total_50_plus
            if !haskey(country_ratios, country)
                country_ratios[country] = []
            end
            push!(country_ratios[country], ratio)
        end
    end
    
    return country_ratios
end
# data
countries = ["United Kingdom", "Taiwan", "D.R. Congo"]
years = 2035:5:2060
base_year = 2035
directory = "../data/demography"
ratios = calculate_age_ratios(countries, years, base_year, directory)
println(ratios)

# +
# DRC
drc_ratios = ratios["D.R. Congo"]
println("D.R. Congo: ", drc_ratios)

# UK
uk_ratios = ratios["United Kingdom"]
println("United Kingdom: ", uk_ratios)

# Taiwan
taiwan_ratios = ratios["Taiwan"]
println("Taiwan: ", taiwan_ratios)
# -

# ### DRC

include("../src/Reff_projection_utils_v2.jl");

# +
years = [2010, 2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060]
s_infant_values = [2.2, 1.9, 1.5, 1.3]
s_vax_values = [0.26, 0.28, 0.22, 0.23]

all_susceptibilities = Dict{Tuple{Float64, Float64}, Dict}()
# get the correct combination of susceptibility by age group
for (s_infant, s_vax) in zip(s_infant_values, s_vax_values)
    s_partvax = (1 + s_vax) / 2  
    susceptibilities = create_susceptibilities_for_years(years, s_infant, s_partvax, s_vax, drc_ratios)
    all_susceptibilities[(s_infant, s_vax)] = susceptibilities
end

for (key, susceptibilities) in all_susceptibilities
    println("For s_infant = $(key[1]), s_vax = $(key[2]):")
    println(susceptibilities)
end
# -

dominant_eigvals_all = []
dominant_eigvals_phys = []
dominant_eigvals_home = []
dominant_eigvals_physhome = []
include("../src/Reff_projection_setup_v2.jl")

# all contact
println("Dominant Eigenvalues: ", dominant_eigvals_all)
println("Used combination: s_infant = $(first_key[1]), s_vax = $(first_key[2])")

# physical contact
println("Dominant Eigenvalues: ", dominant_eigvals_phys)
println("Used combination: s_infant = $(second_key[1]), s_vax = $(second_key[2])")

# home contact
println("Dominant Eigenvalues: ", dominant_eigvals_home)
println("Used combination: s_infant = $(third_key[1]), s_vax = $(third_key[2])")

# physical home contact
println("Dominant Eigenvalues: ", dominant_eigvals_physhome)
println("Used combination: s_infant = $(forth_key[1]), s_vax = $(forth_key[2])")

# +
# Rₑ projection with community cotacts only
# model weights
w1 = 0.22360
w2 = 0.38384
w3 = 0.21433
w4 = 0.17816

β_all = 0.82/dominant_eigvals_all[2]
β_phys = 0.82/dominant_eigvals_phys[2]
β_home = 0.82/dominant_eigvals_home[2]
β_physhome = 0.82/dominant_eigvals_physhome[2]

# model averaging
weighted_avg = []
for i in 1:length(dominant_eigvals_all)
    weighted_value = w1 * β_all * dominant_eigvals_all[i] + w2 * β_phys * dominant_eigvals_phys[i] + w3 * β_home * dominant_eigvals_home[i] + w4 * β_physhome * dominant_eigvals_physhome[i]
    push!(weighted_avg, weighted_value)
end
println("Weighted averages: ", weighted_avg)
# -

dominant_eigvals_drc=weighted_avg
println("R₀ (DRC): ", dominant_eigvals_drc)

β_all

# ### POLYMOD UK

# +
using RCall

@rimport socialmixr as smr

function get_contact_matrix(country::String, age_limits::Vector{Int})
    @rput country age_limits 
    R"""
    result <- socialmixr::contact_matrix(
        socialmixr::polymod, 
        countries = country, 
        age.limits = age_limits
    )
    """
    contact_matrix_uk = rcopy(R"result$matrix") 
    participants = rcopy(R"result$participants")  
    return contact_matrix_uk, participants
end

# get uk contact matrix
country = "United Kingdom"
age_limits = [0, 5, 10, 15, 20, 30, 40, 50]
contact_matrix_uk, participants = get_contact_matrix(country, age_limits)

println("Contact Matrix: ", contact_matrix_uk)
println("Participants: ", participants)


# +
all_susceptibilities = Dict{Tuple{Float64, Float64}, Dict}()
# get the correct combination of susceptibility by age group
for (s_infant, s_vax) in zip(s_infant_values, s_vax_values)
    s_partvax = (1 + s_vax) / 2  
    susceptibilities = create_susceptibilities_for_years(years, s_infant, s_partvax, s_vax, uk_ratios)
    all_susceptibilities[(s_infant, s_vax)] = susceptibilities
end

for (key, susceptibilities) in all_susceptibilities
    println("For s_infant = $(key[1]), s_vax = $(key[2]):")
    println(susceptibilities)
end
sorted_keys = sort(collect(keys(all_susceptibilities)), by = x -> x[1], rev = true)
years = [2010, 2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060] 

# all contact
first_key = sorted_keys[1]
susceptibilities = all_susceptibilities[first_key] 

# +
dominant_eigvals_uk = compute_dominant_eigenvalues(contact_matrix_uk, susceptibilities, years, β_all)

println("Rₑ (UK): ", dominant_eigvals_uk)
# -

maximum(abs.(eigvals(contact_matrix_uk))) 

# ### Taiwan

# +
using RCall

R"taiwan_survey <- readRDS('../data/taiwan.rds')"

@rimport socialmixr as smr

function get_contact_matrix(age_limits::Vector{Int})
    @rput country age_limits 
    R"""
    result <- socialmixr::contact_matrix(
        survey = taiwan_survey,  
        age.limits = age_limits
    )
    """
    contact_matrix = rcopy(R"result$matrix") 
    participants = rcopy(R"result$participants")  
    return contact_matrix, participants
end

# get thailand contact matrix
age_limits = [0, 5, 10, 15, 20, 30, 40, 50]
contact_matrix_taiwan, participants = get_contact_matrix(age_limits)

println("Contact Matrix: ", contact_matrix_taiwan)
println("Participants: ", participants)
# -

contact_matrix_taiwan

# +
all_susceptibilities = Dict{Tuple{Float64, Float64}, Dict}()
# get the correct combination of susceptibility by age group
for (s_infant, s_vax) in zip(s_infant_values, s_vax_values)
    s_partvax = (1 + s_vax) / 2  
    susceptibilities = create_susceptibilities_for_years(years, s_infant, s_partvax, s_vax, taiwan_ratios)
    all_susceptibilities[(s_infant, s_vax)] = susceptibilities
end

for (key, susceptibilities) in all_susceptibilities
    println("For s_infant = $(key[1]), s_vax = $(key[2]):")
    println(susceptibilities)
end
sorted_keys = sort(collect(keys(all_susceptibilities)), by = x -> x[1], rev = true)
years = [2010, 2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060] 

# all contact
first_key = sorted_keys[1]
susceptibilities = all_susceptibilities[first_key] 

# +
dominant_eigvals_taiwan = compute_dominant_eigenvalues(contact_matrix_taiwan, susceptibilities, years, β_all)

println("Rₑ (Taiwan): ", dominant_eigvals_taiwan)
# -

# ### Plot

# +
default_color = palette(:auto)[1]
colours = [:royalblue, :firebrick, :forestgreen]
#[:royalblue, 1, 2, :firebrick, RGBA(red(default_color) * 0.8, green(default_color) * 0.8, blue(default_color) * 0.8, 1.0), RGBA(65/255 * 0.8, 105/255 * 0.8, 225/255 * 0.8, 1.0)]  
bottom_margin = 10 * Plots.PlotMeasures.mm

plot([2015, 2015], [0.79, 0.85], color=:black, xlabel="year", ylabel="reproduction number", label="", ylim=(0,1.3), title="", size=(600,300), bottom_margin=bottom_margin)
scatter!([2015], [0.82], color=:black, marker=:black, markersize=4, markerstrokewidth=0, label=nothing)
plot!([2013, 2017], [0.82, 0.82], color=:black, label=nothing)
annotate!(2015, 0.7, text("Tshuapa", "Helvetica", :black, 7))

plot!(years, dominant_eigvals_drc, color=colours[1], label="DRC", linewidth=2)
plot!(years, dominant_eigvals_uk, color=colours[2], label="UK", linewidth=1.2)
plot!(years, dominant_eigvals_taiwan, color=colours[3], label="Taiwan", linewidth=1.2)

hline!([1], color=RGBA(0.5, 0.5, 0.5, 0.5), linestyle=:dash, label="")
