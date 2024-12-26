using JLD2
using CSV
using DataFrames
using StatsBase
using StatsPlots
using Colors

# utils
function collapseplot(pyramid::Pyramid; kwargs...)
    joinedpyramid=deepcopy(pyramid)
    joinedpyramid.cases = [vcat(joinedpyramid.cases...)]
    plot(joinedpyramid;kwargs...)
end

# plot functions
function Plots.plot(p::Pyramid; color = :black, ascertainment=[1], ylim=(0,0.5), kwargs...)
    cases_errors = [begin
            posterior_agedist = rand(Dirichlet(cases.+1),10000).-normalize(cases,1)
            low =  quantile.(collect.(eachrow(posterior_agedist)),0.025).|> abs
            upp =  quantile.(collect.(eachrow(posterior_agedist)),0.975).|> abs
            (low.*ascertainment,upp.*ascertainment)
        end
    for cases in p.cases]
    xrepeat = length(p.cases|>first)÷length(p.ageinterval)
    xticks = (1:length(p.cases|>first),repeat(makeagegrouplabel(p.ageinterval),xrepeat))
    dataplot=plot(normalize.([p_i .* ascertainment for p_i in p.cases],1), markershape = :circle, ylim=ylim, label = p.graphlabel ,seriestype=:scatter,color=color,
    xticks = xticks,xlabel="age",ylabel="proportion",xrotation=ifelse(xrepeat>1,45,0), kwargs...)
    plot!(dataplot,normalize.([p_i .* ascertainment for p_i in p.cases],1), yerror = cases_errors|>permutedims,seriestype=:scatter,ms=0,label=:none,color=color)
end
function Plots.plot(p::Pyramid, multipanel::Bool; color = :black, ascertainment=[1],kwargs...)
    if multipanel
        splitcategories = tuple.(p.casecategories)
        plot.(aggregatecategories.(Ref(p),splitcategories), xticks=(1:length(ageinterval),makeagegrouplabel(ageinterval)); ascertainment=ascertainment, kwargs...)
    else plot(p,kwargs...) end
end
function Plots.plot(base::Plots.Plot, cm::Union{ContactMatrix,AbstractArray{<:ContactMatrix}}; kwargs...)
    ageinterval = first(cm).ageinterval
    linealpha = [ones(length(ageinterval)-1);0]
    plot!(deepcopy(base),dominanteigvec.(cm);linealpha=linealpha, kwargs...)
end

Plots.heatmap(cm::ContactMatrix)=heatmap(collapseblockmat(ngm(cm)),ticks=(1:length(cm.ageinterval),makeagegrouplabel(cm.ageinterval)),xrotation=45)


# output functions

function output_fit(
        pyramid::Pyramid;
        zmb_skeleton::NamedTuple,
        drc_skeleton,
        dataplots,
        ascertainment=1,
        preview = false)
    
    zmb_fit = deepcopy(zmb_skeleton)
    parameters  = [getindex.(Ref(cmt.parameters), [:s_infant,:s_vax]) for cmt in zmb_fit]
    zmb_opt = estimateparameters!(zmb_fit,pyramid,parameters)

    zmb_plot = plot(dataplots, zmb_fit|>collect; label = ["all contacts" "physical only" "at home" "physical & home"],legend=(0.1,0.96),ylim=(0,0.4), color = [1 1 2 2],linestyle = [:solid :dash],
    title = "empirical matrix (Zimbabwe)")

    drc_fit = deepcopy(drc_skeleton)
    parameters  = [getindex.(Ref(cmt.parameters), [:s_infant,:s_vax]) for cmt in drc_fit]
    drc_opt=estimateparameters!(drc_fit,pyramid, parameters)

    drc_plot=plot(dataplots, drc_fit|>collect;label = ["all contacts" "at home"],legend=(0.1,0.96),ylim=(0,0.4),color = [:royalblue :firebrick],
title = "synthetic matrix (DRC)")
    if preview
        zmb_plot|>display
        drc_plot|>display
    end
   (zmb_fit=zmb_fit,zmb_plot=zmb_plot,drc_fit=drc_fit,drc_plot=drc_plot) 
end

function output_fit(
        pyramids::AbstractArray{<:Pyramid};
        zmb_skeleton::AbstractArray,
        drc_skeleton::AbstractArray,
        dataplots,
        ascertainment=1,
        preview = false
        )
    zmb_copy = deepcopy.(zmb_skeleton)
    zmb_fit = vectoriseNamedTuple(zmb_copy)
    parameters  = [getindex.(Ref(cmt_vec[1].parameters), [:s_infant,:s_vax]) for cmt_vec in zmb_fit]
    zmb_opt = estimateparameters!(zmb_fit,pyramids,parameters,sync = [:s_infant,:s_vax,:s_partvax])

    zmb_plot = plot.(dataplots, zmb_copy.|>collect;label = ["all contacts" "physical only" "at home" "physical & home"],legend=(0.1,0.96),ylim=(0,0.4), color = [1 1 2 2],linestyle = [:solid :dash],
    title = "empirical matrix (Zimbabwe)")
    

    drc_copy = deepcopy.(drc_skeleton)
    drc_fit = vectoriseNamedTuple(drc_copy)
    parameters  = [getindex.(Ref(cmt_vec[1].parameters), [:s_infant,:s_vax]) for cmt_vec in drc_fit]
    drc_opt=estimateparameters!(drc_fit,pyramids, parameters,sync=[:s_infant,:s_vax,:s_partvax])

    drc_plot=plot.(dataplots, drc_copy.|>collect; label = ["all contacts" "at home"],legend=(0.1,0.96),ylim=(0,0.4),color = [:royalblue :firebrick],
title = "synthetic matrix (DRC)")
    if preview
        zmb_plot[end]|>display
        drc_plot[end]|>display
    end
   (zmb_fit=zmb_fit,zmb_plot=zmb_plot,drc_fit=drc_fit,drc_plot=drc_plot) 
end


    
function output_sexual_fit(
        pyramid::Pyramid;
        zmb_sexual_skeleton::NamedTuple,
        zmb_ref,
        drc_sexual_skeleton,
        drc_ref,
        dataplots,
        estkeys = [Symbol.(["1_" "2_"],"addmat",4:7)|>vec;:addmat_v;:addmat_w],
        load = nothing,
        ascertainment=1,
        preview = false
    )
    
    if !isnothing(load)
        # load existing data
        zmb_fit,drc_fit  = load.zmb_fit, load.drc_fit
    else
        zmb_fit = deepcopy(zmb_sexual_skeleton)
        overwriteparameters!.(zmb_fit|>collect,zmb_ref|>collect)
        parameters  = [getindex.(Ref(cmt.parameters), estkeys) for cmt in zmb_fit]
        @time zmb_opt = estimateparameters!(zmb_fit,pyramid,parameters;modifier! = propmix!)
        
        drc_fit = deepcopy(drc_sexual_skeleton)
        overwriteparameters!.(drc_fit|>collect,drc_ref|>collect)
        parameters  = [getindex.(Ref(cmt.parameters), estkeys) for cmt in drc_fit]
        @time drc_opt = estimateparameters!(drc_fit,pyramid,parameters;modifier! = propmix!)
    end

    zmb_plot = plot(dataplots, zmb_fit|>collect; ascertainment=ascertainment,label = ["all contacts" "physical only" "at home" "physical & home"],legend=(0.1,0.96),ylim=(0,0.3), color = [1 1 2 2],linestyle = [:solid :dash],
    title = "empirical matrix (Zimbabwe)")

    drc_plot=plot(dataplots, drc_fit|>collect; ascertainment=ascertainment, label = ["all contacts" "at home"],legend=(0.1,0.96),ylim=(0,0.3),color = [:royalblue :firebrick],
title = "synthetic matrix (DRC)")
    
    if preview 
        zmb_plot|>display
        drc_plot|>display
    end
   (zmb_fit=zmb_fit,zmb_plot=zmb_plot,drc_fit=drc_fit,drc_plot=drc_plot) 
end

function output_sexual_fit(
        pyramid::AbstractArray{<:Pyramid};
        zmb_sexual_skeleton::AbstractArray,
        zmb_ref,
        drc_sexual_skeleton::AbstractArray,
        drc_ref,
        dataplots,
        estkeys = [Symbol.(["1_" "2_"],"addmat",4:7)|>vec;:addmat_v;:addmat_w],
        load = nothing,
        ascertainment=1,
        preview = false,
        sync=true
    )

    if !isnothing(load)
        # load existing data
        zmb_fit,drc_fit  = load.zmb_fit, load.drc_fit
    else
        if sync==true sync = estkeys end
        zmb_copy = deepcopy.(zmb_sexual_skeleton)
        zmb_fit = vectoriseNamedTuple(zmb_copy)
        
        for (fit, ref) in zip(zmb_fit,zmb_ref) overwriteparameters!.(fit, ref) end
        
        parameters  = [
            reduce(vcat,[convert(Vector{Scalar}, getindex.(Ref(cmt_vec[x].parameters), x == 1 ? estkeys : setdiff(estkeys,sync))) for x in 1:length(cmt_vec)])
        for cmt_vec in zmb_fit] # include parameters that are in sync only for the first ContactMatrix
        
        @time zmb_opt = estimateparameters!(zmb_fit,pyramid,parameters;modifier! = propmix!,sync=sync)
        
        drc_copy = deepcopy.(drc_sexual_skeleton)
        drc_fit = vectoriseNamedTuple(drc_copy)
        
        for (fit, ref) in zip(drc_fit,drc_ref) overwriteparameters!.(fit, ref) end

        parameters  = [
            reduce(vcat,[convert(Vector{Scalar},getindex.(Ref(cmt_vec[x].parameters), x == 1 ? estkeys : setdiff(estkeys,sync))) for x in 1:length(cmt_vec)])
        for cmt_vec in drc_fit] # include parameters that are in sync only for the first ContactMatrix
        
        @time drc_opt = estimateparameters!(drc_fit,pyramid,parameters;modifier! = propmix!,sync=sync)
    end

    zmb_plot = plot.(dataplots, zmb_copy.|>collect; ascertainment=ascertainment,label = ["all contacts" "physical only" "at home" "physical & home"],legend=(0.1,0.96),ylim=(0,0.3), color = [1 1 2 2],linestyle = [:solid :dash],
    title = "empirical matrix (Zimbabwe)")

    drc_plot=plot.(dataplots, drc_copy.|>collect; ascertainment=ascertainment, label = ["all contacts" "at home"],legend=(0.1,0.96),ylim=(0,0.3),color = [:royalblue :firebrick],
title = "synthetic matrix (DRC)")
    
    if preview 
        zmb_plot.|>display
        drc_plot.|>display
    end
   (zmb_fit=zmb_fit,zmb_plot=zmb_plot,drc_fit=drc_fit,drc_plot=drc_plot) 
end

function output_validate(
        pyramid::Pyramid;
        zmb_skeleton,
        zmb_ref,
        drc_skeleton,
        drc_ref,
        dataplots,
        preview = false
        )

    zmb_fit = deepcopy(zmb_skeleton)
    overwriteparameters!.(zmb_fit|>collect, zmb_ref|>collect)
    @show likelihood.(pyramid,zmb_fit|>collect)
    for cmt in zmb_fit
        if :opt in keys(cmt.misc)
            opt = cmt.misc[:opt]
            cmt.misc[:opt]=(minimizer = opt.minimizer, minimum = -likelihood(pyramid,cmt),  hessian = opt.hessian,result=opt.result)
        end
    end
    zmb_plot = plot(dataplots, zmb_fit|>collect; label = ["all contacts" "physical only" "at home" "physical & home"],legend=(0.1,0.96),ylim=(0,0.4), color = [1 1 2 2],linestyle = [:solid :dash],
    title = "empirical matrix (Zimbabwe)")
    
    drc_fit = deepcopy(drc_skeleton)
    overwriteparameters!.(drc_fit|>collect, drc_ref|>collect)
    @show likelihood.(pyramid,drc_fit|>collect)
    for cmt in drc_fit
        if :opt in keys(cmt.misc)
            opt = cmt.misc[:opt]
            cmt.misc[:opt]=(minimizer = opt.minimizer, minimum = -likelihood(pyramid,cmt),  hessian = opt.hessian,result=opt.result)
        end
    end
    
    drc_plot=plot(dataplots, drc_fit|>collect; label = ["all contacts" "at home"],legend=(0.1,0.96),ylim=(0,0.4),color = [:royalblue :firebrick],
title = "synthetic matrix (DRC)")
    if preview
        zmb_plot|>display
        drc_plot|>display
    end
   (zmb_fit=zmb_fit,zmb_plot=zmb_plot,drc_fit=drc_fit,drc_plot=drc_plot)
end

# summarise results

function parameterestimates(fit)
    parnames = first(fit).parameters|>sort|>keys|>collect
    res = [begin
            estimates = cmt.parameters|>sort|>values|>collect.|>getindex
            min_se = hcat(cmt.misc[:opt].minimizer, cmt.misc[:opt].hessian|>inv|>diag.|>abs.|>sqrt)
            ll = cmt.misc[:opt].minimum
     (estimates,min_se,ll)       
        end for cmt in fit] 
    (parnames=parnames, estimates = first.(res), min_se = getindex.(res,2),ll = last.(res))
end
function CI(fit)
    est = parameterestimates(fit)
    (CI=[[min_se min_se[:,1].+min_se[:,2].*[-1.96 1.96]] for min_se in est.min_se], ll=est.ll)
end

function eigenanalysis(fit)
    res=[begin
            ngmat = ngm(cmt)
            eigvec = groupsplit(eigvecs(collapseblockmat(ngmat))[:,end].|>Real,4)
            eigcases = (ngmat*eigvec)./sum(sum.(ngmat*eigvec))
            eigval0= eigvals(ngmat[2:2:4,2:2:4]|>collapseblockmat)[end]|>Real
           (eigcases, ngmat, eigval0)
        end for cmt in fit]
    (eigcases = first.(res), eigval0=last.(res), eigval = dominanteigval.(fit|>collect),ngm=getindex.(res,2))
end
function vaccineRmap(cmt, R0adjfactor; targetidx = 5:6, R0baseline = 0.82, ve = 0.86, poprange = (0:0.25:50)./100, fswrange = (0:0.5:100)./100)
       addmat_w = cmt.parameters[:addmat_w][]
        addmat5_7 = getindex.(Ref(cmt.parameters),Symbol.("2","_addmat",5:7))
        cmt.susceptibility[1][targetidx].=Ref(fill(1.0))
        Rmap = [begin
                propmix!(Pyramid([],[],Int[],[]), cmt)
                cmt.susceptibility[1][targetidx[1]].=(1-cvg*ve)
                for x in cmt.addmat[:,3] x.*=[ones(3);fill(1-relw*ve,3);ones(2)]' end
                cmt.addmat[4,1].*=[ones(3);fill(1-relw*ve,3);ones(2)]
                dominanteigval(cmt)/R0adjfactor*R0baseline
            end for cvg in poprange, relw in fswrange]
        propmix!(Pyramid([],[],Int[],[]), cmt)
        cmt.susceptibility[1][targetidx].=Ref(cmt.parameters[:s_baseline]) 
        redmap=(1 .-Rmap./(dominanteigval(cmt)/R0adjfactor*R0baseline))   
        (Rmap=Rmap,redmap=redmap)
end

savefigname(filename::AbstractString; save = true) = plt->
begin 
    @suppress if save savefig(plt,filename)|>display end
    @suppress display(plt)
end 
    