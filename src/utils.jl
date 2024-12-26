function groupsum(v::AbstractVector, breaks::AbstractVector)
    ids = range.(breaks[begin:end-1],breaks[begin+1:end].-1)
    sum.(getindex.(Ref(v),[ids;[breaks[end]:lastindex(v)]]))
end
function groupmean(v::AbstractVector, breaks::AbstractVector)
    ids = range.(breaks[begin:end-1],breaks[begin+1:end].-1)
    mean.(getindex.(Ref(v),[ids;[breaks[end]:lastindex(v)]]))
end
function groupsplit(v::AbstractVector, npieces::Integer)
    if !iszero(length(v)%npieces) throw("length is not multiple of the number of pieces to split into") end
    len = length(v)÷npieces
    getindex.(Ref(v),range.(1:len:length(v),len:len:length(v)))
end
function splitsum(v::AbstractVector, npieces::Integer, breaks::AbstractVector)
    groupsum(groupsplit(v,npieces),breaks)
end
    
function vectoriseNamedTuple(a::NamedTuple, b::NamedTuple)
    if(keys(a)!=keys(b)) throw("keys do not match") end
    NamedTuple((key,[a[Symbol(key)],b[Symbol(key)]]) for key in keys(a))
end
function vectoriseNamedTuple(v::AbstractArray{<:NamedTuple})
    if !(keys.(v).|>Set|>allequal) throw("keys do not match") end
    NamedTuple(((key,getfield.(v,key)) for key in keys(first(v))))
end


collapseblockmat(bm)= collapseblockmat(reduce(vcat, (reduce(hcat, row) for row in eachrow(bm))))
collapseblockmat(bm::AbstractArray{<:Scalar})= getindex.(bm)
collapseblockmat(bm::AbstractArray{<:Number})= bm
