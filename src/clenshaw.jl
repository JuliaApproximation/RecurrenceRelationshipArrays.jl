###
# Operator clenshaw
###

Base.@propagate_inbounds function clenshaw_next!(n, A::AbstractFill, ::Zeros, C::Ones, x::AbstractMatrix, c, bn1::AbstractMatrix{T}, bn2::AbstractMatrix{T}) where T
    muladd!(getindex_value(A), x, bn1, -one(T), bn2)
    view(bn2,band(0)) .+= c[n]
    bn2
end

Base.@propagate_inbounds function clenshaw_next!(n, A::AbstractVector, ::Zeros, C::AbstractVector, x::AbstractMatrix, c, bn1::AbstractMatrix{T}, bn2::AbstractMatrix{T}) where T
    muladd!(A[n], x, bn1, -C[n+1], bn2)
    view(bn2,band(0)) .+= c[n]
    bn2
end


# allow special casing first arg, for ChebyshevT in ClassicalOrthogonalPolynomials
Base.@propagate_inbounds function _clenshaw_first!(A, ::Zeros, C, X, c, bn1, bn2) 
    muladd!(A[1], X, bn1, -C[2], bn2)
    view(bn2,band(0)) .+= c[1]
    bn2
end


Base.@propagate_inbounds function _clenshaw_first!(A, ::Zeros, C, X, c, f::AbstractVector, bn1, bn2) 
    muladd!(A[1], X, bn1, -C[2], bn2)
    bn2 .+= c[1] .* f
    bn2
end

Base.@propagate_inbounds function clenshaw_next!(n, A::AbstractVector, B::AbstractVector, C::AbstractVector, x::AbstractMatrix, c, bn1::AbstractMatrix{T}, bn2::AbstractMatrix{T}) where T
    # bn2 .= B[n] .* bn1 .- C[n+1] .* bn2
    lmul!(-C[n+1], bn2)
    LinearAlgebra.axpy!(B[n], bn1, bn2)
    muladd!(A[n], x, bn1, one(T), bn2)
    view(bn2,band(0)) .+= c[n]
    bn2
end

# Operator * f Clenshaw

Base.@propagate_inbounds function clenshaw_next!(n, A, B, C, X::AbstractMatrix, c, f::AbstractVector, bn1::AbstractVector{T}, bn2::AbstractVector{T}) where T
    bn2 .= B[n] .* bn1 .- C[n+1] .* bn2 .+ c[n] .* f
    muladd!(A[n], X, bn1, one(T), bn2)
    bn2
end

Base.@propagate_inbounds function _clenshaw_first!(A, B, C, X, c, bn1, bn2) 
    lmul!(-C[2], bn2)
    LinearAlgebra.axpy!(B[1], bn1, bn2)
    muladd!(A[1], X, bn1, one(eltype(bn2)), bn2)
    view(bn2,band(0)) .+= c[1]
    bn2
end

Base.@propagate_inbounds function _clenshaw_first!(A, B, C, X, c, f::AbstractVector, bn1, bn2) 
    bn2 .= B[1] .* bn1 .- C[2] .* bn2 .+ c[1] .* f
    muladd!(A[1], X, bn1, one(eltype(bn2)), bn2)
    bn2
end


# FillArrays
Base.@propagate_inbounds function clenshaw_next!(n, A::AbstractFill, ::Zeros, C::Ones, X::AbstractMatrix, c, f::AbstractVector, bn1::AbstractVector{T}, bn2::AbstractVector{T}) where T
    muladd!(getindex_value(A), X, bn1, -one(T), bn2)
    bn2 .+= c[n] .* f
    bn2
end

Base.@propagate_inbounds function clenshaw_next!(n, A, ::Zeros, C, X::AbstractMatrix, c, f::AbstractVector, bn1::AbstractVector{T}, bn2::AbstractVector{T}) where T
    muladd!(A[n], X, bn1, -C[n+1], bn2)
    bn2 .+= c[n] .* f
    bn2
end


_clenshaw_op(::AbstractBandedLayout, Z, N) = BandedMatrix(Z, (N-1,N-1))

function clenshaw(c::AbstractVector, A::AbstractVector, B::AbstractVector, C::AbstractVector, X::AbstractMatrix)
    N = length(c)
    T = promote_type(eltype(c),eltype(A),eltype(B),eltype(C),eltype(X))
    @boundscheck check_clenshaw_recurrences(N, A, B, C)
    m = size(X,1)
    m == size(X,2) || throw(DimensionMismatch("X must be square"))
    N == 0 && return zero(T)
    bn2 = _clenshaw_op(MemoryLayout(X), Zeros{T}(m, m), N)
    bn1 = _clenshaw_op(MemoryLayout(X), c[N]*Eye{T}(m), N)
    _clenshaw_op!(c, A, B, C, X, bn1, bn2)
end

function clenshaw(c::AbstractVector, A::AbstractVector, B::AbstractVector, C::AbstractVector, X::AbstractMatrix, f::AbstractVector)
    N = length(c)
    T = promote_type(eltype(c),eltype(A),eltype(B),eltype(C),eltype(X))
    @boundscheck check_clenshaw_recurrences(N, A, B, C)
    m = size(X,1)
    m == size(X,2) || throw(DimensionMismatch("X must be square"))
    m == length(f) || throw(DimensionMismatch("Dimensions must match"))
    N == 0 && return [zero(T)]
    bn2 = zeros(T,m)
    bn1 = Vector{T}(undef,m)
    bn1 .= c[N] .* f
    _clenshaw_op!(c, A, B, C, X, f, bn1, bn2)
end

function _clenshaw_op!(c, A, B, C, X, bn1, bn2)
    N = length(c)
    N == 1 && return bn1
    @inbounds begin
        for n = N-1:-1:2
            bn1,bn2 = clenshaw_next!(n, A, B, C, X, c, bn1, bn2),bn1
        end
        bn1 = _clenshaw_first!(A, B, C, X, c, bn1, bn2)
    end
    bn1
end

function _clenshaw_op!(c, A, B, C, X, f::AbstractVector, bn1, bn2)
    N = length(c)
    N == 1 && return bn1
    @inbounds begin
        for n = N-1:-1:2
            bn1,bn2 = clenshaw_next!(n, A, B, C, X, c, f, bn1, bn2),bn1
        end
        bn1 = _clenshaw_first!(A, B, C, X, c, f, bn1, bn2)
    end
    bn1
end


"""
    Clenshaw(c, A, B, C, X, p0=1)

represents the operator `a(X)` where a is a polynomial
where `a` is represented by coefficients `c` in a basis dictated by the recurrence coefficients
`A`, `B`, and `C`. `p0` is the initial condition in the recurrence relationship.
"""
struct Clenshaw{T, Coefs<:AbstractVector, AA<:AbstractVector, BB<:AbstractVector, CC<:AbstractVector, Jac<:AbstractMatrix} <: AbstractBandedMatrix{T}
    c::Coefs
    A::AA
    B::BB
    C::CC
    X::Jac
    p0::T
end

Clenshaw(c::AbstractVector{T}, A::AbstractVector, B::AbstractVector, C::AbstractVector, X::AbstractMatrix{T}, p0=one(T)) where T = 
    Clenshaw{T,typeof(c),typeof(A),typeof(B),typeof(C),typeof(X)}(c, A, B, C, X, p0)

function Clenshaw(c::AbstractVector, A::AbstractVector, B::AbstractVector, C::AbstractVector, X::AbstractMatrix, p0...)
    T = promote_type(eltype(c), eltype(X))
    Clenshaw(convert(AbstractVector{T}, c), A, B, C, convert(AbstractMatrix{T},X), p0...)
end

Clenshaw(c::Number, A, B, C, X, p) = Clenshaw([c], A, B, C, X, p)

copy(M::Clenshaw) = M
size(M::Clenshaw) = size(M.X)
axes(M::Clenshaw) = axes(M.X)
bandwidths(M::Clenshaw) = (length(M.c)-1,length(M.c)-1)

Base.array_summary(io::IO, C::Clenshaw{T}, inds::Tuple{Vararg{OneToInf{Int}}}) where T =
    print(io, Base.dims2string(length.(inds)), " Clenshaw{$T} with $(length(C.c)) degree polynomial")

struct ClenshawLayout <: AbstractLazyBandedLayout end
MemoryLayout(::Type{<:Clenshaw}) = ClenshawLayout()
sublayout(::ClenshawLayout, ::Type{<:NTuple{2,AbstractUnitRange{Int}}}) = ClenshawLayout()
sublayout(::ClenshawLayout, ::Type{<:Tuple{AbstractUnitRange{Int},Union{Slice,AbstractInfUnitRange{Int}}}}) = LazyBandedLayout()
sublayout(::ClenshawLayout, ::Type{<:Tuple{Union{Slice,AbstractInfUnitRange{Int}},AbstractUnitRange{Int}}}) = LazyBandedLayout()
sublayout(::ClenshawLayout, ::Type{<:Tuple{Union{Slice,AbstractInfUnitRange{Int}},Union{Slice,AbstractInfUnitRange{Int}}}}) = LazyBandedLayout()
sub_materialize(::ClenshawLayout, V) = BandedMatrix(V)

function _BandedMatrix(::ClenshawLayout, V::SubArray{<:Any,2})
    M = parent(V)
    kr,jr = parentindices(V)
    b = bandwidth(M,1)
    jkr = max(1,min(first(jr),first(kr))-b÷2):max(last(jr),last(kr))+b÷2
    # relationship between jkr and kr, jr
    kr2,jr2 = kr.-first(jkr).+1,jr.-first(jkr).+1
    lmul!(M.p0, clenshaw(M.c, M.A, M.B, M.C, layout_getindex(M.X,jkr, jkr))[kr2,jr2])
end

function getindex(M::Clenshaw{T}, kr::AbstractUnitRange, j::Integer) where T
    b = bandwidth(M,1)
    jkr = max(1,min(j,first(kr))-b÷2):max(j,last(kr))+b÷2
    # relationship between jkr and kr, jr
    kr2,j2 = kr.-first(jkr).+1,j-first(jkr)+1
    f = [Zeros{T}(j2-1); one(T); Zeros{T}(length(jkr)-j2)]
    lmul!(M.p0, clenshaw(M.c, M.A, M.B, M.C, M.X[jkr, jkr], f)[kr2])
end

getindex(M::Clenshaw, k::Int, j::Int) = M[k:k,j][1]

function getindex(S::Symmetric{T,<:Clenshaw}, k::Integer, jr::AbstractUnitRange) where T
    m = max(jr.start,jr.stop,k)
    return Symmetric(getindex(S.data,1:m,1:m),Symbol(S.uplo))[k,jr]
end
function getindex(S::Symmetric{T,<:Clenshaw}, kr::AbstractUnitRange, j::Integer) where T
    m = max(kr.start,kr.stop,j)
    return Symmetric(getindex(S.data,1:m,1:m),Symbol(S.uplo))[kr,j]
end
function getindex(S::Symmetric{T,<:Clenshaw}, kr::AbstractUnitRange, jr::AbstractUnitRange) where T
    m = max(kr.start,jr.start,kr.stop,jr.stop)
    return Symmetric(getindex(S.data,1:m,1:m),Symbol(S.uplo))[kr,jr]
end

transposelayout(M::ClenshawLayout) = LazyBandedLayout()
# TODO: generalise for layout, use Base.PermutedDimsArray
Base.permutedims(M::Clenshaw{<:Number}) = transpose(M)


function materialize!(M::MatMulVecAdd{<:ClenshawLayout,<:AbstractPaddedLayout,<:AbstractPaddedLayout})
    α,A,x,β,y = M.α,M.A,M.B,M.β,M.C
    length(y) == size(A,1) || throw(DimensionMismatch("Dimensions must match"))
    length(x) == size(A,2) || throw(DimensionMismatch("Dimensions must match"))
    x̃ = paddeddata(x);
    m = length(x̃)
    b = bandwidth(A,1)
    jkr=1:m+b
    p = [x̃; zeros(eltype(x̃),length(jkr)-m)];
    Ax = lmul!(A.p0, clenshaw(A.c, A.A, A.B, A.C, A.X[jkr, jkr], p))
    _fill_lmul!(β,y)
    resizedata!(y, last(jkr))
    v = view(paddeddata(y),jkr)
    LinearAlgebra.axpy!(α, Ax, v)
    y
end




##
# Banded dot is slow
###

LinearAlgebra.dot(x::AbstractVector, A::Clenshaw, y::AbstractVector) = dot(x, mul(A, y))

