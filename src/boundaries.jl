abstract type AbstractBoundaryCondition end
struct ZeroPadded <: AbstractBoundaryCondition end
struct Periodic <: AbstractBoundaryCondition end
struct PeriodicExtended <: AbstractBoundaryCondition
    t_bounds::Tuple{Int,Int}
end