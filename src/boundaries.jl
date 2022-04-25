abstract type AbstractBoundaryCondition end
struct ZeroPadded <: AbstractBoundaryCondition end
struct Periodic <: AbstractBoundaryCondition end
struct PeriodicExtended <: AbstractBoundaryCondition
    boundary::Int
end