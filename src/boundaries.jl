abstract type AbstractBoundaryCondition end
struct ZeroPadded <: AbstractBoundaryCondition end
struct Periodic <: AbstractBoundaryCondition end
struct PeriodicExtended <: AbstractBoundaryCondition
    boundary::Int
end

function get_raster_size(raster, boundary::AbstractBoundaryCondition)
    size(raster)
end
function get_raster_size(raster, boundary::PeriodicExtended)
    (size(raster)[1:end-1]..., size(raster)[end] - 2*boundary.boundary)
end