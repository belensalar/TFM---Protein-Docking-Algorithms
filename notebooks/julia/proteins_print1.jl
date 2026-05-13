using BioStructures
using LinearAlgebra
using PlotlyJS
using DataFrames
using StaticArrays
using Rotations
using ImageTransformations
using CoordinateTransformations
using Interpolations
using FFTW
using StatsBase

#= Julia implementation of the algorithm described in "Docking Unbound Proteins Using Shape Complementarity,
Desolvation, and Electrostatics", Rong Chen and Zhiping Weng (2002). Implements the shape complementarity scoring
function only. Data is from https://github.com/piercelab/antibody_benchmark.
=#

function get_plot_data(residues_dict::Dict)
    interface_l, interface_r, complex = [], [], []

    for designation in ('l', 'r')
        bound_residues = residues_dict[string(designation, "_b")]
        unbound_residues = residues_dict[string(designation, "_u")]

        for bound in (true, false)
            residues = if bound bound_residues else unbound_residues end
            for res in values(residues)
                x, y, z = res["c_a_coords"][1], res["c_a_coords"][2], res["c_a_coords"][3]

                if bound
                    d = if designation == 'l' "ligand" else "receptor" end
                    res["is_interface"] && (d = string(d, "_interface"))
                    
                    complex_plot_dict = Dict("x_coord" => x, "y_coord" => y, "z_coord" => z, "color_val" => d)
                    push!(complex, complex_plot_dict)
                end

                res["is_interface"] || continue
                plot_dict = Dict("x_coord" => x, "y_coord" => y, "z_coord" => z, "color_val" => res["res_num"] % 5)                
                if designation == 'l' push!(interface_l, plot_dict) else push!(interface_r, plot_dict) end
            end
        end
    end
    (interface_l, interface_r, complex)
end

function plot_3d(plot_dict::Vector, output_fname::String, figures_path::String)
    df = DataFrame(plot_dict)
    plt = plot(df, x=:x_coord, y=:y_coord, z=:z_coord, color=:color_val, type="scatter3d", mode="markers")
    open(string(figures_path, '/', output_fname, ".html"), "w") do io
        PlotlyBase.to_html(io, plt.plot)
    end
end

function apply_transformation(residues_dict::Dict, rotation::Matrix{Float64}, translation::Vector)#, midpt::Vector = [0, 0, 0])
    copy = Dict()
    midpt = calc_midpoint([val["c_a_coords"] for val in values(residues_dict)])
    for k in keys(residues_dict)
        val = residues_dict[k]
        coords = val["c_a_coords"]
        new_coords = coords - midpt
        new_coords = [dot(new_coords, rotation[1,:]), dot(new_coords, rotation[2,:]), dot(new_coords, rotation[3,:])]
        new_coords += translation + midpt
        copy[k] = Dict("is_interface" => val["is_interface"], "c_a_coords" => new_coords, "res_num" => val["res_num"])
    end
    copy
end

function get_rotation(B::Vector)
    A = [1, 0, 0]
    G = [dot(A,B) -norm(cross(A,B)) 0; norm(cross(A,B)) dot(A,B) 0; 0 0 1]
    Fi = [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ]
    rotation_matrix = Fi*G*inv(Fi)
    rotation_matrix
end

function calc_midpoint(coords::Vector)
    m = reduce(vcat,transpose.(coords))
    0.5 * [minimum(m[:,1]) + maximum(m[:,1]), minimum(m[:,2]) + maximum(m[:,2]), minimum(m[:,3]) + maximum(m[:,3])]
end


function calc_interface_rmsd(proteins_dict::Dict)
    sum, n = 0, 0
    
    for designation in ('l', 'r')
        bound_residues = proteins_dict[string(designation, "_b")]
        unbound_residues = proteins_dict[string(designation, "_u")]

        for res_key in keys(bound_residues)
            res_bound = bound_residues[res_key]
            (haskey(unbound_residues, res_key) && res_bound["is_interface"]) || continue

            res_unbound = unbound_residues[res_key]
            dist2 = norm(res_bound["c_a_coords"] - res_unbound["c_a_coords"]) ^ 2
            sum += dist2
            n += 1
        end
    end
    rmsd = sqrt(sum / n)
end

function calc_dist_matrix(atoms::Vector)
    out = zeros(length(atoms), length(atoms))

    for k in 1:length(atoms)
        @inbounds out[k,k] = 0.0
        for j in 1:(k-1) 
        @inbounds out[j,k] = norm(atoms[j].coords - atoms[k].coords)
        end
    end
    Symmetric(out)
end

function get_sphere_points(n::Int)
    dl = pi * (3.0 - sqrt(5.0))
    dz = 2.0 / n

    longitude = 0
    z = 1 - dz / 2

    coords = zeros(n, 3)
    for k in 1:n
        r = sqrt((1 - z * z))
        coords[k, 1] = cos(longitude) * r
        coords[k, 2] = sin(longitude) * r
        coords[k, 3] = z
        z -= dz
        longitude += dl
    end

    coords
end

# based on https://biopython.org/docs/dev/api/Bio.PDB.SASA.html#Bio.PDB.SASA.ShrakeRupley
function get_asa(atoms::Vector, radius_list::Vector, dist_matrix::Symmetric{Float64, Matrix{Float64}}, sphere_points = nothing)
    if sphere_points == nothing
        sphere_points = get_sphere_points(100)
    end

    probe_radius = 1.4

    max_radius = 2.0 * ( maximum(radius_list) + probe_radius )
    neighbors = [ findall(<=(max_radius + 0.01), dist_matrix[i,:]) for i in 1:(size(dist_matrix)[1])]

    atom_asa = zeros(length(atoms))

    for atom_idx in 1:length(atoms)
        atom_radius = radius_list[atom_idx]
        neighbor_indexes = neighbors[atom_idx]

        atom_sphere_points = [atoms[atom_idx] + ( (atom_radius + probe_radius) * sphere_points[i,:]) for i in 1:(size(sphere_points)[1])]
        
        points_inaccessible = [false for n in 1:n_points]

        for neighbor_idx in neighbor_indexes
            (atom_idx == neighbor_idx) && continue

            neighbor_atom = atoms[neighbor_idx]
            neighbor_radius = radius_list[neighbor_idx]

            neighbor_contains_points = [ (norm(neighbor_atom - pt) < (neighbor_radius + probe_radius)) for pt in atom_sphere_points]
            points_inaccessible = [ (neighbor_contains_points[pt_idx] || points_inaccessible[pt_idx]) for pt_idx in 1:n_points]
        end
        
        percent_accessible = 1.0 - ( sum(points_inaccessible) / n_points )

        surface_area = 4 * pi * ( (atom_radius + probe_radius) ^2 )
        atom_asa[atom_idx] = percent_accessible * surface_area
    end
    atom_asa
end

function scan_6d(rotation_vectors::Matrix{Float64}, receptor_volume::Array{ComplexF64, 3}, ligand_volume::Array{ComplexF64, 3})
    top_scores = []
    vol_size = size(receptor_volume)
    n_voxels = vol_size[1] * vol_size[2] * vol_size[3]
    dft_lc = fft(receptor_volume)
    n_vectors = size(rotation_vectors)[1]

