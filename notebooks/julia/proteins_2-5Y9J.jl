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

    for i in 1:n_vectors
        if i % 100 == 0
            println(string("Calculating scores (rotation ", i, " out of ", n_vectors, ')'))
        end
        v = rotation_vectors[i,:]
        rotation = recenter(transpose(get_rotation(v)), center(ligand_volume))
        ligand_rotated = warp(ligand_volume, rotation, axes(ligand_volume), method=BSpline(Constant()), fillvalue=0)

        ift_rc = n_voxels * ifft(ligand_rotated)
        inv_val = ifft(ift_rc .* dft_lc)
        inv_val = fftshift(inv_val)
        ssc = real(inv_val) .- imag(inv_val)

        cutoff = percentile(vec(ssc), 99.9)
        top_indexes = findall(>(cutoff), ssc)

        # this is probably slowing this function down significantly...
        top_scores = [top_scores; [(i, idx, ssc[idx]) for idx in top_indexes]]
    end

    top_scores
end

# -----------------------------------------LOAD DATA---------------------------------------------------------------

println("Loading data...")

proteins_dict = Dict()
protein_complex_pdb_id = ARGS[1]
figures_path = mkpath(string(protein_complex_pdb_id, "_figures"))

for designation in ("l_u", "r_u", "l_b", "r_b")
    protein = read(string("./protein_data/", protein_complex_pdb_id, '_', designation, ".pdb"), PDBFormat)
    length(protein.models) == 1 || throw(DomainError(protein.models), "pdb files should only contain one model")
    proteins_dict[designation] = protein[1]
end

# make sure the bound and unbound proteins have the same chain IDs
keys(proteins_dict["l_u"].chains) == keys(proteins_dict["l_b"].chains) || throw(DomainError(keys(proteins_dict["l_u"].chains), "bound and unbound l proteins have different chain IDs"))
keys(proteins_dict["r_u"].chains) == keys(proteins_dict["r_b"].chains) || throw(DomainError(keys(proteins_dict["r_u"].chains), "bound and unbound r proteins have different chain IDs"))

# all atoms in the unbound complex
ligand_atoms, receptor_atoms = collectatoms(proteins_dict["l_u"]), collectatoms(proteins_dict["r_u"])
ligand_coords, receptor_coords = [atom.coords for atom in ligand_atoms], [atom.coords for atom in receptor_atoms]

# spatial midpoints of unbound proteins
ligand_midpt, receptor_midpt = calc_midpoint(ligand_coords), calc_midpoint(receptor_coords)
ligand_coords, receptor_coords = [c - ligand_midpt for c in ligand_coords], [c - receptor_midpt for c in receptor_coords]

# gather all amino acids into a dictionary where the keys are (chain ID, amino acid ID) so the same amino acids can be compared
# between the bound and unbound versions of the proteins
all_residues_dict = Dict()

for k in keys(proteins_dict)
    protein = proteins_dict[k]
    residues_dict = Dict()
    for chain_id in keys(protein.chains)
        for res in collectresidues(protein[chain_id])
            # ignore non-protein molecules
            res.het_res && continue
            coords = [atom.coords - receptor_midpt for atom in values(res.atoms)]
            residues_dict[(chain_id, res.number)] = Dict("coords" => coords, "res_num" => res.number, "is_interface" => false, "c_a_coords" => res.atoms[" CA "].coords - receptor_midpt )
        end
    end
    all_residues_dict[k] = residues_dict
end


println("----------------------------------------------------------------------------")

# -----------------------------------------EXTRACT INTERFACE, PLOT DATA AND CONFIRM PAPER IRMSD----------------------------------

println("Extracting interface and calculating I-RMSD...")

# determine if each amino acid is part of the interface of the complex
interface_max_dist = 10.0

for res_key_l in keys(all_residues_dict["l_b"])
    atom_l_coords = all_residues_dict["l_b"][res_key_l]["coords"]

    for res_key_r in keys(all_residues_dict["r_b"])
        atom_r_coords = all_residues_dict["r_b"][res_key_r]["coords"]
        is_interface = any([norm(coord2 - coord1) <= interface_max_dist for coord2 in atom_r_coords for coord1 in atom_l_coords])
        
        if is_interface
            all_residues_dict["r_b"][res_key_r]["is_interface"] = true
            haskey(all_residues_dict["r_u"], res_key_r) && (all_residues_dict["r_u"][res_key_r]["is_interface"] = true)

            all_residues_dict["l_b"][res_key_l]["is_interface"] = true
            haskey(all_residues_dict["l_u"], res_key_l) && (all_residues_dict["l_u"][res_key_l]["is_interface"] = true)
        end
    end
end

# plot interaface atoms of bound and unbound proteins together and plot the bound complex
interface_l, interface_r, complex = get_plot_data(all_residues_dict)

plot_3d(interface_l, string("ligand_compare_", protein_complex_pdb_id), figures_path)
plot_3d(interface_r, string("receptor_compare_", protein_complex_pdb_id), figures_path)
plot_3d(complex, string("complex_", protein_complex_pdb_id), figures_path)

# calculate interface rmsd
irmsd = calc_interface_rmsd(all_residues_dict)
println(string("The I-RMSD for ", protein_complex_pdb_id, " is ", irmsd))

println("----------------------------------------------------------------------------")

# -----------------------------------------CREATE 3D VOLUMES-------------------------------------------------------

println("Discretizing proteins into 3D volumes...")

# rotate unbound ligand protein by arbitrary matrix
init_rotation_vec = [0.646, 0.676, 0.355]
init_rotation_vec /= norm(init_rotation_vec)
inv_rotation_vec = [init_rotation_vec[1], -1.0*init_rotation_vec[2], -1.0*init_rotation_vec[3]]
init_rotation = get_rotation(init_rotation_vec)

println(string("Ligand rotation vector: ", round.(inv_rotation_vec, sigdigits=3)))

all_residues_dict["l_u"] = apply_transformation(all_residues_dict["l_u"], init_rotation, receptor_midpt-ligand_midpt)

ligand_coords = [[dot(coord, init_rotation[1,:]), dot(coord, init_rotation[2,:]), dot(coord, init_rotation[3,:])] for coord in ligand_coords]

# Van der Waal radius of each atom type
atom_radius_dict = Dict("C"=> 1.700, "N"=> 1.550, "O"=> 1.520, "S"=> 1.800, "NA"=> 2.270, "MG"=> 1.730)

# https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1992.pdf page 28, section B ii
ligand_radius_list = [atom_radius_dict[strip(atom.name[1:2])] for atom in ligand_atoms]
receptor_radius_list = [atom_radius_dict[strip(atom.name[1:2])] for atom in receptor_atoms]

# distance matrix for surface area calculation
ligand_dist_matrix, receptor_dist_matrix = calc_dist_matrix(ligand_atoms), calc_dist_matrix(receptor_atoms)

# calculate accessible surface area
n_points = 1000
sphere_points = get_sphere_points(n_points)

# calculate surface area
ligand_atoms_asa = get_asa(ligand_coords, ligand_radius_list, ligand_dist_matrix, sphere_points)
receptor_atoms_asa = get_asa(receptor_coords, receptor_radius_list, receptor_dist_matrix, sphere_points)

ligand_max_dist, receptor_max_dist = maximum(ligand_dist_matrix), maximum(receptor_dist_matrix)

volume_size = ligand_max_dist + receptor_max_dist + 2.0
volume_dim = 128
center_val = Int(volume_dim // 2)
voxel_size = volume_size / (volume_dim - 1)
vol_start = -0.5 * volume_size

println(string("voxel size: ", voxel_size))

receptor_volume, ligand_volume = zeros(volume_dim, volume_dim, volume_dim), zeros(volume_dim, volume_dim, volume_dim)
vols = Dict('l' => ligand_volume, 'r' => receptor_volume)

asa_min = 1.0

# receptor -- iteration 1: [core: sqrt(1.5)*r, surface: sqrt(0.8)*r], iteration 2: [core: 0.0 (skip), surface: r + 3.4]
rr_iter1 = [receptor_atoms_asa[i] < asa_min ? sqrt(1.5)*receptor_radius_list[i] : sqrt(0.8)*receptor_radius_list[i] for i in 1:length(receptor_atoms_asa)]
rr_iter2 = [receptor_atoms_asa[i] < asa_min ? 0.0 : receptor_radius_list[i]+3.4 for i in 1:length(receptor_atoms_asa)]

# ligand -- iteration 1: [core: sqrt(1.5)*r, surface: 0.0 (skip)], iteration 2: [core: 0.0 (skip), surface: 1.0 * r]
lr_iter1 = [ligand_atoms_asa[i] < asa_min ? sqrt(1.5)*ligand_radius_list[i] : 0.0 for i in 1:length(ligand_atoms_asa)]
lr_iter2 = [ligand_atoms_asa[i] < asa_min ? 0.0 : ligand_radius_list[i] for i in 1:length(ligand_atoms_asa)]

for (iter, designation, radius_list) in [(1, 'r', rr_iter1), (2, 'r', rr_iter2), (1, 'l', lr_iter1), (2, 'l', lr_iter2)]
    atoms = if designation == 'r' receptor_coords else ligand_coords end
    set_val = if iter == 1 2.0 else 1.0 end

    for atom_idx in 1:length(atoms)
        atom, atom_r = atoms[atom_idx], radius_list[atom_idx]

        vol_idxs = [1 + floor( ( atom[i] - vol_start ) / voxel_size ) for i in 1:3]
        n_neighbor_voxels = ceil(atom_r / voxel_size)
        start_end_idx = [(max(1, vol_idxs[i] - n_neighbor_voxels), min(vol_idxs[i] + n_neighbor_voxels, volume_dim)) for i in 1:3]

        for x in start_end_idx[1][1]:start_end_idx[1][2]
            for y in start_end_idx[2][1]:start_end_idx[2][2]
                for z in start_end_idx[3][1]:start_end_idx[3][2]
                    ( (iter == 2) && (designation == 'r') && (vols[designation][Int(x),Int(y),Int(z)] == 2.0) ) && continue

                    voxel_coords = [x, y, z]
                    voxel_center = [vol_start + voxel_size * (voxel_coords[i] - 0.5) for i in 1:3]
                    center_dist = norm(voxel_center - atoms[atom_idx])

                    (center_dist <= atom_r) && (vols[designation][Int(x),Int(y),Int(z)] = set_val)
                end
            end
        end
    end
end

change_idxs = []

for x in 2:volume_dim-1
    for y in 2:volume_dim-1
        for z in 2:volume_dim-1
            (ligand_volume[x,y,z] != 2.0) && continue
            # ambiguous in ZDOCK paper, but using 1-connected and 2-connected neighbors (neighbors that share a corner or a face)
            neighbor_vals = [ligand_volume[x1,y1,z1] for z1 in z-1:z+1 for y1 in y-1:y+1 for x1 in x-1:x+1 if (x1==x || y1==y || z1==z)]
            (sum([val == 0 for val in neighbor_vals]) >= 2) && (push!(change_idxs, (x,y,z)))
        end
    end
end

for (x, y, z) in change_idxs
    ligand_volume[x,y,z] = 1.0
end

# plot volume slices
layout = Layout(width=800, height=800, scene_aspectratio=attr(x=1, y=1))
p1 = plot(heatmap(z=vols['l'][:,:,center_val]), layout)
p2 = plot(heatmap(z=vols['r'][:,:,center_val]), layout)
plt = [p1 p2]
relayout!(plt)

open(string(string(figures_path, "/vol_slices_", protein_complex_pdb_id, ".html")), "w") do io
    PlotlyBase.to_html(io, plt.plot)
end

println("----------------------------------------------------------------------------")

# -----------------------------------------SCAN 6D SPACE USING FFT---------------------------------------------------------

println("Scanning 6D search space...")

# convert volumes to complex numbers
rho = 9im
ligand_volume =  convert(Array{ComplexF64, 3}, ifelse.(ligand_volume .> 1, rho, ligand_volume))
receptor_volume = convert(Array{ComplexF64, 3}, ifelse.(receptor_volume .> 1, rho, receptor_volume))

n_rotations = 500
retain_n_scores = 2000

rotation_vectors = get_sphere_points(n_rotations)
top_scores = scan_6d(rotation_vectors, receptor_volume, ligand_volume)
# retain top n scores
top_scores = (sort(top_scores, by=last, rev=true))[1:retain_n_scores]

println("----------------------------------------------------------------------------")

# -----------------------------------------EVALUATE CANDIDATE TRANSFORMATIONS---------------------------------------------------------

println("Evaluating candidate transformations...")

cutoff_a = 2.5

l_u_copy = deepcopy(all_residues_dict["l_u"])
close_irmsds = []

for score_rank in 1:length(top_scores)
    rotation_idx, t, score = top_scores[score_rank]
    translation = voxel_size * ([t[1],t[2],t[3]] - [center_val+1,center_val+1,center_val+1])
    all_residues_dict["l_u"] = apply_transformation(l_u_copy, get_rotation(rotation_vectors[rotation_idx,:]), translation)
    i_rmsd = calc_interface_rmsd(all_residues_dict)

    (i_rmsd < cutoff_a) && push!(close_irmsds, (score_rank, rotation_vectors[rotation_idx,:], score, i_rmsd))
end


close_irmsds = sort(close_irmsds, by=last, rev=true)

println(string(length(close_irmsds), " out of ", retain_n_scores, " with I-RMSD < ", cutoff_a))

for candidate_idx in 1:length(close_irmsds)
    rank, rotation_vector, score, i_rmsd = close_irmsds[candidate_idx]

    println(string("\n-----Candidate ", candidate_idx, "-----"))
    println(string("Score rank: ", rank))
    println(string("Score: ", score))
    println(string("I-RMSD: ", i_rmsd))
    println(string("Rotation vector: ", round.(rotation_vector, sigdigits=3)))

    theta = acos(dot(rotation_vector, inv_rotation_vec)) * 180 / pi
    println(string("Angle formed with correct inverse vector (degrees): ", theta))
end
