"""
Calculate matching score script.

This script demonstrates how to use inverse_msmd package to calculate
matching scores based on interaction profiles.
"""
from inverse_msmd.utils.bio_utils import SuperImposer, PDB
from gridData import Grid
import numpy as np

# GAMMA=0.003
# pairs = [("5C84", "EZZ")]
GAMMA=0.00

def calculate_matching_score(atoms_of_interest, profile_residues, probe_center, profiles):
    """
    Calculate matching score based on interaction profiles.
    
    Parameters
    ----------
    atoms_of_interest : list
        List of atoms to consider
    profile_residues : list
        List of residue names with profiles
    probe_center : np.ndarray
        Center coordinates of the probe
    profiles : dict
        Dictionary of Grid objects for each residue
        
    Returns
    -------
    float
        Log matching score
    """
    log_score = 0
    for atom in atoms_of_interest:
        resname = atom.get_parent().get_resname()
        if resname not in profile_residues:
            continue # ほんとうはこんなものあってはいけない、、、でもないか。GLYはどうしてもここにくる
        # if resname == "ASP":
        #     continue # E16について、spline補完のせいで負の値が出てしまったので一旦除外
        coord = atom.get_coord()
        distance = np.linalg.norm(probe_center - coord)
        weight = np.exp(-GAMMA * (distance**2))
        log_score += np.log(max(profiles[resname].interpolated([coord[0]], [coord[1]], [coord[2]])[0], profiles[resname].grid.min())) * weight
    return log_score


# Main execution
if __name__ == "__main__":
    matching_ids = "A08 A08_1 A08_2 A17 A57 A58 A58_1 B23 B39 E21 E21_1 E21_2 E22 E23 E24 E24_1 E24_2 E25 E25_1 E26 E27 E28 E28_1 E41 E42 E43 E44 E45 E46 E47 E48 E49 E50 E51 E52".split(" ")
    for matching in matching_ids:
        probe_id = matching[:3]
        protein = PDB.get_structure(f"4hw3_aligned_to_{matching}.pdb")
        probe   = PDB.get_structure(f"/home/7/16D30270/workspace/9998_share/msmd/probe/{probe_id}/{probe_id}.pdb")
        probe_center = PDB.get_attr(probe, "coord").mean(axis=0)
        # print(probe_center)
        atoms_of_interest = [a for a in protein.get_atoms() if a.get_name() == "CB"]

        profile_residues = "ALA CYS GLU ASP PHE HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split(" ")
        profiles = {res: Grid(f"../profile_generator/interaction_profile/{probe_id}_{res}_profile.dx") for res in profile_residues}
        # print(profiles)
        score = calculate_matching_score(atoms_of_interest, profile_residues, probe_center, profiles)
        print(matching, f"{score:.2f}")