# >>>>>>>>> copied code from AbTrimerFibril_TREM2dock_tryingnewmols.py <<<<<<<<<<<<<
import IMP
import IMP.atom
import IMP.isd
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.io
import IMP.pmi.io.crosslink


def create_fibril_transforms(z_displacement, len_protofilament):

    transforms = []

    for t in range(1, len_protofilament):
        if t % 2 == 0:  # even
            trans_pos = IMP.algebra.Transformation3D(IMP.algebra.Vector3D(0, 0, ((t / 2) * z_displacement)))
            transforms.append(trans_pos)
            print("Transforms:", t, " | ", trans_pos)
        else: # odd
            trans_neg = IMP.algebra.Transformation3D(IMP.algebra.Vector3D(0, 0, ((-1 * t) / 2) * z_displacement))
            transforms.append(trans_neg)
            print("Transforms:", t, " | ", trans_neg)
    return transforms


# -----------------------------CREATE SYSTEM AND STATE ------------------------------------
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
bs = IMP.pmi.macros.BuildSystem(mdl)
st = s.create_state()

# ##Input Data####
pdb_dir = "../data/pdb/"
fasta_dir = "../data/fasta/"
xl_data_BS3 = "../data/xl/BS3_inter.csv"    # Experimental cross-links set-1 24 A long.
xl_data_NHSF = "../data/xl/NHSF_inter.csv"  # Experimental cross-links set-2 20 A long.
xl_weight = 1.0

# Store the FASTA sequences in a dictionary
sequences = IMP.pmi.topology.Sequences(fasta_dir + "TREM2_ABfibril.fasta")

# _____ ACTIVE Abeta Trimer_______________________________________________________________
# Create a molecule from 1 trimer (chain A) of the amyloid fibril (2LMP)
# & copy it to build the ACTIVE trimer

Abeta = st.create_molecule("Abeta", sequence=sequences["Abeta40"], chain_id='A')
AbetaG = Abeta.create_copy('G')
AbetaM = Abeta.create_copy('M')
AbA = Abeta.add_structure(pdb_dir + "Ab_tri_AGM.pdb", chain_id='A', offset=0)
AbG = AbetaG.add_structure(pdb_dir + "Ab_tri_AGM.pdb", chain_id='G', offset=0)
AbM = AbetaM.add_structure(pdb_dir + "Ab_tri_AGM.pdb", chain_id='M', offset=0)

# Add representation to the monomers in the ACTIVE Abeta trimer
# Chain A monomer
Abeta.add_representation(AbA, resolutions=[1], color=0.1)
Abeta.add_representation(Abeta.get_non_atomic_residues(), resolutions=[1], color=0.1)
# Chain G monomer
AbetaG.add_representation(AbG, resolutions=[1], color=0.1)
AbetaG.add_representation(AbetaG.get_non_atomic_residues(), resolutions=[1], color=0.1)
# Chain M monomer
AbetaM.add_representation(AbM, resolutions=[1], color=0.1)
AbetaM.add_representation(AbetaM.get_non_atomic_residues(), resolutions=[1], color=0.1)

# ______Fibril extension____

# Create the transforms needed to build up the fibril from symmetry z-displacement translation

fibril_transforms = create_fibril_transforms(z_displacement=4.8, len_protofilament=10)
print('Transforms  || ', len(fibril_transforms), " || ", fibril_transforms)

# Create a clone of the trimer fibril unit for every transform
# Collect all trimer fibril units created in a list

fibril_mols = [Abeta, AbetaG, AbetaM]
protofil_mols = [[Abeta],[AbetaG],[AbetaM]]

chains = ['A', 'G', 'M']
additional_chains = ['B', 'C', 'D', 'E', 'F', 'H', 'I', 'J', 'K', 'L', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z','B', 'C', 'D', 'E', 'F', 'H', 'I', 'J', 'K', 'L', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
count =0
for i in range(0, len(chains)):
    for nc in range(len(fibril_transforms)):
        clone = fibril_mols[i].create_clone(additional_chains[count])  # chains[nc])
        protofil_mols[i].append(clone)
        fibril_mols.append(clone)
        count += 1

# ______TREM2 ECD__________________________________________________________________________
# Create a molecule for 1 TREM2 Extra Cellular Domain (based on 5UD7, ch. A)
trem2 = st.create_molecule('TREM2', sequence=sequences["TREM2"])

T = trem2.add_structure(pdb_dir + "TREM2_A.pdb",
                        chain_id='A',
                        offset=0)  # offset -18
trem2.add_representation(trem2.get_atomic_residues(), resolutions=[1], color=0.5)

trem2_clones = [trem2]
# clone two more molecules of TREM2 such that the stiochiometry is maintained
for i in range(0, 2):
    clone = trem2.create_clone(additional_chains[count])
    count +=1
    trem2_clones.append(clone)

# _____ Build the system & degrees of freedom _______________________________________________
root_hier = s.build()
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)

# ----------- DEFINE THE SYSTEM RIGID BODIES & Set up connectivity restraints -----------------------------
# Lists for collecting molecules/particles
TREM2_mol = []
active_tri_mol = []
fib_ext_mol = []
shuffle_exclude_rbs = []

# Lists useful for collecting restraints
output_objects = []
sample_objects = []
crs = []

# get a dictionary of the molecules in the system w/ KEY = molecule name
moldict = st.get_molecules()
print(moldict)

for molname in moldict:
    for mol in moldict[molname]:
        if 'TREM2' in molname:
            atomic = mol.get_atomic_residues()
            dof.create_rigid_body(atomic,
                                  max_trans=1.0,
                                  max_rot=0.5,
                                  resolution='all')
            TREM2_mol.append(mol)
            # print("******************** TREM2 	", dof.get_rigid_bodies())

        elif 'Abeta' in molname:
            crA = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol, scale=2.0)
            crA.add_to_model()
            output_objects.append(crA)
            sample_objects.append(crA)
            crs.append(crA)

            atomic = mol.get_atomic_residues()
            dof.create_rigid_body(atomic,
                                  max_trans=0,
                                  max_rot=0,
                                  resolution='all')
            dof.create_flexible_beads(mol.get_non_atomic_residues(),
                                      max_trans=1.0,
                                      resolution=1)
            active_tri_mol.append(mol)
            #print(" ******************* Abeta 	", dof.get_rigid_bodies())

# ______Composite Restraint for Fibril___________________________________________________________________

# print("TREM2_mol 		||	", TREM2_mol)
# print("active_tri_mol 		||	", active_tri_mol)
# print('********** just TREM2 :', dof.get_rigid_bodies()[:-1])
shuffle_exclude_rbs = dof.get_rigid_bodies()[:-1]

# Constrain the fibril trimer unit copies with z-trans. symmetry ### This is very important
for i in range(0, len(chains)):
    for t in range(0, len(fibril_transforms)):
        threefold_trans = fibril_transforms[t]
        dof.constrain_symmetry(protofil_mols[i][0], protofil_mols[i][t + 1], threefold_trans)
mdl.update()

print("|||||  All sym constraints added  |||||")

# Write a single-frame RMF to view the system
out = IMP.pmi.output.Output()
out.init_rmf("symmetry_test1.rmf3", hierarchies=[root_hier])
out.write_rmf("symmetry_test1.rmf3")

# -----------------------------ADD RESTRAINTS ------------------------------------

# ______External Barrier Restraint__________________________________________________________________________
# get the center of mass of the fibril
Fibril_sel = IMP.atom.Selection(root_hier, molecule="Abeta")
fib_particles = Fibril_sel.get_selected_particles()
if len(fib_particles) == 0:
    print("COM not set up. Cannot select protein %s)" & (Abeta))
print(fib_particles)

fibril_COM = IMP.atom.CenterOfMass.setup_particle(IMP.Particle(mdl), fib_particles)
fibril_COM_coor = IMP.core.XYZ(fibril_COM).get_coordinates()

# set up an an external barrier sphere with radius equal to the length of the fibril from the center of mass to
# exclude docking from fibril ends
eb = IMP.pmi.restraints.basic.ExternalBarrier(hierarchies=root_hier, radius=50, center=fibril_COM_coor, label='barrier')
eb.add_to_model()
output_objects.append(eb)

# ______Excluded Volume Restraint__________________________________________________________________________

ev1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=TREM2_mol, resolution=1)
ev1.add_to_model()
ev1.set_label('TREM2')
output_objects.append(ev1)

ev2 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=TREM2_mol,
                                                              other_objects=[fibril_mols],
                                                              resolution=5)
ev2.add_to_model()
ev2.rs.set_weight(10.0)
ev2.set_label('all')
output_objects.append(ev2)

# ______Cross-linking Restraint__________________________________________________________________________
xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()

# Set up XL restraint for BS3 XLs
crosslink_restraints = []
xldb_BS3 = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_BS3.create_set_from_file(file_name=xl_data_BS3, converter=xldbkc)

xlrB = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=root_hier,  # Must pass the root hierarchy to the system
    CrossLinkDataBase=xldb_BS3,  # The crosslink database.
    length=24,  # The crosslinker plus side chain length
    resolution=1,  # The resolution at which to evaluate the crosslink
    slope=0.02,
    label='BS3',  # This adds a linear term to the scoring function to bias crosslinks towards each other
    weight=xl_weight)  # Scaling factor for the restraint score.
xlrB.add_to_model()
output_objects.append(xlrB)
crosslink_restraints.append(xlrB)

xldb_NHSF = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_NHSF.create_set_from_file(file_name=xl_data_NHSF, converter=xldbkc)
xlrN = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=root_hier,  # Must pass the root hierarchy to the system
    CrossLinkDataBase=xldb_NHSF,  # The cross-link database.
    length=20,  # The cross-linker plus side chain length
    resolution=1,  # The resolution at which to evaluate the cross-link
    slope=0.02,
    label='NHSF',  # This adds a linear term to the scoring function to bias cross-links towards each other
    weight=xl_weight)  # Scaling factor for the restraint score.
xlrN.add_to_model()
output_objects.append(xlrN)
crosslink_restraints.append(xlrN)

# ----------------------------- SAMPLING ------------------------------------
for n in range(10): # run 10 independent runs over #no of replicas suggested by mpirun
    outdir = "testsystem3/run"
    run_num = n + 1
    output = '/output'
    global_output_directory = outdir + str(run_num) + output
    print(global_output_directory)

    IMP.pmi.tools.shuffle_configuration(root_hier, excluded_rigid_bodies=shuffle_exclude_rbs,
                                        max_translation=150)  # bounding_box = bb)
    dof.optimize_flexible_beads(500)
    out = IMP.pmi.output.Output()

    rex = IMP.pmi.macros.ReplicaExchange0(mdl,
                                          root_hier=root_hier,
                                          monte_carlo_sample_objects=dof.get_movers(),
                                          output_objects=output_objects,
                                          crosslink_restraints=[crosslink_restraints],
                                          monte_carlo_temperature=1.0,
                                          number_of_best_scoring_models=0,
                                          monte_carlo_steps=10,
                                          number_of_frames=50000,
                                          global_output_directory=global_output_directory)
    rex.execute_macro()
    n += 1

out = IMP.pmi.output.Output()
out.init_rmf("singleframe.rmf3", hierarchies=[root_hier])
out.write_rmf("singleframe.rmf3")
