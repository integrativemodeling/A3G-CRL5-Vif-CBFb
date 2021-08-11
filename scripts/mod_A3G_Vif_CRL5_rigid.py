############################################
# Modeling the A3G-VCBC complex
# with only three rbs.
#
# iecheverria - Sali Lab - UCSF
############################################
import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.mmcif
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.em
#import IMP.pmi.restraints.residue_binding
from IMP.pmi.io.crosslink import CrossLinkDataBaseKeywordsConverter

import random
import numpy as np
import glob
import sys
import os
from sys import exit
#from sys import argv

###################### SYSTEM SETUP #####################

include_amg_distance_restraint = True
top_dir = '../'

mdl = IMP.Model()

###############################
# Species 1
###############################
top_spec1 = 'top_A3G_CRL5_3rb.dat'
reader_spec1 = IMP.pmi.topology.TopologyReader(top_spec1,
                                               pdb_dir = top_dir+'data/',
                                               fasta_dir = top_dir+'data/')

bs = IMP.pmi.macros.BuildSystem(mdl,
                                resolutions=[1,10])

##############################
# Generate mmcif file
##############################
    
if '--mmcif' in sys.argv:
    import ihm
    # Record the modeling protocol to an mmCIF file
    title = """Characterization of a A3G-Vif HIV-1-CRL5-CBFb structure using a 
             cross-linking mass spectrometry pipeline for integrative modeling 
             of host-pathogen complexes"""
    authors = ['Robyn M Kaake', 'Ignacia Echeverria', 'Seung Joong Kim', 'John Von Dollen',
               'Nicholas M Chesarino', 'Yuqing Feng', 'Clinton Yu', 'Hai Ta', 'Linda Chelico', 
               'Lan Huang', 'John Gross', 'Andrej Sali', 'Nevan J Krogan']
    po = IMP.pmi.mmcif.ProtocolOutput(open('A3G-CRL5-Vif_complex_rigid.cif', 'w'))
    po.system.title = ('Integrative structure determination of the A3G-CRL5-Vif complex')
    bs.system.add_protocol_output(po)
    # Add publication
    po.system.citations.append(ihm.Citation(pmid='000',
                                            title=title,
                                            journal='Molecular & Celullar Proteomics',
                                            authors=authors,
                                            volume='NA',
                                            page_range='NA',
                                            year='2021',
                                            doi='10.1016/j.mcpro.2021.100132'))

##############################
# Build state
##############################
    
bs.add_state(reader_spec1)

hier,  dof = bs.execute_macro(max_rb_trans=3.0,
                                          max_rb_rot=0.03)
mols = bs.get_molecules()[0]


##############################
# Connectivity
##############################
output_objects = [] # keep a list of functions that need to be reported
sample_objects = []
rmf_restraints = []

crs = []
for molname in mols:
    for mol in mols[molname]:
        copy_n = IMP.atom.Copy(mol.get_hierarchy()).get_copy_index()
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.set_label(mol.get_name()+'.'+str(copy_n))
        cr.add_to_model()
        output_objects.append(cr)
        crs.append(cr)


##############################
# Excluded Volume
##############################
evr1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols.values(),
                                                               resolution=10)
evr1.add_to_model()
evr1.set_weight(1.0)
output_objects.append(evr1)

##############################
# Cross-links
##############################
# INITIALIZE DB    

cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Protein1")
cldbkc.set_protein2_key("Protein2")
cldbkc.set_residue1_key("AbsPos1")
cldbkc.set_residue2_key("AbsPos2")
cldbkc.set_unique_id_key("Id")
#cldbkc.set_psi_key("Score")

# XLs RESTRAINT
cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.create_set_from_file(top_dir+"data/Interlinks_A3G_Vif_CRL5_unique_modeling.csv")

xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                                                                            CrossLinkDataBase=cldb,
                                                                            resolution=1.0,
                                                                            length=26.0,
                                                                            slope=0.015)
xl1.add_to_model()
xl1.set_weight(1.5)

rmf_restraints.append(xl1)
output_objects.append(xl1)
dof.get_nuisances_from_restraint(xl1)

#####################################################
# Distance restraints for A3G-Vif
#####################################################
if include_amg_distance_restraint == True:
    br1 = IMP.pmi.restraints.residue_binding.ResidueBindingRestraint(hier,
                                                                    ('A3G',126,132,'Vif'),
                                                                    label = 'A124')
    br1.add_to_model()
    br1.set_weight(8.0)
    output_objects.append(br1)

    br2 = IMP.pmi.restraints.residue_binding.ResidueBindingRestraint(hier,
                                                                    ('Vif',40,45,'A3G'),
                                                                    label = 'V26')
    br2.add_to_model()
    br2.set_weight(8.0)
    output_objects.append(br2)

    
    print('br', br1.get_output(), br2.get_output())
   
##############################
# Shuffle
##############################

IMP.pmi.tools.shuffle_configuration(hier,
                                    max_translation=60)

dof.optimize_flexible_beads(200)
############################# SAMPLING ##############################
# Run replica exchange Monte Carlo sampling

num_frames = 60000
if '--mmcif' in sys.argv or '--test' in sys.argv:
    num_frames=5

rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=hier,                          
                                    crosslink_restraints=rmf_restraints,           
                                    monte_carlo_sample_objects=dof.get_movers(),
                                    replica_exchange_maximum_temperature=3.0,
                                    global_output_directory="output/",
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_frames=num_frames,
                                    number_of_best_scoring_models=0)

rex.execute_macro()

##############################
# Generate mmcif
##############################  
if '--mmcif' in sys.argv:
    import ihm.cross_linkers
    import ihm.dumper
    import ihm.format
    import ihm.location
    import ihm.representation
    import ihm.startmodel
    import ihm.dataset
    import ihm.protocol
    import ihm.analysis
    import ihm.model
    import ihm.restraint
    import ihm.geometry
    import ihm.reference
    import make_archive
    
    fname = '../data/Interlinks_A3G_Vif_CRL5_unique_modeling.csv'

    # Add publication
    po.system.title = "Integrative structure determination of the A3G-CRL5-Vif complex"
    #po.system.citations.append(ihm.Citation.from_pubmed_id(32719457))
    
    s = po.system
    print("restraint datasets:", [r.dataset for r in s.restraints])
    
    # Datasets for XL-MS restraint
    for r in s.restraints:
        if isinstance(r, ihm.restraint.CrossLinkRestraint):
            r.linker = ihm.cross_linkers.dsso
            print("XL-MS dataset at:", r.dataset.location.path)
            print("Details:", r.dataset.location.details)
    
    # Correct number of output models to account for multiple runs
    protocol = po.system.orphan_protocols[-1]
    protocol.steps[-1].num_models_end = 937225

    # Get last protocol in the file
    protocol = po.system.orphan_protocols[-1]
    # State that we filtered 
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.ClusterStep(
                            feature='RMSD', num_models_begin=3000000,
                            num_models_end=937225))
    
    
    # Create an ensemble for the cluster
    e = po._add_simple_ensemble(analysis.steps[-1],
                                name="Cluster 0", num_models=857561,
                                drmsd=9.27, num_models_deposited=1,
                                localization_densities={}, ensemble_file=None)
    
    # Add the model from RMF
    rh = RMF.open_rmf_file_read_only(
        '../results/model_rigid/clustering_rigid/cluster.0/cluster_center_model.rmf3')
    IMP.rmf.link_hierarchies(rh, [hier])
    IMP.rmf.load_frame(rh, RMF.FrameID(0))
    del rh
    model = po.add_model(e.model_group)
    
    # Add localization densities
    # Look up the ihm.AsymUnit corresponding to a PMI component name
    for asym in po.asym_units:
        name = asym.split('.')[0]
        fname = f'../results/model_rigid/clustering_rigid/cluster.0/LPD_{name}_orie.mrc'
        print('fname', fname)
        loc = ihm.location.OutputFileLocation(fname)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=po.asym_units[asym])
        # Add to ensemble
        e.densities.append(den)

    
    # Add uniprot of proteins
    lpep = ihm.LPeptideAlphabet()
    d = 'Construct used for crystallization'
    sd_vif = [ihm.reference.SeqDif(33, lpep['R'], lpep['K'], details=d),
              ihm.reference.SeqDif(39, lpep['V'], lpep['F'], details=d),
              ihm.reference.SeqDif(41, lpep['K'], lpep['R'], details=d),
              ihm.reference.SeqDif(48, lpep['N'], lpep['H'], details=d),
              ihm.reference.SeqDif(50, lpep['R'], lpep['K'], details=d),
              ihm.reference.SeqDif(51, lpep['I'], lpep['V'], details=d),
              ihm.reference.SeqDif(61, lpep['E'], lpep['D'], details=d),
              ihm.reference.SeqDif(64, lpep['I'], lpep['L'], details=d),
              ihm.reference.SeqDif(65, lpep['I'], lpep['V'], details=d),
              ihm.reference.SeqDif(67, lpep['R'], lpep['T'], details=d),
              ihm.reference.SeqDif(74, lpep['I'], lpep['T'], details=d),
              ihm.reference.SeqDif(77, lpep['K'], lpep['R'], details=d),
              ihm.reference.SeqDif(92, lpep['N'], lpep['R'], details=d),
              ihm.reference.SeqDif(95, lpep['H'], lpep['S'], details=d),
              ihm.reference.SeqDif(105, lpep['H'], lpep['Q'], details=d),
              ihm.reference.SeqDif(110, lpep['Y'], lpep['H'], details=d),
              ihm.reference.SeqDif(125, lpep['I'], lpep['L'], details=d),
              ihm.reference.SeqDif(127, lpep['E'], lpep['Q'], details=d),
              ihm.reference.SeqDif(151, lpep['K'], lpep['T'], details=d),
              ihm.reference.SeqDif(153, lpep['V'], lpep['L'], details=d),
              ihm.reference.SeqDif(154, lpep['V'], lpep['I'], details=d),
              ihm.reference.SeqDif(155, lpep['A'], lpep['T'], details=d),
              ihm.reference.SeqDif(156, lpep['S'], lpep['P'], details=d),
              ihm.reference.SeqDif(157, lpep['T'], lpep['K'], details=d),
              ihm.reference.SeqDif(158, lpep['R'], lpep['K'], details=d),
              ihm.reference.SeqDif(170, lpep['V'], lpep['T'], details=d),
              ihm.reference.SeqDif(173, lpep['R'], lpep['H'], details=d),
              ihm.reference.SeqDif(174, lpep['W'], lpep['H'], details=d)]
    
    # name : (uniprot id, mutations, [[db_begin, db_end, entity_begin, entity_end]]
    Uniprot={'A3G.0': ('Q9HC16',[],[]),
             'Vif.0':  ('Q90QQ9',sd_vif,[[1,174,1,174]]),
             'CBFB.0': ('Q13951',[],[]),
             'EloB.0': ('Q15370',[],[]),
             'EloC.0': ('Q15369',[],[]),
             'CUL5.0': ('Q93034',[],[]),
             'Rbx2.0': ('Q9UBF6',[],[])}
    
    for prot, (entry, sd, limits) in Uniprot.items():
        print(prot, entry, sd, limits)
        ref = ihm.reference.UniProtSequence.from_accession(entry)
        for seg in limits:
            ref.alignments.append(ihm.reference.Alignment(
                db_begin=seg[0], db_end=seg[1], entity_begin=seg[2], entity_end=seg[3], seq_dif=sd))
            
        po.asym_units[prot].entity.references.append(ref)
    
    
    # Point to the raw mass spec data and peaklists used to derive the crosslinks.
    l = ihm.location.PRIDELocation('PXD025391',
                                   details='All raw mass spectrometry files and '
                                   'peaklists used in the study')
    xl1.dataset.add_primary(ihm.dataset.MassSpecDataset(location=l))
        
    # Replace local links with DOIs
    repos = []
    for subdir, zipname in make_archive.ARCHIVES.items():
        print('subdir', subdir)
        repos.append(ihm.location.Repository(
            doi="10.5281/zenodo.5176959", root="../%s" % subdir,
            url="https://zenodo.org/record/10.5281/files/%s.zip" % zipname,
            top_directory=os.path.basename(subdir)))
    
    #po.system.update_locations_in_repositories(repos)
    
    
    po.flush()




