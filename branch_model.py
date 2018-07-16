from ete3 import EvolTree
from Bio import SeqIO
import re
import sys
from ete3 import NCBITaxa


# first argument, where your trees are
# second argument, where your alignments are
# third argument, where your control file is
# fourth argument, where your null control file is


# read each line in tree
ncbi = NCBITaxa()
record_dict = SeqIO.index(sys.argv[2], 'fasta')  #### second argument ###
run = 0
out = open('output.txt','w')
out.write('tree #\ttree likelihood\tS\D\tMRCA\tAge\tbl_1\tdN_1\tdS_1\tw_1\tleaves_1\tbl_2\tdN_2\tdS_2\tw_2\tleaves_2\tw_diff\n')
with open(sys.argv[1]) as f: #### first argument ####
    for line in f:
        run = run + 1
        tree = EvolTree(line, format= 1)
        tree.workdir = '.'
        internal_node_names=[]
        for node in tree.traverse():
            # split up species name and protein name, based on conventions I had set up earlier
            if node.is_leaf():
                tiplist = node.name.rsplit('__', 1)
                node.name = tiplist[0]
                node.add_features(foot=tiplist[-1])
            else:
                internal_node_names.append(node.name)
        if all(node == '' for node in internal_node_names):
            print "no calibration points in tree, skipping"
            continue

        for node in tree.traverse():
            # switch internal node names with calibration times for PAML
            if node.is_leaf() == False:
                if node.name == 'AV':
                    node.name = '@4.248'
                if node.name == 'J':
                    node.name = '@3.225'
                if node.name == 'I':
                    node.name = '@2.506'
                if node.name == 'H':
                    node.name = '@1.729'
                if node.name == 'C':
                    node.name = '@1.268'
                if node.name == 'B':
                    node.name = '@0.965'
                if node.name == 'F':
                    node.name = '@1.048'
                if node.name == 'E':
                    node.name = '@0.356'
                if node.name == 'BY':
                    node.name = '@4.094'
                if node.name == 'BX':
                    node.name = '@3.354'
                if node.name == 'BV':
                    node.name = '@2.998'
                if node.name == 'BU':
                    node.name = '@2.645'
                if node.name == 'BT':
                    node.name = '@0.719'
                if node.name == 'BS':
                    node.name = '@0.237'
                if node.name == 'BR':
                    node.name = '@0.546'
                if node.name == 'BQ':
                    node.name = '@0.15'
                if node.name == 'BP':
                    node.name = '@1.847'
                if node.name == 'BO':
                    node.name = '@1.728'
                if node.name == 'BN':
                    node.name = '@0.76'
                if node.name == 'BL':
                    node.name = '@0.892'
                if node.name == 'L':
                    node.name = '@0.603'
                if node.name == 'K':
                    node.name = '@0.694'
                if node.name == 'N':
                    node.name = '@0.871'
                if node.name == 'BK':
                    node.name = '@0.839'
                if node.name == 'Y':
                    node.name = '@0.76'
                if node.name == 'O':
                    node.name = '@0.615'
                if node.name == 'AA':
                    node.name = '@0.727'
                if node.name == 'T':
                    node.name = '@0.717'
                if node.name == 'X':
                    node.name = '@0.614'
                if node.name == 'W':
                    node.name = '@0.581'
                if node.name == 'V':
                    node.name = '@0.528'
                if node.name == 'Z':
                    node.name = '@0.699'
                if node.name == 'R':
                    node.name = '@0.542'
                if node.name == 'Q':
                    node.name = '@0.485'
                if node.name == 'P':
                    node.name = '@0.455'
                if node.name == 'BJ':
                    node.name = '@0.757'
                if node.name == 'AU':
                    node.name = '@0.708'
                if node.name == 'AB':
                    node.name = '@0.479'
                if node.name == 'AS':
                    node.name = '@0.644'
                if node.name == 'AR':
                    node.name = '@0.323'
                if node.name == 'AL':
                    node.name = '@0.465'
                if node.name == 'AK':
                    node.name = '@0.359'
                if node.name == 'AJ':
                    node.name = '@0.27'
                if node.name == 'AF':
                    node.name = '@0.14'
                if node.name == 'AH':
                    node.name = '@0.205'
                if node.name == 'BH':
                    node.name = '@0.682'
                if node.name == 'AW':
                    node.name = '@0.551'
                if node.name == 'BG':
                    node.name = '@0.65'
                if node.name == 'BF':
                    node.name = '@0.374'
                if node.name == 'BE':
                    node.name = '@0.264'
                if node.name == 'BD':
                    node.name = '@0.131'
                if node.name == 'BC':
                    node.name = '@0.105'
                if node.name == 'BA':
                    node.name = '@0.171'
                if node.name == 'AZ':
                    node.name = '@0.104'
                if node.name == 'AY':
                    node.name = '@0.082'

        # species names are now in feature 'foot'
        def get_species_name(node):
            clean = node.foot
            return clean
        # get speciation and duplication events based on species of leaves, make empty alignment
        tree.set_species_naming_function(get_species_name)
        tree.get_descendant_evol_events()
        align = ""

        # get MRCA for internal nodes, alignment from leaves
        for node in tree.traverse():
            if node.is_leaf() == False:
                leaf_names = []
                leaves = node.get_leaves()
                for leaf in leaves:
                    name = ncbi.get_name_translator([leaf.foot])
                    name_val = name.values()
                    name_val1 = name_val[0]
                    name_val2 = name_val1[0]
                    leaf_names.append(name_val2)
                nodetree = ncbi.get_topology(leaf_names)
                kids = nodetree.get_leaves()
                MRCA = nodetree.get_common_ancestor(kids)
                clade = ncbi.get_taxid_translator([MRCA.name])
                node.add_features(common_ancestor=clade.values()[0])

            else:
                align = align + str(record_dict[node.name].format("fasta"))

        align = align.replace('/','')
        tree.link_to_alignment(align)

        # run null model (M0) to get K and BL priors
        M0_ctrl = open(sys.argv[4]).read() #### null string ###
        tree.run_model('XX.null' + '_' + "{0:02d}".format(run), ctrl_string=M0_ctrl)
        null_model = tree.get_evol_model('XX.null' + '_' + "{0:02d}".format(run))
        null_tree = EvolTree(null_model._tree.write(format=1),format=1)
        null_tree.workdir = '.'
        null_tree.link_to_alignment(align)


        # branch model
        string = open(sys.argv[3]).read() #### control string ###
        string = re.sub('( kappa = ).+', '\g<1>' + str(null_model.stats['kappa']), string)
        null_tree.run_model('XX.' + "{0:02d}".format(run), ctrl_string=string)
        current_model = null_tree.get_evol_model('XX.' + "{0:02d}".format(run))


        for node in tree.traverse():
            out = open('output.txt', 'a')
            if node.is_leaf() == False:
                number = str(run)
                likelihood = str(current_model.lnL)
                EvolEvent = node.evoltype
                CA = node.common_ancestor
                Age = str(node.get_distance(node.get_leaves()[0]))
                br_length1 = str(node.get_distance(node.get_children()[0]))
                dn1 = str(current_model.branches[node.get_children()[0].node_id]['dN'])
                ds1 = str(current_model.branches[node.get_children()[0].node_id]['dS'])
                w1 = str(current_model.branches[node.get_children()[0].node_id]['w'])
                leaves = []
                for leaf in (node.get_children()[0]).get_leaves():
                    leaves.append(leaf.name)
                    leaves1 = ", ".join(leaves)
                br_length2 = str(node.get_distance(node.get_children()[1]))
                dn2 = str(current_model.branches[node.get_children()[1].node_id]['dN'])
                ds2 = str(current_model.branches[node.get_children()[1].node_id]['dS'])
                w2 = str(current_model.branches[node.get_children()[1].node_id]['w'])
                leaves = []
                for leaf in (node.get_children()[1]).get_leaves():
                        leaves.append(leaf.name)
                leaves2 = ", ".join(leaves)
                contrast = str(abs(float(w1)-float(w2)))
                out.write(number + '\t' + likelihood + '\t' + EvolEvent + '\t' + CA + '\t' + Age + '\t' + br_length1 + '\t' + dn1 + '\t' + ds1 + '\t' + w1 + '\t' + leaves1 + '\t' + br_length2 + '\t' + dn2 + '\t' + ds2 + '\t' + w2 + '\t' + leaves2 + '\t' + contrast + '\n')
            out.close()
