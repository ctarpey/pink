os.chdir("Y:\WORK\WAPLES\Stacks_mapping\Python")

import cPickle as pickle
import numpy

import Stacks_SQL_class
from models_to_test import models_to_test



# Load data created with export_sql.pl 

mykiss_fams = Stacks_SQL_class.Stacks_SQL(path = "Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_fams_RW_v2.txt", 
    name = 'chum',
    num_parents = 14, 
    num_offspring = 497, 
    parent_offspring = (77, 75, 78, 19, 40, 20, 20, 20, 18, 18, 30, 33, 43, 6)
    )
    
# Add TaqMan Data
#chum_fams.add_external_genotypes(filepath = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/genotypes/taqMan_01H.tsv", parent = 'CMUW10X_0001')
#chum_fams.add_external_genotypes(filepath = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/genotypes/taqMan_08H.tsv", parent = 'CMUW10X_0008')
#chum_fams.add_external_genotypes(filepath = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/genotypes/taqMan_09H.tsv", parent = 'CMUW10X_0009')

# Kick out Loci (within families) with too many missing genotypes
too_many_missing_01 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_1')
too_many_missing_02 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_2')
too_many_missing_03 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_3')
too_many_missing_04 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_4')
too_many_missing_05 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_5')
too_many_missing_06 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_6')
too_many_missing_07 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_7')
too_many_missing_08 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_8')
too_many_missing_09 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_9')
too_many_missing_10 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_10')
too_many_missing_11 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_11')
too_many_missing_12 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_12')
too_many_missing_13 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_13')
too_many_missing_14 = mykiss_fams.remove_missing_within_cross(max_miss = .25, parent = 'parent_14')


# Kick out loci (within families) with too many alleles
#too_many_alleles_01 = chum_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'CMUW10X_0001')
too_many_alleles_01 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_1')
too_many_alleles_02 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_2')
too_many_alleles_03 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_3')
too_many_alleles_04 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_4')
too_many_alleles_05 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_5')
too_many_alleles_06 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_6')
too_many_alleles_07 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_7')
too_many_alleles_08 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_8')
too_many_alleles_09 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_9')
too_many_alleles_10 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_10')
too_many_alleles_11 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_11')
too_many_alleles_12 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_12')
too_many_alleles_13 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_13')
too_many_alleles_14 = mykiss_fams.remove_max_alleles_within_cross(max_alleles = 4, parent = 'parent_14')

# Purge CatIDs that failed previous tests in all families
purged = mykiss_fams.purge_absent_catIDs()

# Calculate model results using models found in in "models_to_test.py"
print("Calculating PSV Model Likelihoods")
mykiss_fams.set_model_results(to_test = models_to_test, epsilon_bound_high = 0.1, epsilon_bound_low = 0.01)


## STOPPED HERE 4/18/2014
## We may want to exclude certain individuals and re-run the above analyses

#mappable_01, stats_01 = chum_fams.find_mappable_markers(parent = 'CMUW10X_0001', models_to_test = models_to_test)
mappable_01, stats_01 = mykiss_fams.find_mappable_markers(parent = 'parent_1', models_to_test = models_to_test)
mappable_02, stats_02 = mykiss_fams.find_mappable_markers(parent = 'parent_2', models_to_test = models_to_test)
mappable_03, stats_03 = mykiss_fams.find_mappable_markers(parent = 'parent_3', models_to_test = models_to_test)
mappable_04, stats_04 = mykiss_fams.find_mappable_markers(parent = 'parent_4', models_to_test = models_to_test)
mappable_05, stats_05 = mykiss_fams.find_mappable_markers(parent = 'parent_5', models_to_test = models_to_test)
mappable_06, stats_06 = mykiss_fams.find_mappable_markers(parent = 'parent_6', models_to_test = models_to_test)
mappable_07, stats_07 = mykiss_fams.find_mappable_markers(parent = 'parent_7', models_to_test = models_to_test)
mappable_08, stats_08 = mykiss_fams.find_mappable_markers(parent = 'parent_8', models_to_test = models_to_test)
mappable_09, stats_09 = mykiss_fams.find_mappable_markers(parent = 'parent_9', models_to_test = models_to_test)
mappable_10, stats_10 = mykiss_fams.find_mappable_markers(parent = 'parent_10', models_to_test = models_to_test)
mappable_11, stats_11 = mykiss_fams.find_mappable_markers(parent = 'parent_11', models_to_test = models_to_test)
mappable_12, stats_12 = mykiss_fams.find_mappable_markers(parent = 'parent_12', models_to_test = models_to_test)
mappable_13, stats_13 = mykiss_fams.find_mappable_markers(parent = 'parent_13', models_to_test = models_to_test)
mappable_14, stats_14 = mykiss_fams.find_mappable_markers(parent = 'parent_14', models_to_test = models_to_test)

def write_stats(filename, stats):
    with open(filename, 'w') as OUTFILE:
        for catID, values in stats.items():
            OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                catID, 
                values.epsilon,
                values.best_model,
                values.best_likelihood,
                values.second_best_model,
                values.second_best_likelihood,
                values.x1_seg_pval,
                values.x2_seg_pval
                )
            )
#from psv_working import stats_08
#write_stats('/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/psv/chum_01.stats', stats_01)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_01.stats', stats_01)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_02.stats', stats_02)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_03.stats', stats_03)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_04.stats', stats_04)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_05.stats', stats_05)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_06.stats', stats_06)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_07.stats', stats_07)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_08.stats', stats_08)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_09.stats', stats_09)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_10.stats', stats_10)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_11.stats', stats_11)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_12.stats', stats_12)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_13.stats', stats_13)
write_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_14.stats', stats_14)

# Conduct filtering on loci and catIDs
# Kick out failed loci
# segregation testing
import mne
import collections
def filter_mappable_segregation(stats, alpha = .05):
    stats_of_locus = dict()
    loci = list()
    raw_p_vals = list()
    for catID in stats:
        if stats[catID].x1_seg_pval == "NA":
            pass
        else:
            raw_p_vals.append(stats[catID].x1_seg_pval)
            loci.append(catID + "_x1")
        if stats[catID].x2_seg_pval == "NA":
            pass
        else:
            raw_p_vals.append(stats[catID].x2_seg_pval)
            loci.append(catID + "_x2")
    reject, pval_corrected = mne.stats.fdr_correction(pvals = raw_p_vals, alpha = alpha, method = 'indep')
    seg_test_fdr = collections.namedtuple('seg_test_fdr', ['reject', 'pval_corrected', 'pval_raw'])
    for locus, rej, pvc, pvr in zip(loci, reject, pval_corrected, raw_p_vals):
        stats_of_locus[locus] = seg_test_fdr(rej, pvc, pvr)
    return(stats_of_locus)
# do segregation testing with FDR    


seg_test_01 = filter_mappable_segregation(stats_01, alpha = 0.05)
seg_test_02 = filter_mappable_segregation(stats_02, alpha = 0.05)
seg_test_03 = filter_mappable_segregation(stats_03, alpha = 0.05)
seg_test_04 = filter_mappable_segregation(stats_04, alpha = 0.05)
seg_test_05 = filter_mappable_segregation(stats_05, alpha = 0.05)
seg_test_06 = filter_mappable_segregation(stats_06, alpha = 0.05)
seg_test_07 = filter_mappable_segregation(stats_07, alpha = 0.05)
seg_test_08 = filter_mappable_segregation(stats_08, alpha = 0.05)
seg_test_09 = filter_mappable_segregation(stats_09, alpha = 0.05)
seg_test_10 = filter_mappable_segregation(stats_10, alpha = 0.05)
seg_test_11 = filter_mappable_segregation(stats_11, alpha = 0.05)
seg_test_12 = filter_mappable_segregation(stats_12, alpha = 0.05)
seg_test_13 = filter_mappable_segregation(stats_13, alpha = 0.05)
seg_test_14 = filter_mappable_segregation(stats_14, alpha = 0.05)


# Print how many *would be* excluded
print("Reject {} out of {} loci in fam_01".format(sum([x.reject for x in seg_test_01.values()]), len([x.reject for x in seg_test_01.values()])))
print("Reject {} out of {} loci in fam_08".format(sum([x.reject for x in seg_test_08.values()]), len([x.reject for x in seg_test_08.values()])))
print("Reject {} out of {} loci in fam_09".format(sum([x.reject for x in seg_test_09.values()]), len([x.reject for x in seg_test_09.values()])))

failed_seg_test_01 = [k for k,v in seg_test_01.items() if v.reject]
failed_seg_test_02 = [k for k,v in seg_test_02.items() if v.reject]
failed_seg_test_03 = [k for k,v in seg_test_03.items() if v.reject]
failed_seg_test_04 = [k for k,v in seg_test_04.items() if v.reject]
failed_seg_test_05 = [k for k,v in seg_test_05.items() if v.reject]
failed_seg_test_06 = [k for k,v in seg_test_06.items() if v.reject]
failed_seg_test_07 = [k for k,v in seg_test_07.items() if v.reject]
failed_seg_test_08 = [k for k,v in seg_test_08.items() if v.reject]
failed_seg_test_09 = [k for k,v in seg_test_09.items() if v.reject]
failed_seg_test_10 = [k for k,v in seg_test_10.items() if v.reject]
failed_seg_test_11 = [k for k,v in seg_test_11.items() if v.reject]
failed_seg_test_12 = [k for k,v in seg_test_12.items() if v.reject]
failed_seg_test_13 = [k for k,v in seg_test_13.items() if v.reject]
failed_seg_test_14 = [k for k,v in seg_test_14.items() if v.reject]

# Epsilon
#hist([x.epsilon for x in stats_01.values()])
#hist([x.epsilon for x in stats_08.values()])
#hist([x.epsilon for x in stats_09.values()])

failed_epsilon_test_01 = [k for k,v in stats_01.items() if v.epsilon > .2]
failed_epsilon_test_02 = [k for k,v in stats_02.items() if v.epsilon > .2]
failed_epsilon_test_03 = [k for k,v in stats_03.items() if v.epsilon > .2]
failed_epsilon_test_04 = [k for k,v in stats_05.items() if v.epsilon > .2]
failed_epsilon_test_05 = [k for k,v in stats_05.items() if v.epsilon > .2]
failed_epsilon_test_06 = [k for k,v in stats_06.items() if v.epsilon > .2]
failed_epsilon_test_07 = [k for k,v in stats_07.items() if v.epsilon > .2]
failed_epsilon_test_08 = [k for k,v in stats_08.items() if v.epsilon > .2]
failed_epsilon_test_09 = [k for k,v in stats_09.items() if v.epsilon > .2] 
failed_epsilon_test_10 = [k for k,v in stats_10.items() if v.epsilon > .2]
failed_epsilon_test_11 = [k for k,v in stats_11.items() if v.epsilon > .2]
failed_epsilon_test_12 = [k for k,v in stats_12.items() if v.epsilon > .2] 
failed_epsilon_test_13 = [k for k,v in stats_13.items() if v.epsilon > .2]
failed_epsilon_test_14 = [k for k,v in stats_14.items() if v.epsilon > .2]

# What is going on with  'Oke_ROA1-209', is family_01??

# missingness per locus
#miss_stats_01 = zip(mappable_01.keys(), [xx.count('-')/float(len(xx)) for xx in mappable_01.values()])
#miss_stats_08 = zip(mappable_08.keys(), [xx.count('-')/float(len(xx)) for xx in mappable_08.values()])
#miss_stats_09 = zip(mappable_09.keys(), [xx.count('-')/float(len(xx)) for xx in mappable_09.values()])


#hist([y for x,y in miss_stats_01])
#hist([y for x,y in miss_stats_08])
#hist([y for x,y in miss_stats_09])

failed_missing_test_01 = [k for k,v in mappable_01.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_02 = [k for k,v in mappable_02.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_03 = [k for k,v in mappable_03.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_04 = [k for k,v in mappable_04.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_05 = [k for k,v in mappable_05.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_06 = [k for k,v in mappable_06.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_07 = [k for k,v in mappable_07.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_08 = [k for k,v in mappable_08.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_09 = [k for k,v in mappable_09.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_10 = [k for k,v in mappable_10.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_11 = [k for k,v in mappable_11.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_12 = [k for k,v in mappable_12.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_13 = [k for k,v in mappable_13.items() if v.count('-')/float(len(v)) > .25]
failed_missing_test_14 = [k for k,v in mappable_14.items() if v.count('-')/float(len(v)) > .25]




def remove_locus(mappable, locus):
    if locus in mappable:
        del mappable[locus]
        return(1)
    else:
        return(0)
    
def remove_catID(mappable, catID):
    removed_count = 0
    if catID+"_x1" in mappable:
        del mappable[ catID+"_x1"]
        removed_count +=1
    if catID+"_x2" in mappable:
        del mappable[ catID+"_x2"]
        removed_count +=1
    return(removed_count)

def remove_catIDs(mappable, catIDs):
    num_removed = 0
    for catID in catIDs:
        num_removed += remove_catID(mappable, catID)
    print("Removed {} loci".format(num_removed))
    
def remove_loci(mappable, loci):
    num_removed = 0
    for locus in loci:
        num_removed += remove_locus(mappable, locus)
    print("Removed {} loci".format(num_removed))  
   
    
remove_catIDs(mappable_01, failed_epsilon_test_01)
remove_catIDs(mappable_02, failed_epsilon_test_02)
remove_catIDs(mappable_03, failed_epsilon_test_03)
remove_catIDs(mappable_04, failed_epsilon_test_04)
remove_catIDs(mappable_05, failed_epsilon_test_05)
remove_catIDs(mappable_06, failed_epsilon_test_06)
remove_catIDs(mappable_07, failed_epsilon_test_07)
remove_catIDs(mappable_08, failed_epsilon_test_08)
remove_catIDs(mappable_09, failed_epsilon_test_09)
remove_catIDs(mappable_10, failed_epsilon_test_10)
remove_catIDs(mappable_11, failed_epsilon_test_11)
remove_catIDs(mappable_12, failed_epsilon_test_12)
remove_catIDs(mappable_13, failed_epsilon_test_13)
remove_catIDs(mappable_14, failed_epsilon_test_14)

remove_catIDs(mappable_01, failed_missing_test_01)
remove_catIDs(mappable_02, failed_missing_test_02)
remove_catIDs(mappable_03, failed_missing_test_03)
remove_catIDs(mappable_04, failed_missing_test_04)
remove_catIDs(mappable_05, failed_missing_test_05)
remove_catIDs(mappable_06, failed_missing_test_06)
remove_catIDs(mappable_07, failed_missing_test_07)
remove_catIDs(mappable_08, failed_missing_test_08)
remove_catIDs(mappable_09, failed_missing_test_09)
remove_catIDs(mappable_10, failed_missing_test_10)
remove_catIDs(mappable_11, failed_missing_test_11)
remove_catIDs(mappable_12, failed_missing_test_12)
remove_catIDs(mappable_13, failed_missing_test_13)
remove_catIDs(mappable_14, failed_missing_test_14)

remove_catIDs(mappable_01, failed_seg_test_01)
remove_catIDs(mappable_02, failed_seg_test_02)
remove_catIDs(mappable_03, failed_seg_test_03)
remove_catIDs(mappable_04, failed_seg_test_04)
remove_catIDs(mappable_05, failed_seg_test_05)
remove_catIDs(mappable_06, failed_seg_test_06)
remove_catIDs(mappable_07, failed_seg_test_07)
remove_catIDs(mappable_08, failed_seg_test_08)
remove_catIDs(mappable_09, failed_seg_test_09)
remove_catIDs(mappable_10, failed_seg_test_10)
remove_catIDs(mappable_11, failed_seg_test_11)
remove_catIDs(mappable_12, failed_seg_test_12)
remove_catIDs(mappable_13, failed_seg_test_13)
remove_catIDs(mappable_14, failed_seg_test_14)

def write_mappable_stats(filename, stats, mappable):
    with open(filename, 'w') as OUTFILE:
        for catID, values in stats.items():
            locally_mappable = sum((catID+"_x1" in mappable, catID+"_x2" in mappable))
            OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                 
                )
            )

write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_01.stats', stats_01, mappable_01)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_02.stats', stats_02, mappable_02)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_03.stats', stats_03, mappable_03)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_04.stats', stats_04, mappable_04)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_05.stats', stats_05, mappable_05)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_06.stats', stats_06, mappable_06)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_07.stats', stats_07, mappable_07)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_08.stats', stats_08, mappable_08)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_09.stats', stats_09, mappable_09)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_10.stats', stats_10, mappable_10)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_11.stats', stats_11, mappable_11)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_12.stats', stats_12, mappable_12)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_13.stats', stats_13, mappable_13)
write_mappable_stats('Y:/DATA/PROJECTS/MYKISS/OmykissDomestication/mykiss_14.stats', stats_14, mappable_14)

# TODO
# Could also filter for difference in likelihoods.

# Write loci to MSTmap format
chum_fams.write_mstmap('CMUW10X_0001', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_01_mstmap.txt', mappable_01)
chum_fams.write_mstmap('CMUW10X_0008', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_08_mstmap.txt', mappable_08)
chum_fams.write_mstmap('CMUW10X_0009', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_09_mstmap.txt', mappable_09)

# Set inital map
chum_fams.set_map(parent = 'CMUW10X_0001', mappable_loci = mappable_01, source_file = None, source_format = None)
chum_fams.set_map(parent = 'CMUW10X_0008', mappable_loci = mappable_08, source_file = None, source_format = None)
chum_fams.set_map(parent = 'CMUW10X_0009', mappable_loci = mappable_09, source_file = None, source_format = None)

# Write loci to Rqtl format
chum_fams.write_Rqtl('CMUW10X_0001', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/rqtl/chum_01_Rqtl.txt', mappable_01)
chum_fams.write_Rqtl('CMUW10X_0008', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/rqtl/chum_08_Rqtl.txt', mappable_08)
chum_fams.write_Rqtl('CMUW10X_0009', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/rqtl/chum_09_Rqtl.txt', mappable_09)


# pickle 'chum_fams' here, as very costly to generate
chum_fams_pickle_file = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/pkl/chum_fams.pkl"
with open(chum_fams_pickle_file, 'wb') as OUTFILE:
    pickle.dump(chum_fams, OUTFILE)

map_01_pickle_file = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/pkl/mappable_01.pkl"
with open(map_01_pickle_file, 'wb') as OUTFILE:
    pickle.dump(mappable_01, OUTFILE)
    
map_08_pickle_file = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/pkl/mappable_08.pkl"
with open(map_08_pickle_file, 'wb') as OUTFILE:
    pickle.dump(mappable_08, OUTFILE)
    
map_09_pickle_file ="/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/pkl/mappable_09.pkl"
with open(map_09_pickle_file, 'wb') as OUTFILE:
    pickle.dump(mappable_09, OUTFILE)
    
## Can't pickle stats_0X files because of namedTuple 
#map_01_stats_pickle_file = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/pkl/stats_01.pkl"
#with open(map_01_stats_pickle_file, 'wb') as OUTFILE:
#    pickle.dump(stats_01, OUTFILE)
#    
#map_08_stats_pickle_file = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/pkl/stats_08.pkl"
#with open(map_08_stats_pickle_file, 'wb') as OUTFILE:
#    pickle.dump(stats_08, OUTFILE)
#    
#map_09_stats_pickle_file = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/pkl/stats_09.pkl"
#with open(map_09_stats_pickle_file, 'wb') as OUTFILE:
#    pickle.dump(stats_09, OUTFILE)
    
##Load Pickle    
INFILE = open(chum_fams_pickle_file, 'rb')
chum_fams = pickle.load(INFILE)
INFILE.close() 


  

########################


frac_missing = numpy.array([xx.count('-')/float(len(xx)) for xx in mappable_08.values()])
frac_pass = frac_missing > .25
[miss_stats[xx] if frac_pass[xx] else None for xx in range(len(frac_missing))]
# Perhaps the code that evaluates which splits catIDs could be made more tolerant, so that one could 'split' into one genotype rather than two

# Print how many *would be* excluded
print("Reject {} out of {} loci in fam_01".format(numpy.sum(numpy.array([xx.count('-')/float(len(xx)) for xx in mappable_01.values()]) > 0.25), len(mappable_01)))
print("Reject {} out of {} loci in fam_08".format(numpy.sum(numpy.array([xx.count('-')/float(len(xx)) for xx in mappable_08.values()]) > 0.25), len(mappable_08)))
print("Reject {} out of {} loci in fam_09".format(numpy.sum(numpy.array([xx.count('-')/float(len(xx)) for xx in mappable_09.values()]) > 0.25), len(mappable_09)))
                            
# Write loci to MSTmap format
chum_fams.write_mstmap('CMUW10X_0001', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_01_mstmap.txt', mappable_01)
chum_fams.write_mstmap('CMUW10X_0008', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_08_mstmap.txt', mappable_08)
chum_fams.write_mstmap('CMUW10X_0009', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_09_mstmap.txt', mappable_09)





# Set inital map
chum_fams.set_map(parent = 'CMUW10X_0001', mappable_loci = mappable_01, source_file = None, source_format = None)
chum_fams.set_map(parent = 'CMUW10X_0008', mappable_loci = mappable_08, source_file = None, source_format = None)
chum_fams.set_map(parent = 'CMUW10X_0009', mappable_loci = mappable_09, source_file = None, source_format = None)

# Write loci to Rqtl format
chum_fams.write_Rqtl('CMUW10X_0001', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/rqtl/chum_01_Rqtl.txt', mappable_01)
chum_fams.write_Rqtl('CMUW10X_0008', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/rqtl/chum_08_Rqtl.txt', mappable_08)
chum_fams.write_Rqtl('CMUW10X_0009', '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/rqtl/chum_09_Rqtl.txt', mappable_09)


#chum fams
chum_fams_pickle_file = "Z:/HOME/Waples/Stacks_mapping/Chum_data/pkl/chum_fams.pkl"
with open(chum_fams_pickle_file, 'wb') as OUTFILE:
    pickle.dump(chum_fams, OUTFILE)

#mappable loci









#
#
#mapped_08 = dict()
#for key, value in mappable_08[0].items():
#    mapped_08[key] = value
#for key, value in mappable_08[1].items():
#    mapped_08[key] = value
#    
##mapped_08 = dict()
##for key, value in switched_mappable_08[0].items():
##    mapped_08[key] = value
##for key, value in switched_mappable_08[1].items():
##    mapped_08[key] = value    
#
#    
##take a map, and examine RF and lod across LGs,
## here i take rqtl map, TODO:use mst map output
#ini_map_filename = "Z:/HOME/Waples/Stacks_mapping/Chum_data/rqtl/chum08_ini_map.tsv"
##ini_map_filename = "Z:/HOME/Waples/Stacks_mapping/Chum_data/rqtl/chum08_switched1_map.tsv"
#loci = list()
#lgs = list()
#loci_on_lg = collections.OrderedDict()
#positions = list() 
#with open (ini_map_filename, 'r') as INFILE:
#    for line in INFILE:
#        locus, lg, pos = line.strip().split()
#        loci.append(locus)
#        lgs.append(lg)
#        if lg in loci_on_lg:
#            loci_on_lg[lg] = loci_on_lg[lg]+1
#        else:
#            loci_on_lg[lg] = 1
#        positions.append(pos)
#ini_map = zip(loci, lgs, positions)
#
#mapped_genotypes = list()
#for locus_tuple in ini_map:
#    locus_name = locus_tuple[0]
#    mapped_genotypes.append(mapped_08[locus_name]) #nested list
#    
#geno_arr = numpy.array(mapped_genotypes)
#my_a = numpy.where(geno_arr == 'a')
#my_b = numpy.where(geno_arr == 'b')
#my_m = numpy.where(geno_arr == '-')
#int_arr = numpy.empty(shape = (geno_arr.shape), dtype = numpy.int)
#int_arr[my_a] = 1 # allele 1
#int_arr[my_b] = 2 # allele 2
#int_arr[my_m] = 0 # missing
#    
#pairs = itertools.combinations(int_arr, 2)
#
#def getR(pairs):
#    for x, y in pairs:
#        mult = x * y
#        yield sum(mult == 2)
#        
#def getNR(pairs):
#    for x, y in pairs:
#        mult = x * y
#        yield (sum(mult == 1) + sum(mult == 4))
#    
#def getM(pairs):
#    for x, y in pairs:
#        mult = x * y
#        yield (sum(mult == 0))
#
#pairs = itertools.combinations(int_arr, 2)
## recombinants
#R = numpy.fromiter(getR(pairs), dtype = numpy.float)
#pairs = itertools.combinations(int_arr, 2)
## non-recombinants
#NR = numpy.fromiter(getNR(pairs), dtype =numpy.float)
#pairs = itertools.combinations(int_arr, 2)
## missing
#M = numpy.fromiter(getM(pairs), dtype =numpy.float)
#
## R_frac
#ml_R_frac = R / (R + NR)
#
## lod score
#Z = log10(
#    (power((1-ml_R_frac), NR) * power(ml_R_frac, R)) / power(.5, (R + NR))
#    )
#
##scipy.spatial.distance.squareform(M)
#rf = scipy.spatial.distance.squareform(ml_R_frac)
#numpy.fill_diagonal(rf, np.nan)
#
#lod = scipy.spatial.distance.squareform(Z)
#numpy.fill_diagonal(lod, np.nan)
#
## to get the mean of a region:
## scipy.stats.nanmean(array_slice.flatten())
#
#
#index_of_lg = collections.OrderedDict()
#index = 0
#for lg, num_loci in loci_on_lg.items():
#    start = index
#    stop = start + num_loci
#    index = stop
#    index_of_lg[lg] = start, stop
#    print lg, start, stop 
#
#    
#def get_mean(pairs, array):
#    for x, y in pairs:
#        lg_a_start = index_of_lg[x][0]
#        lg_a_stop = index_of_lg[x][1]
#        lg_b_start = index_of_lg[y][0]
#        lg_b_stop = index_of_lg[y][1]
#        mean = numpy.mean(array[lg_a_start:lg_a_stop, lg_b_start:lg_b_stop].flatten())
#        yield (mean)
#        
#def get_sum(pairs, array):
#    for x, y in pairs:
#        lg_a_start = index_of_lg[x][0]
#        lg_a_stop = index_of_lg[x][1]
#        lg_b_start = index_of_lg[y][0]
#        lg_b_stop = index_of_lg[y][1]
#        my_sum = numpy.sum(array[lg_a_start:lg_a_stop, lg_b_start:lg_b_stop].flatten())
#        yield (my_sum)        
#        
#        
#ordered_lgs = index_of_lg.keys()
#lgs_longer_than_1 = list()
#for xx in ordered_lgs:
#    if loci_on_lg[xx] > 1:
#        lgs_longer_than_1.append(xx)
#
## all lgs
##pairs_of_lgs = itertools.combinations(ordered_lgs, 2)        
##mean_rf = numpy.fromiter(get_mean(pairs_of_lgs, rf), dtype = numpy.float)
##pairs_of_lgs = itertools.combinations(ordered_lgs, 2)
##mean_lod = numpy.fromiter(get_mean(pairs_of_lgs, lod), dtype = numpy.float)
##pairs_of_lgs = itertools.combinations(ordered_lgs, 2)
##sum_lod = numpy.fromiter(get_sum(pairs_of_lgs, lod), dtype = numpy.float)
#
##lgs with more than 1 locus
#pairs_of_lgs = itertools.combinations(lgs_longer_than_1, 2)        
#mean_rf = numpy.fromiter(get_mean(pairs_of_lgs, rf), dtype = numpy.float)
#pairs_of_lgs = itertools.combinations(lgs_longer_than_1, 2)
#mean_lod = numpy.fromiter(get_mean(pairs_of_lgs, lod), dtype = numpy.float)
#pairs_of_lgs = itertools.combinations(lgs_longer_than_1, 2)
#sum_lod = numpy.fromiter(get_sum(pairs_of_lgs, lod), dtype = numpy.float)
#
##change to square_matrix
#sq_sum_lod = scipy.spatial.distance.squareform(sum_lod)
#upper_tri = numpy.triu_indices(len(lgs_longer_than_1))
#sq_sum_lod[upper_tri] = 0
#
#where_pairs = numpy.where(sq_sum_lod>15000)
#
#   
#    
#
##Select LG that appears in most switched pairs, if tied pick one
## Add that LG to switch list
## remove all pairs contianing that LG from the initial list
## repeat until no pairs remain, return switch list
#
#def find_LGs_to_switch(where_pairs):
#    matched_pairs = zip(*where_pairs)
#    LGs_to_switch = list()
#    while len(matched_pairs) > 0:
#        count_switched = collections.defaultdict(int)
#        for ii, jj in matched_pairs:
#            count_switched[ii] += 1
#            count_switched[jj] += 1
#        sorted_switched = sorted(count_switched.items(),key=operator.itemgetter(1), reverse = True)
#        most_common_LG = sorted_switched.pop()[0]
#        LGs_to_switch.append(most_common_LG)
#        to_remove = list()
#        for index in range(len(matched_pairs)):
#            xx, yy = matched_pairs[index]
#            if most_common_LG == xx or  most_common_LG == yy:
#                to_remove.append((xx, yy))
#        for ii in to_remove:    
#            matched_pairs.remove(ii)
#    return(LGs_to_switch)
#    
#LGs_to_switch = find_LGs_to_switch(where_pairs)
#
## Now to turn this list of LGs into a a list of locus names
#list_of_loci_on_LG = collections.defaultdict(list)                                                            
#for locus, LG, pos in ini_map:
#    #convert to int
#    LG = int(LG)
#    list_of_loci_on_LG[LG].append(locus)
#    
#loci_to_switch = list()
#for LG in LGs_to_switch:
#    loci_to_switch += list_of_loci_on_LG[LG]
#    
#loci_to_switch    
#    
#    
#def switch_alleles(loci_to_switch, mappable):
#    x1, x2 = mappable
#    for jj in loci_to_switch:
#        if jj in x1:
#            original_genotypes = x1[jj]
#            switched_genotpyes = ['a' if xx is 'b' else 'b' if xx is 'a' else xx for xx in original_genotypes]
#            x1[jj] = switched_genotpyes
#    for jj in loci_to_switch:
#        if jj in x2:
#            original_genotypes = x2[jj]
#            switched_genotpyes = ['a' if xx is 'b' else 'b' if xx is 'a' else xx for xx in original_genotypes]
#            x2[jj] = switched_genotpyes
#    return(x1, x2)
#
#switched_mappable_08 = switch_alleles(loci_to_switch, mappable_08)
#    
#chum_fams.write_Rqtl('CMUW10X_0008', 'Z:/HOME/Waples/Stacks_mapping/Chum_data/psv/chum_08_switched_Rqtl.txt', switched_mappable_08[0], switched_mappable_08[1])
#chum_fams.write_mstmap('CMUW10X_0008', 'Z:/HOME/Waples/Stacks_mapping/Chum_data/psv/chum_08_switched_mstmap.txt', switched_mappable_08[0], switched_mappable_08[1])
#    
#    
#    
#        
#len(find_LGs_to_switch(where_pairs))
#
#
#
#print(len(find_LGs_to_switch(where_pairs, False)))
#
#
#
#
#selected_pairs = zip(*numpy.where(sq_sum_lod>10000))
#for 
#
#
#
#
#
#
#max(sum_lod)
#min(sum_lod)
#numpy.where(sum_lod>1000)
#
#scatter(mean_rf, sum_lod)
#                
#scipy.spatial.distance.squareform(sum_lod)
#
#
####################    
#
#        
#    
#    
#
#
##get pairwise R fraction
#mapr = mapped_08[0].values() # just singly mappable loci for now
##nl = [[ii for ii in xx] for xx in mapr.values()]
#
#my_arr = numpy.array(mapr)
#my_a = numpy.where(my_arr == 'a')
#my_b = numpy.where(my_arr == 'b')
#my_m = numpy.where(my_arr == '-')
#int_arr = numpy.empty(shape = (my_arr.shape), dtype = numpy.int)
#int_arr[my_a] = 1 # allele 1
#int_arr[my_b] = 2 # allele 2
#int_arr[my_m] = 0 # missing
#
## remove values where all missing, should likely be done earlier
#aa = numpy.array([sum(x)> 0 for x in int_arr])
#cleaned_arr = int_arr[aa]
#
#
#pairs = itertools.combinations(cleaned_arr, 2)
#
#def getR(pairs):
#    for x, y in pairs:
#        mult = x * y
#        yield sum(mult == 2)
#        
#def getNR(pairs):
#    for x, y in pairs:
#        mult = x * y
#        yield (sum(mult == 1) + sum(mult == 4))
#    
#def getM(pairs):
#    for x, y in pairs:
#        mult = x * y
#        yield (sum(mult == 0))
#
#pairs = itertools.combinations(cleaned_arr, 2)
## recombinants
#R = numpy.fromiter(getR(pairs), dtype = numpy.float)
#pairs = itertools.combinations(cleaned_arr, 2)
## non-recombinants
#NR = numpy.fromiter(getNR(pairs), dtype =numpy.float)
#pairs = itertools.combinations(cleaned_arr, 2)
## missing
#M = numpy.fromiter(getM(pairs), dtype =numpy.float)
#
#ml_R_frac = R / (R + NR)
#
#Z = log10(
#    (power((1-ml_R_frac), NR) * power(ml_R_frac, R)) / power(.5, (R + NR))
#    )
#switched_Z = -log10(
#    (power(1-ml_R_frac, R) * power((ml_R_frac), NR)) / power(.5, (R + NR))
#    )
#
#hexbin(Z, switched_Z,
#mincnt = 1)
#matplotlib.pyplot.xlim([-0.04, 1.04])
##matplotlib.pyplot.ylim([-0.3, max(Z)*1.04])
#cb = matplotlib.pyplot.colorbar()
#cb.set_label('count')
#matplotlib.pyplot.show()
#
#
#
#
## convert to square, check diagonal
#scipy.spatial.distance.squareform(R).shape
#scipy.spatial.distance.squareform(NR).shape
#scipy.spatial.distance.squareform(M).shape
#
#
## really slow
#R = numpy.empty(shape = (len(int_arr), len(int_arr)), dtype = numpy.int) #recombinants
#NR = numpy.empty(shape = (len(int_arr), len(int_arr)), dtype = numpy.int) # non
#for ii in range(len(int_arr)):
#    for jj in range(len(int_arr)):
#        if jj >= ii:
#            mult = int_arr[ii] * int_arr[jj]
#            R[ii, jj] = sum(mult == 2)
#            NR[ii, jj] = sum(mult == 1) + sum(mult == 4)
#    
#
#
#
#
#
#chum_fams.find_psvs(parent = 'CMUW10X_0001', outfile = 'Z:/HOME/Waples/Stacks_mapping/Chum_data/psv/chum_01.psvlikes.tsv')
#chum_fams.find_psvs(parent = 'CMUW10X_0008', outfile = 'Z:/HOME/Waples/Stacks_mapping/Chum_data/psv/chum_08.psvlikes.tsv')
#chum_fams.find_psvs(parent = 'CMUW10X_0009', outfile = 'Z:/HOME/Waples/Stacks_mapping/Chum_data/psv/chum_09.psvlikes.tsv')