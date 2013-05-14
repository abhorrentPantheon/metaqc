''' metaqc.gcqc.proc.py
'''

import os
import re
import sys
import time

sys.path.append('/x/PyMS/')
sys.path.append('/x/metaqc/')

# From common_functions
from pylib.Function import config_reader

# From current directory
from gcqc.Class import Compound
from gcqc.Function import file2matrix, target_reader, get_peak_list, \
    qc_db_update_badfile, qc_db_update, db_plotter, report_gen, \
    qc_report_badfile, qc_report_tbs, qc_report, symm_test

from math import fabs
from pycdf import CDF

# From original pyms (i.e. svn)
from pyms.Peak.Class import Peak

#==========================================================================
#    Parse arguments from command line 
#==========================================================================
## FUTUREPY: argparse is used in python3 (also available >= 2.7)
## Tested. Differences to note are that argparse defaults to type=string and
##    action='store'. nargs defaults to the next argument for -x or --xxx 
##    and the argument itself in the case of positional arguments. Specifying
##    nargs here will cause the parser to store the result as a list of 
##    strings, not as a string.
##    Further, the positional argument is different. It is automatically stored
##    as a string with the name specified in the 'name or flag' parameter. To
##    display this differently in help, use metavar to change this.
#import argparse
#parser = argparse.ArgumentParser( \
#    description='Gas Chromatography Quality Control (GCQC) Report creator')
#parser.add_argument('-d', '--database', #nargs=1, 
#    dest='db_filename', metavar='FILE',
#    help='Specify the filename of the database to use')
#parser.add_argument('-m', '--method', #nargs=1,
#    type=int, dest='method_val', metavar='[7|15|25]',
#    help='Specify either 7, 15 or 25 degree method')
#parser.add_argument('-t', '--targets', #nargs=1,
#    dest='target_filename', metavar='FILE.csv',
#    help='Compounds to use for QC method')
#parser.add_argument('-c', '--config', #nargs=1,
#    dest='config_file', metavar='FILE.cfg',
#    help='''Configuration file to use. Specifying other options
#    on the command line overrides settings in the config file''')
#parser.add_argument('andi_file', #nargs=1,
#    metavar='GCMS_file.cdf', 
#    help='GCMS data file to process (a .cdf file)')

## To display help:
##args = parser.print_help()

## To run from cmd line:
#args = parser.parse_args()

#==========================================================================

## Only valid for python2 (deprecated in python3)
from optparse import OptionParser
#
usage = 'usage: %prog [options] GCMS_file.cdf'
parser = OptionParser(usage, description=\
    'Gas Chromatography Quality Control (GCQC) Report creator')
parser.add_option('-d', '--database', nargs=1, action='store', 
    type='string', dest='db_filename', metavar='FILE',
    help='Specify the filename of the database to use')
parser.add_option('-m', '--method', nargs=1, action='store',
    type='string', dest='method_val', metavar='[7|15|25|TBS]',
    help='Specify either 7, 15 or 25 degree method')
parser.add_option('-t', '--targets', nargs=1, action='store',
    type='string', dest='target_filename', metavar='FILE.csv',
    help='Compounds to use for QC method')
parser.add_option('-c', '--config', nargs=1, action='store',
    type='string', dest='config_file', metavar='FILE.cfg',
    help=''.join(['Configuration file to use. Specifying other options ', 
    'on the command line overrides settings in the config file']))
#TODO: Add verbose flag?

#==========================================================================
#    Specifying options
#==========================================================================
#
#    To display help:
#
#(options, args) = parser.parse_args(['-h'])

#
#    Methods for options (within an IDE, single run, cut/paste, etc)
#
#(options, args) = parser.parse_args( [\
## If using alkanes:
##    ['-d', 'gcqc_bot7.db', '-t', 'gcqc_cpd_list_7.csv', '-m', '7', 
##    'datafiles/botany_7/201211C_DD_QC5 CHECKROBIN_MA7C.CDF'])
## Or use the defaults and let it pick the methods based on .cdf file
## (must specify database and .cdf at a minimum)
#    '-d', '/x/metaqc/metaqc_gcqc.db',
#    '-c', '/x/metaqc/default.cfg', 
#    ### Test using a blank file (all fail)
#    'datafiles/tmp/reagent_148-16508.CDF' ])
#    ### Test using a real file
#    'datafiles/sepsis/iqc_161-17207.CDF'])
##    To use actual command line options (normal use):
(options, args) = parser.parse_args()
#==========================================================================

# TODO: modify for optparse/argparse
# Check that an argument (.cdf file name) was provided
### NB: andi_file must be specified as full path (it is by picker)
if len(args) == 0:
    print ' !! No .cdf file specified; exiting.'
    raise RuntimeError
elif len(args) > 1:
    print ' !! More than one argument (.cdf file) specified; exiting.'
    raise RuntimeError
else:
    andi_file = args[0]

# Get instrument name from file path
# Bio21
if 'GCMS01' in andi_file:
    instr_name = 'GC01'
elif 'GCMS02' in andi_file:
    instr_name = 'GC02'
# Botany
#if 'robin' in andi_file:
#    instr_name = 'robin'
#elif 'batman' in andi_file:
#    instr_name = 'batman'

# Get values in config file if specified
if options.config_file != None:
    config_set = config_reader(options.config_file)
    c_ind = config_set[0]['general']
    l_ind = config_set[0]['test_names']
    v_ind = config_set[0]['pf_vals']
    
    config_vals = config_set[c_ind]
    report_labels = config_set[l_ind]
    report_vals = config_set[v_ind]
    
    # Reformat strings for report vals
    for key in report_vals.keys():
        report_vals[key] = float(report_vals[key])    
    
    # Try to convert method_temp to int, if it exists
    try:
        config_vals['method_temp'] = int(config_vals['method_temp'])
    except KeyError:
        pass
    
# If config file not used, set defaults here.
else:
    config_vals = {'target_filename':'', 'db_file':'', 'report_folder':''}
    
    report_labels = {'ala_i':'Alanine intensity', \
        'ar_r':'Asparagine/Ribose ratio', \
        'av_r':'Alanine/Valine ratio', \
        'gr_r':'Glutamate/Ribose ratio', \
        'mlt_i':'Maltotriose intensity', \
        'pr_r':'Putrescine/Ribitol ratio', \
        'rbtl_i':'Ribitol intensity', \
        'rbose_i':'Ribose intensity', \
        'rt_lock':'Retention Time Lock delta'}
    
    report_vals = {'ala_i_fail':1250000.0, \
        'ala_i_pass':1500000.0, \
        'ar_r_fail':10.0, \
        'ar_r_pass':15.0, \
        'av_r_fail':0.5, \
        'av_r_pass':0.6, \
        'gr_r_fail':10.0, \
        'gr_r_pass':15.0, \
        'mlt_i_fail':5000.0, \
        'mlt_i_pass':7500.0, \
        'pr_r_fail':50.0, \
        'pr_r_pass':70.0, \
        'rbtl_i_fail':50000.0, \
        'rbtl_i_pass':120000.0, \
        'rbose_i_fail':25000.0, \
        'rbose_i_pass':55000.0}

# If values specified on cmd line, override config values
if options.db_filename == None: 
    # Check if a database file name was provided in config file
    if config_vals['db_file'] == '':
        print ' !! No database file specified; exiting.'
        raise RuntimeError
    else:
        db_file = config_vals['db_file']
else:
    db_file = options.db_filename

# Determine maximum time value in file 
# Used to expand time window for early-eluting compounds
cdf_time = CDF(andi_file).var('scan_acquisition_time')
time_list = cdf_time.get().tolist()
max_time = max(time_list)

#
#    Method file used in experiment file
#
# Check method_val for a valid number
method_filename = CDF(andi_file).attributes()['external_file_ref_0']

if options.method_val in ['7', '15', '25', 'TBS']:
    try:
        method_temp = int(options.method_val)
    except:
        method_temp = options.method_val
    
# If no method specified, try to get info from file
else:
    try:
        if config_vals['method_temp'] in ['7', '15', '25', 'TBS']:
            try:
                method_temp = int(config_vals['method_temp'])
            except:
                method_temp = config_vals['method_temp']
    except:
        # TODO: Ensure this covers everyone's method filenames
        # TODO: Scan for TBS first, if so, break and use that
        if 'TBS' in method_filename:
            method_temp = 'TBS'
        else:
            try:
                method_temp = int(re.sub('[^\\d]*(\\d+)[^\\d]*', \
                    '\\1', method_filename))
                # If this produces a number, check if it is valid
                if method_temp not in [7, 15, 25]:
                    ### 25 deg method run time = 18.6m = 1116s
                    if max_time < 1250:
                        method_temp = 25
                        #rt_window = 0.06
                    ### 15 deg method run time = 24.3m = 1460s
                    elif 1250 < max_time < 1650:
                        method_temp = 15
                        #rt_window = 0.10
                    ### 7 deg method run time = 41.0m = 2461s
                    else:
                        method_temp = 7
                        #rt_window = 0.15       
            
            # If this doesn't produce an integer, raises ValueError
            # Attempt to determine run length, and set value based on that
            except ValueError:
                ### 25 deg method run time = 18.6m = 1116s
                if max_time < 1250:
                    method_temp = 25
                    #rt_window = 0.06
                ### 15 deg method run time = 24.3m = 1460s
                elif 1250 < max_time < 1650:
                    method_temp = 15
                    #rt_window = 0.10
                ### 7 deg method run time = 41.0m = 2461s
                else:
                    method_temp = 7
                    #rt_window = 0.15

if config_vals['report_folder'] == '':
    report_folder = './reports'
    if not os.path.isdir(report_folder):
        os.mkdir(report_folder)
else:
    report_folder = config_vals['report_folder']
    if not os.path.isdir(report_folder):
        os.makedirs(report_folder)

#
#    Compounds for QC tests
#
# If no target file specified, use these as defaults
# Compound class specified as ('Name', RT, MST1)
# TODO: Leave these here or only have in a config?
if options.target_filename == None:
    if config_vals['target_filename'] == '':
        if method_temp == 25:
            ala2 = Compound('Alanine 2TMS', 6.9, 190)               # main
            gly2 = Compound('Glycine 2TMS', 7.031, 204)
            val2 = Compound('Valine 2TMS', 7.557, 144)              # main
            c13val = Compound('13C,15N-Valine', 7.557, 149)
            ser2 = Compound('Serine 2TMS', 7.806, 116)
            leu2 = Compound('Leucine 2TMS', 7.868, 158)
            ile2 = Compound('Isoleucine 2TMS', 8.006, 158)
            thr2 = Compound('Threonine 2TMS', 8.015, 117)
            pro1 = Compound('Proline 1TMS', 8.068, 142)
            gly3 = Compound('Glycine 3TMS', 8.095, 248)
            succ2 = Compound('Succinate 2TMS', 8.112, 247)
            norleu2 = Compound('Norleucine 2TMS', 8.134, 158)
            ser3 = Compound('Serine 3TMS', 8.334, 204)
            thr3 = Compound('Threonine 3TMS', 8.474, 218)
            asp2 = Compound('Aspartate 2TMS', 9.176, 245)
            met2 = Compound('Methionine 2TMS', 9.232, 293)
            pyrglu2 = Compound('Pyroglutamate 2TMS', 9.268, 156)
            cys3 = Compound('Cysteine 3TMS', 9.3823, 294)
            phe1 = Compound('Phenylalanine', 9.458, 120)
            glu3 = Compound('Glutamate 3TMS', 9.663, 246)           # main
            phe2 = Compound('Phenylalanine 2TMS', 9.778, 192)
            rbose = Compound('Ribose', 9.886, 307)                  # main
            asn3 = Compound('Asparagine 3TMS', 9.92, 188)           # main
            rbtl = Compound('Ribitol', 10.092, 319)                 # main
            fucs1 = Compound('Fucose MX1', 10.1374, 117)
            fucs2 = Compound('Fucose MX2', 10.2099, 117)
            putr4 = Compound('Putrescine 4TMS', 10.261, 214)        # main
            isocitr = Compound('Isocitric acid', 10.542, 245)
            citr = Compound('Citric acid', 10.543, 183)
            orn = Compound('Ornithine', 10.571, 142)
            fruc1 = Compound('Fructose MX1', 10.737, 307)
            lys3b = Compound('Lysine 3TMS pk2', 10.767, 258)
            fruc2 = Compound('Fructose MX2', 10.778, 307)
            gluc1 = Compound('Glucose MX1', 10.857, 319)
            tyr2 = Compound('Tyrosine 2TMS', 10.939, 179)
            gluc2 = Compound('Glucose MX2', 10.954, 319)
            mann = Compound('Mannitol', 10.99, 319)                 # main
            lys4 = Compound('Lysine 4TMS', 11.007, 174)
            c13sorb = Compound('13C Sorbitol', 11.007, 323)
            his = Compound('Histidine', 11.06, 154)
            tyr3 = Compound('Tyrosine 3TMS', 11.106, 280)
            mino6 = Compound('Myoinositol 6TMS', 11.658, 305)
            rbos5p1 = Compound('Ribose-5-phosphate MX1', 11.668, 357)
            rbos5p2 = Compound('Ribose-5-phosphate MX2', 11.695, 357)
            trp2 = Compound('Tryptophan 2TMS', 12.176, 231)
            trp3 = Compound('Tryptophan 3TMS', 12.2, 305)
            fruc6p = Compound('Fructose-6-Phosphate', 12.383, 217)
            cys4 = Compound('Cysteine 4TMS', 12.408, 411)
            gluc6p1 = Compound('Glucose-6-Phosphate MX1', 12.435, 299)
            gluc6p2 = Compound('Glucose-6-Phosphate MX2', 12.515, 299)
            trhl8 = Compound('Trehalose 8TMS', 13.774, 361)
            lano = Compound('Lanosterol', 16.6592, 498)
            mlt1 = Compound('Maltotriose MX1', 17.333, 361)         # main
            mlt2 = Compound('Maltotriose MX2', 17.602, 361)
            
            cpd_list = [ala2, gly2, val2, c13val, ser2, leu2, ile2, thr2, \
                pro1, gly3, succ2, norleu2, ser3, thr3, asp2, met2, pyrglu2, \
                cys3, phe1, glu3, phe2, rbose, asn3, rbtl, fucs1, fucs2, \
                putr4, isocitr, citr, orn, fruc1, lys3b, fruc2, gluc1, tyr2, \
                gluc2, mann, lys4, c13sorb, his, tyr3, mino6, rbos5p1, \
                rbos5p2, trp2, trp3, fruc6p, cys4, gluc6p1, gluc6p2, trhl8, \
                lano, mlt1, mlt2]
        
        elif method_temp == 15:
            pyr1 = Compound('Pyruvate MX1', 8.3326, 174)
            ala2 = Compound('Alanine 2TMS', 8.7989, 190)
            gly2 = Compound('Glycine 2TMS', 8.9729, 204)
            val2 = Compound('Valine 2TMS', 9.8352, 144)
            c13val = Compound('13C,15N-Valine', 9.8352, 149)
            ser2 = Compound('Serine 2TMS', 10.2328, 116)
            leu2 = Compound('Leucine 2TMS', 10.3421, 158)
            ile2 = Compound('Isoleucine 2TMS', 10.547, 158)
            pro1 = Compound('Proline 1TMS', 10.6311, 142)
            gly3 = Compound('Glycine 3TMS', 10.686, 248)
            succ2 = Compound('Succinate 2TMS', 10.7392, 247)
            norleu2 = Compound('Norleucine 2TMS', 10.7741, 158)
            ser3 = Compound('Serine 3TMS', 11.11, 204)
            thr3 = Compound('Threonine 3TMS', 11.3372, 218)
            met1 = Compound('Methionine 1TMS', 11.6908, 104)
            asp2 = Compound('Aspartate 2TMS', 11.7728, 245)
            asp3 = Compound('Aspartate 3TMS', 12.4947, 232)
            met2 = Compound('Methionine 2TMS', 12.5428, 293)
            pyrglu2 = Compound('Pyroglutamate 2TMS', 12.6012, 156)
            cys3 = Compound('Cysteine 3TMS', 12.8003, 294)
            akg = Compound('alpha-Ketoglutaric acid', 12.9473, 198)
            pro2 = Compound('Proline 2TMS', 13.0583, 142)
            orn3 = Compound('Ornithine 3TMS', 13.2506, 162)
            glu3 = Compound('Glutamate 3TMS', 13.2889, 246)
            phe2 = Compound('Phenylalanine 2TMS', 13.4268, 192)
            rbose = Compound('Ribose', 13.6712, 307)
            asn3 = Compound('Asparagine 3TMS', 13.6986, 188)
            rbtl = Compound('Ribitol', 14.0276, 319)
            fucs1 = Compound('Fucose MX1', 14.0585, 117)
            fucs2 = Compound('Fucose MX2', 14.1684, 117)
            putr4 = Compound('Putrescine 4TMS', 14.2325, 214)
            gln4 = Compound('Glutamine 4TMS', 14.4637, 347)
            isocitr = Compound('Isocitric acid', 14.734, 245)
            citr = Compound('Citric acid', 14.735, 183)
            orn = Compound('Ornithine', 14.7612, 142)
            lys3b = Compound('Lysine 3TMS pk2', 15.061, 258)
            fruc1 = Compound('Fructose MX1', 15.0708, 307)
            fruc2 = Compound('Fructose MX2', 15.1354, 307)
            gluc1 = Compound('Glucose MX1', 15.2544, 319)
            tyr2 = Compound('Tyrosine 2TMS', 15.3323, 208)
            gluc2 = Compound('Glucose MX2', 15.4021, 205)
            mann = Compound('Mannitol', 15.4776, 319)
            c13sorb = Compound('13C Sorbitol', 15.4805, 323)
            lys4 = Compound('Lysine 4TMS', 15.4833, 174)
            his = Compound('Histidine', 15.5199, 154)
            tyr3 = Compound('Tyrosine 3TMS', 15.6235, 280)
            mino6 = Compound('Myoinositol 6TMS', 16.5539, 305)
            rbos5p1 = Compound('Ribose-5-phosphate MX1', 16.6065, 357)
            rbos5p2 = Compound('Ribose-5-phosphate MX2', 16.644, 357)
            trp2 = Compound('Tryptophan 2TMS', 17.3412, 231)
            trp3 = Compound('Tryptophan 3TMS', 17.4145, 305)
            cys4 = Compound('Cysteine 4TMS', 17.7921, 411)
            gluc6p1 = Compound('Glucose-6-Phosphate MX1', 17.8562, 299)
            gluc6p2 = Compound('Glucose-6-Phosphate MX2', 17.979, 299)
            trhl8 = Compound('Trehalose 8TMS', 20.0569, 361)
            mlt1 = Compound('Maltotriose MX1', 24.2288, 361)
            
            cpd_list = [pyr1, ala2, gly2, val2, c13val, ser2, leu2, ile2, \
                pro1, gly3, succ2, norleu2, ser3, thr3, met1, asp2, asp3, \
                met2, pyrglu2, cys3, akg, pro2, orn3, glu3, phe2, rbose, \
                asn3, rbtl, fucs1, fucs2, putr4, gln4, isocitr, citr, orn, \
                lys3b, fruc1, fruc2, gluc1, tyr2, gluc2, mann, c13sorb, \
                lys4, his, tyr3, mino6, rbos5p1, rbos5p2, trp2, trp3, \
                cys4, gluc6p1, gluc6p2, trhl8, mlt1]
        
        elif method_temp == 7:
            ala2 = Compound('Alanine 2TMS', 7.6079, 190)            # main
            gly2 = Compound('Glycine 2TMS', 7.9569, 204)
            val2 = Compound('Valine 2TMS', 9.64837, 144)            # main
            c13val = Compound('13C,15N-Valine', 9.64837, 149)
            ser2 = Compound('Serine 2TMS', 10.443, 116)
            leu2 = Compound('Leucine 2TMS', 10.6887, 158)
            ile2 = Compound('Isoleucine 2TMS', 11.0975, 158)
            thr2 = Compound('Threonine 2TMS', 11.1255, 117)
            pro1 = Compound('Proline 1TMS', 11.215, 142)
            gly3 = Compound('Glycine 3TMS', 11.3427, 248)
            succ2 = Compound('Succinate 2TMS', 11.5179, 247)
            norleu2 = Compound('Norleucine 2TMS', 11.5823, 158)
            ser3 = Compound('Serine 3TMS', 12.3110, 204)
            thr3 = Compound('Threonine 3TMS', 12.7748, 218)
            asp2 = Compound('Aspartate 2TMS', 13.5456, 245)
            gln3 = Compound('Glutamine 3TMS', 14.4139, 229)
            putr3 = Compound('Putrescine 3TMS', 15.0551, 304)
            asp3 = Compound('Aspartate 3TMS', 15.1949, 232)
            met2 = Compound('Methionine 2TMS', 15.2210, 293)
            pyrglu2 = Compound('Pyroglutamate 2TMS', 15.2934, 156)
            cys3 = Compound('Cysteine 3TMS', 15.797, 294)
            phe1 = Compound('Phenylalanine', 16.1, 120)
            akg = Compound('alpha-Ketoglutaric acid', 16.1591, 198)
            orn3 = Compound('Ornithine 3TMS', 16.7809, 162)
            glu3 = Compound('Glutamate 3TMS', 16.8722, 246)         # main
            phe2 = Compound('Phenylalanine 2TMS', 17.029, 192)
            rbose = Compound('Ribose', 17.5, 307)                   # main
            asn3 = Compound('Asparagine 3TMS', 17.7046, 188)        # main
            lys3a = Compound('Lysine 3TMS pk1', 18.2672, 230)
            rbtl = Compound('Ribitol', 18.4521, 319)                # main
            fucs1 = Compound('Fucose MX1', 18.4796, 117)
            fucs2 = Compound('Fucose MX2', 18.6808, 117)
            putr4 = Compound('Putrescine 4TMS', 18.7296, 214)       # main
            gln4 = Compound('Glutamine 4TMS', 19.3203, 347)
            orn4 = Compound('Ornithine 4TMS', 19.9296, 420)
            citr = Compound('Citric acid', 19.93015, 183)
            icitr = Compound('Isocitric acid', 19.9374, 245)
            lys3b = Compound('Lysine 3TMS pk2', 20.50565, 258)
            fruc1 = Compound('Fructose MX1', 20.6538, 307)
            fruc2 = Compound('Fructose MX2', 20.7974, 307)
            gluc1 = Compound('Glucose MX1', 21.0234, 319)
            tyr2 = Compound('Tyrosine 2TMS', 21.04575, 208)
            gluc2 = Compound('Glucose MX2', 21.30165, 205)
            lys4 = Compound('Lysine 4TMS', 21.46197, 174)
            mann = Compound('Mannitol', 21.46375, 319)              # main
            his = Compound('Histidine', 21.48895, 154)
            c13sorb = Compound('13C Sorbitol', 21.55495, 323)
            tyr3 = Compound('Tyrosine 3TMS', 21.727, 280)
            mino6 = Compound('Myoinositol 6TMS', 23.7003, 305)
            rbos5p1 = Compound('Ribose-5-phosphate MX1', 23.9253, 357)
            rbos5p2 = Compound('Ribose-5-phosphate MX2', 23.9842, 357)
            trp2 = Compound('Tryptophan 2TMS', 25.2276, 231)
            trp3 = Compound('Tryptophan 3TMS', 25.4253, 305)
            cys4 = Compound('Cysteine 4TMS', 26.3286, 411)
            fruc6p = Compound('Fructose-6-Phosphate', 26.45, 217)
            gluc6p1 = Compound('Glucose-6-Phosphate MX1', 26.534, 299)
            gluc6p2 = Compound('Glucose-6-Phosphate MX2', 26.767, 299)
            trhl8 = Compound('Trehalose 8TMS', 31.17497, 361)
            lano = Compound('Lanosterol', 37.06965, 498)
            mlt1 = Compound('Maltotriose MX1', 38.3011, 361)        # main
            mlt2 = Compound('Maltotriose MX2', 38.62, 361)
            
            cpd_list = [ala2, gly2, val2, c13val, ser2, leu2, ile2, thr2, \
                pro1, gly3, succ2, norleu2, ser3, thr3, asp2, gln3, putr3, \
                asp3, met2, pyrglu2, cys3, phe1, akg, orn3, glu3, phe2, \
                rbose, asn3, lys3a, rbtl, fucs1, fucs2, putr4, gln4, orn4, \
                citr, icitr, lys3b, fruc1, fruc2, gluc1, tyr2, gluc2, lys4, \
                mann, his, c13sorb, tyr3, mino6, rbos5p1, rbos5p2, trp2, \
                trp3, cys4, fruc6p, gluc6p1, gluc6p2, trhl8, lano, mlt1, mlt2]
        
        elif method_temp == 'TBS':
            ala2 = Compound('Alanine 2TBS', 8.7861, 232)
            gly = Compound('Glycine', 9.0313, 218)
            val2 = Compound('Valine 2TBS', 9.8751, 288)
            leu2 = Compound('Leucine 2TBS', 10.2167, 200)
            ser2 = Compound('Serine 2TBS', 10.3374, 276)
            succ2 = Compound('Succinate 2TBS', 10.8134, 289)
            oxla = Compound('Oxaloacetate', 10.833, 259)
            ile = Compound('Isoleucine', 10.509, 200)
            norleu = Compound('Norleucine', 10.653, 200)
            succ = Compound('Succinate', 10.8134, 289)
            pro = Compound('Proline', 10.9166, 184)
            phe1 = Compound('Phenylalanine', 11.2194, 222)
            pyr = Compound('Pyroglutamate', 12.4185, 300)
            met2 = Compound('Methionine 2TBS', 12.5286, 292)
            ser3 = Compound('Serine 3TBS', 12.6237, 390)
            glu2 = Compound('Glutamate 2TBS', 12.7282, 318)
            thr = Compound('Threonine', 12.8645, 404)
            asn2 = Compound('Asparagine 2TBS', 13.2486, 201)
            phe2 = Compound('Phenylalanine 2TBS', 13.4963, 308)
            oxlb = Compound('Oxaloacetate TBS', 13.542, 419)
            asp3 = Compound('Aspartate 3TBS', 13.8741, 418)
            cys = Compound('Cysteine', 14.2664, 378)
            glu3 = Compound('Glutamate 3TBS', 14.7142, 432)
            orn3 = Compound('Ornithine 3TBS', 14.7737, 184)
            asn3 = Compound('Asparagine 3TBS', 14.9419, 417)
            lys3 = Compound('Lysine 3TBS', 15.458, 300)
            gln3 = Compound('Glutamine 3TBS', 15.757, 431)
            his3 = Compound('Histidine 3TBS', 16.8273, 196)
            citr = Compound('Citrate', 16.8642, 459)
            icitr = Compound('Isocitrate', 16.937, 591)
            tyr3 = Compound('Tyrosine 3TBS', 17.1094, 302)
            trp2 = Compound('Tryptophan 2TBS', 17.516, 302)
            
            cpd_list = [ala2, gly, val2, leu2, ser2, succ2, oxla, ile, \
                norleu, succ, pro, phe1, pyr, met2, ser3, glu2, thr, asn2, \
                phe2, oxlb, asp3, cys, glu3, orn3, asn3, lys3, gln3, his3, \
                citr, icitr, tyr3, trp2]
            
    else:
        cpd_matrix = file2matrix(config_vals['target_filename'])
        cpd_list = target_reader(cpd_matrix)
    
else:
    cpd_matrix = file2matrix(options.target_filename)
    cpd_list = target_reader(cpd_matrix)

if method_temp == 25:
    config_vals.update(rt_lock_t=11.0)
    report_vals.update(rt_lock_pass=0.025)
    report_vals.update(rt_lock_fail=0.05)
elif method_temp == 15:
    config_vals.update(rt_lock_t=15.5)
    report_vals.update(rt_lock_pass=0.0325)
    report_vals.update(rt_lock_fail=0.075)
elif method_temp == 7:
    config_vals.update(rt_lock_t=21.5)
    report_vals.update(rt_lock_pass=0.05)
    report_vals.update(rt_lock_fail=0.10)
elif method_temp == 'TBS':
    # TODO: not pass here
    pass

# Get date of experiment from andi_file attributes
expr_date = CDF(andi_file).attributes()['experiment_date_time_stamp']
formatted_date = time.strptime(expr_date[:-5], '%Y%m%d%H%M%S')
report_date = time.strftime('%Y%m%d_%H%M', formatted_date)

#==========================================================================
#    File processing
#==========================================================================
# Generate peak list
print ' -> Processing file', os.path.basename(andi_file)
# To change parameters:
#peak_list = get_peak_list(andi_file,method_temp,t=500)
peak_list = get_peak_list(andi_file, report_folder)

# If iQC file has failed horribly, there will be no peaks
if len(peak_list) == 0:
    print ' !! Apparently there are no peaks in this file. Please'
    print ' !! verify manually before submitting an error report.'
    qc_report_badfile(andi_file, report_folder, report_date)
    qc_db_update_badfile(db_file, andi_file, instr_name, report_date)
    raise RuntimeError

print ' -> Getting area and intensity values for compounds...'

all_sym = []
cpd_sym = []
for peak in peak_list:
    [lb_peak, ap_peak, rb_peak] = peak.get_pt_bounds()
    try:
        sym_val = float(lb_peak)/float(rb_peak)
    except:
        sym_val = 'NA'
    all_sym.append(sym_val)

for cpd in cpd_list:
    #print cpd.get_name()
    ion = cpd.get_ion()
    rt = cpd.get_rt()*60
    possibles = []
    
    for peak in peak_list:
        # Expand window for early eluting compounds
        if peak.get_rt() < 0.2 * max_time:
            if -10 < peak.get_rt()-rt < 10:
                possibles.append(peak)
        else:
            if -5 < peak.get_rt()-rt < 5:
                possibles.append(peak)
        
    # Get values
    try:
        probable = possibles[0]
        for peak in possibles:
            if peak.get_ion_area(ion) > probable.get_ion_area(ion):
                probable = peak
        
        # Determine whether the peak is symmetrical
        [lb_peak, ap_peak, rb_peak] = probable.get_pt_bounds()
        
        try:
            sym_val = float(lb_peak)/float(rb_peak)
        except:
            sym_val = 'NA'
        cpd_sym.append(sym_val)
        
        delta_val = fabs(float(probable.get_rt()) - rt)
        
        cpd.set_area(probable.get_ion_area(ion))
        cpd.set_intensity(probable.get_int_of_ion(ion))
        cpd.set_symmetry(sym_val)
        cpd.set_delta(delta_val)
        # Set rt to current rt (not that of the expected)
        cpd.set_rt(probable.get_rt()/60)
        
        # TODO: Calculate confidence of identity
        # Metric: 1) Closeness to expected/calculated RT
        #         2) Expected ions present (?)
        #         3) Ion ratios (?)
      
    except:
        cpd.set_area(0.0)
        cpd.set_intensity(0.0)
        cpd.set_delta('NA')
        cpd.set_symmetry('NA')
        cpd_sym.append('NA')
#==========================================================================

#
#     Append values to database for this file
#
print ' -> Updating qc database...'
qc_db_update(cpd_list, db_file, andi_file, instr_name, report_date)

#
#    Create images for report (test histories)
#
if method_temp in [7, 15, 25]:
    print ' -> Generating test history plots...'
    db_plotter(db_file, report_vals, os.path.basename(andi_file), \
        instr_name, report_folder)
elif method_temp == 'TBS':
    # TODO: Create db_plotter_tbs
#    db_plotter_tbs(db_file, report_vals, os.path.basename(andi_file), \
#        instr_name, report_folder)
    pass

#
#    Perform tests
#
if method_temp in [7, 15, 25]:
    #
    #    Generate the report for this run
    #
    # TODO: Check if this should be in report gen or here (here?)
    print ' -> Performing tests...'
    (test_report, watch_points, maint_reqs) = report_gen(\
        cpd_list, config_vals, report_vals, report_folder)
    symm_report = symm_test(all_sym, cpd_sym)
    params = []
    params.append(''.join(['GC-MS data file', ' & ', '\ldots'*5, ' & \\verb|',\
        os.path.basename(andi_file), '| \\\\']))
    params.append(''.join(['Instrument name', ' & ', '\ldots'*5, ' & ', 
        instr_name, '\\\\']))
    params.append(''.join(['CDF file timestamp', ' & ', '\ldots'*5, \
        ' & \\verb|', report_date, '| \\\\']))
    params.append(''.join(['Database file', ' & ', '\ldots'*5, ' & \\verb|',\
        os.path.basename(db_file), '| \\\\']))
    params.append(''.join(['Method temperature ramp', ' & ', \
        '\ldots'*5, ' & ', str(method_temp), ' \\\\']))
    params.append(''.join(['Method file', ' & ', '\ldots'*5, ' & \\verb|',\
        method_filename, '| \\\\']))
    params.append(''.join(['Retention Time Locking', ' & ', '\ldots'*5, \
        ' & ', str(config_vals['rt_lock_t']), ' min \\\\']))
    if options.target_filename != None:
        params.append(''.join(['Target compound file', ' & ', 
            '\ldots'*5, ' & ', options.target_filename, ' \\\\']))
    elif config_vals['target_filename'] != '':
        params.append(''.join(['Target compound file', ' & ',\
            '\ldots'*5, ' & ', config_vals['target_filename'], ' \\\\']))
    else:
        params.append(''.join(['Target compound file', ' & ', '\ldots'*5, \
            ' & ', 'None', ' \\\\']))
    
    params.append('\end{tabular*}')
    params.append('\\'*8)
    params.append('\\begin{tabular*}{0.9\\textwidth}{lcl@{ - }r}')
    params.append('Cutoff values & & Fail & Pass\\\\')
    
    for arg in sorted(report_labels.keys(), reverse=True):
        pass_val = report_vals[''.join([arg, '_pass'])]
        fail_val = report_vals[''.join([arg, '_fail'])]
        # Tidy up formatting
        if pass_val > 100:
            pass_val = int(pass_val) 
        else:
            pass_val = '%0.2f' % (pass_val) 
            
        if fail_val > 100:
            fail_val = int(fail_val)
        else:
            fail_val = '%0.2f' % (fail_val)
            
        params.append(''.join([
            # '\hspace{0.7cm}',
            report_labels[arg],
            ' & ', '\ldots'*5,' & ', \
            str(fail_val), ' & ', str(pass_val), \
            ' \\\\']))       
    
    print ' -> Creating test report...'
    qc_report(test_report, watch_points, maint_reqs, \
        symm_report, params, andi_file, report_folder, timecode=report_date)
    
elif method_temp == 'TBS':
#    (test_report, watch_points, maint_reqs) = report_gen_tbs(\
#        cpd_list, config_vals, report_vals, report_folder).
    qc_report_tbs(andi_file, report_folder, timecode=report_date)

print ' -> Done! Exiting.'

# EOF
