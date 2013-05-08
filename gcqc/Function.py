''' iQC.gcqc.Function.py
'''

# Requirements:
#     MA logo file (MA_logo_a3.png) in reports dir
#     Fancy header style file (fancyhdr.sty) in reports dir

import ConfigParser
import csv
import os
import re
import sqlite3
import string
import sys
import time

# TODO: Edit peak_top_ion_areas (don't forget def line) 
# TODO: Insert get_int_of_ion function in pyms.Peak.Class
sys.path.append('/x/PyMS/')
sys.path.append('/x/iQC/')

import matplotlib.pyplot as plt

from gcqc.Class import Compound
from numpy import arange, mean, median, std
from pycdf import CDF

from pyms.GCMS.IO.ANDI.Function import ANDI_reader
from pyms.GCMS.Function import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Noise.Analysis import window_analyzer
from pyms.Baseline.TopHat import tophat
from pyms.Peak.Class import Peak
from pyms.Peak.Function import peak_top_ion_areas
from pyms.Display.Class import Display

from pyms.Deconvolution.BillerBiemann.Function import BillerBiemann, \
    rel_threshold, num_ions_threshold

from pyms.Experiment.Class import Experiment
from pyms.Experiment.IO import store_expr
from pyms.Utils.IO import dump_object

from pylib.Function import config_reader, str2fname

# TODO: Fix this section. Not yet implemented.
def config_writer(filename=None):
    # NB: This will strip all comments, leaving only sections and options
    # Retain the 'default.cfg' file for reference 
    # Alternatively, write an ascii file straight up.
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    
    if filename == None:
        config_time = time.strftime('%Y%m%d_%H%M')
        filename = ''.join(['settings_', config_time, '.cfg'])
        
    with open(filename, 'wb') as configfile:
        config.write(configfile)

def file2matrix(filename):
    """"
    @summary:    Convert a .csv file to a matrix (list of lists)
    
    @param filename:    Filename (.csv) to convert
    @type filename:    StringType
    
    @return:    Data matrix
    @rtype:    ListType (List of lists)
    """
    
    with open(filename) as fp:
        reader = csv.reader(fp, delimiter=",",quotechar="\"")
        matrix = []
        for row in reader:
            newrow = []
            for each in row:
                try:
                    each = float(each)
                except:
                    pass
                newrow.append(each)
            matrix.append(newrow)
    
    return matrix

def target_reader(target_matrix):
    """
    @summary: Create a list of compounds from a matrix
    
    @param target_matrix:    Imported .csv file (file2matrix)
    @type target_matrix:    ListType (List of lists)
    
    @return:    List of compounds (as class Compound)
    @rtype:    ListType
    """
    # Find where Name, RI, MSTs and RT are
    nm_ind = target_matrix[0].index('Name')
    rt_ind = target_matrix[0].index('RT')
    #ri_ind = target_matrix[0].index('RI')
    i1_ind = target_matrix[0].index('MST1')
    target_list = []
    # Don't use the header row
    for line in target_matrix[1:]:
        # Compound(name, rt, ion); headers (Name, RI, MST1, MST2, RT)
        target = Compound(line[nm_ind],line[rt_ind],line[i1_ind])
        #target.set_ri(line[ri_ind])
        target_list.append(target)
    
#    for cpd in target_list:
#        print cpd.get_name(), " - ", cpd.get_rt()
    return target_list

def plot_tic(andi_data, fig_name):
    """
    @summary:    Save a plot of the tic
    
    @param andi_file:    Name of .cdf file
    @type andi_file:    pyms.GCMS.Class.GCMS_data object
    @param fig_name:    Name of output image file
    @type fig_name:    StringType
    """
    
    tic = andi_data.get_tic()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    intensity_list = tic.get_intensity_array()
    time_list = tic.get_time_list()
    mtime_list = []
    for time in time_list:
        mtime_list.append(time/60)
    
    plt.plot(mtime_list, intensity_list,color='r')
    plt.ylabel(('Intensity').title())
    plt.xlabel(('Time').title())
    
    fig.savefig(fig_name, dpi=300, transparent=True)

def get_peak_list(in_file, report_folder, \
    points=5, scans=2, n=3, t=3000, r=2):
    """
    @summary:    Create a peak list from a .cdf file
    
    @param in_file:    Input filename
    @type in_file:    StringType
    @param report_folder:    Location of reports
    @type report_folder:    StringType
    @param points:    Number of points used to determine a maxima
    @type points:    IntType
    @param scans:    Minimum number of points used to determine a peak
    @type scans:    IntType
    @param n:    Number of ions required to create a peak
    @type n:    IntType
    @param t:    Threshold value - below this value is considered noise
    @type t:    IntType
    @param r:    Relative intensity of ion as a percentage of  
                 highest intensity ion for the peak
    @type r:    IntType
    
    @return:    Peak list
    @rtype:    ListType
    """
    ##################################################################
    #
    # The first parameter 'points' is the window over which PyMS looks
    # for a local maximum. It should be approximately a third of the 
    # points per peak value and should be odd.
    #
    # The second parameter 'scans' determines how many adjacent points
    # a peak can be detected on and still be considered a single peak
    #
    # 'n' and 't' work together. The peak must have n number of Ions
    # with intensity above the value of 't' to remain in the list.
    #
    # r is a percentage. In an individual peak, if the intensity of any
    # ion is less than r% of the intensity of the highest intensity Ion 
    # in the peak, it is discarded from the peak
    #
    ###################################################################
    
    #==================================================================
    #    This needs to be checked - Do we need to analyse the entire
    #    datafile, or is it easier to only look for a peak in the 
    #    particular region of the compounds of interest?
    #==================================================================
    
    # Read in data file
    data = ANDI_reader(in_file)
    # Save the tic as an image
    base_name = os.path.splitext(os.path.basename(in_file))[0]
    tic_dir_path =  os.path.join( report_folder, 'images', base_name)
    
    if not os.path.isdir(tic_dir_path):
        os.makedirs(tic_dir_path)
    
    ticfig_name = os.path.join(tic_dir_path, \
        ''.join([base_name, '_tic.png']))
    plot_tic(data, ticfig_name)
    
    print ' -> Building intensity matrix...'
    # build integer intensity matrix
    im = build_intensity_matrix_i(data)
    # ignore TMS ions and use same mass range for all experiments
    # Performed at start to reduce further processing
    im.crop_mass(50,550)
    im.null_mass(73)
    im.null_mass(147)
    im.null_mass(207)
    
    # get the size of the intensity matrix
    n_scan, n_mz = im.get_size()
    
    print ' -> Smoothing data...'
    # smooth data
    for ii in range(n_mz):
        ic = im.get_ic_at_index(ii)
        ic1 = savitzky_golay(ic, points)
        ic_smooth = savitzky_golay(ic1, points)
        ic_base = tophat(ic_smooth, struct='1.5m')
        im.set_ic_at_index(ii, ic_base)
    
    # do peak detection on pre-trimmed data
    
    print ' -> Getting peak objects...'
    # get the list of Peak objects
    pl = BillerBiemann(im, points, scans)
    
    # trim by relative intensity
    apl = rel_threshold(pl, r)
    
    # trim by threshold
    peak_list = num_ions_threshold(apl, n, t)
    
    print '      Number of Peaks found:', len(peak_list)
    
    print ' -> Determining ion peak areas...'
    for peak in peak_list:
        # If the compound is ribitol or ribose, take the top ten ions
        # TODO: Verify if alanine also needs to be fiddled with
        # NB This value is for the 25 deg method
        # Old method - too complex when dealing with > 8 cpds
#        if method_temp == 25:
#            if 580 <= peak.get_rt() <= 630:
#                try:
#                    areas_dict, pt_bounds = peak_top_ion_areas( \
#                        im, peak, n_top_ions=10, pt_bounds=True)
#                except ZeroDivisionError:
#                    areas_dict = peak_top_ion_areas( \
#                        im, peak, n_top_ions=10, pt_bounds=False)
#                    pt_bounds = [0,0,0]
#            # Otherwise use the defaults
#            else:
#                try:
#                    areas_dict, pt_bounds = peak_top_ion_areas( \
#                        im, peak, pt_bounds=True)
#                except ZeroDivisionError:
#                    areas_dict = peak_top_ion_areas( \
#                        im, peak, n_top_ions=10, pt_bounds=False)
#                    pt_bounds = [0,0,0]
#        elif method_temp == 15:
#            if 800 <= peak.get_rt() <= 870:
#                [...]
#        elif method_temp == 7:
#            if 1020 <= peak.get_rt() <= 1120:
#                [...]

        # TODO: Add exception for citrate, isocitrate etc to get 
        #       top however many necessary
        try:
            areas_dict, pt_bounds = peak_top_ion_areas( \
                im, peak, n_top_ions=10, pt_bounds=True)
        except ZeroDivisionError:
            areas_dict = peak_top_ion_areas( \
                im, peak, n_top_ions=10, pt_bounds=False)
            pt_bounds = [0,0,0]
        
        peak.set_ion_areas(areas_dict)
        peak.set_pt_bounds(pt_bounds)
    
    # Return a peak list
    return peak_list

def get_cpd_data(cpd_name, col_name, input_matrix):
    '''
    @summary:    Retrieve list of area/intensity values for a compound
    
    @param cpd_name:    Name of compound
    @type cpd_name:    StringType
    @param col_name:    Area or Intensity of compound
    @type col_name:    StringType
    @param input_matrix:    Data to retrieve list from
    @type input_matrix:    ListType (List of lists)
    
    @return:    Compound name, area/intensity
    @rtype:    ListType
    '''
    # db file as [file_name, compound, area, intensity]
    #TODO: Update this; use index?
    if string.lower(col_name) == 'area':
        col_num = 2
    elif string.lower(col_name) == 'intensity':
        col_num = 3
    
    out_data = []
    for r in input_matrix:
        if string.lower(r[1]) == string.lower(cpd_name):
            out_data.append([r[0],r[col_num]])
    
    return out_data

def plot_history_old(input_list, list_type, plot_title, \
        filename, cutoffs, report_folder):
    '''
    @summary:    Plot barchart of QC history and save file
    
    @param input_list:    Date & (area/intensity/ratio) list
    @type input_list:    ListType (List of lists)
    @param list_type:    Area/intensity/ratio
    @type list_type:    StringType
    @param plot_title:    Ronseal
    @type plot_title:    StringType
    @param filename:    Name of .CDF file being processed
    @type filename:    StringType
    @param cutoffs:    Cut-off values for test
    @type cutoffs:    TupleType (2-tuple)
    '''
    n_list = len(input_list)
    dates = []
    vals = []
    for ln in input_list:
        dates.append(ln[0])
        vals.append(ln[1])
    
    # Replace any None or NA values with zero
    for val in vals:
        try:
            val = float(val)
        except TypeError:
            val = 0.0
    
    x_ind = arange(n_list)
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_xticklabels(dates,rotation='vertical')
    ax.set_xticks(x_ind)
#    ax.set_xlim(-0.75, n_list-0.25)  # see below
    ax.axhline(cutoffs[0],color='b')
    ax.axhline(cutoffs[1],color='b')
    ax.bar(x_ind,vals,color='r',align='center')
    # ax.set_xlim needs to be called after ax.bar to plot zero values
    ax.set_xlim(-0.75, n_list-0.25)
    #ax.xticks(range(len(dates))+1,)
    #ax.yticks(arange(0,max(ints)+1,1000)
    plt.xlabel('Experiment Date')
    plt.ylabel((list_type).title())
    plt.title((plot_title).title())
    
    # Bottom shifted down by 0.4 (empirical)
    fig.subplots_adjust(bottom=0.4)
#    plt.show()
    # Create a new directory for the images
    # Cut the '.cdf' from the filename (case-insensitive)
    #fname = re.sub('(.*[^\.])\.[Cc][Dd][Ff]', '\\1', filename)
    fname = os.path.splitext(filename)[0]
    dir_path = ''.join([report_folder, 'images/', fname,'/'])
    # Ensure the existence of 'images' directory
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
    
    img_name = ''.join([dir_path, str2fname(plot_title)])
    # Output image as file
    fig.savefig(img_name,dpi=300,transparent=True)

def plot_history(input_list, list_type, plot_title, \
        filename, cutoffs, report_folder):
    '''
    @summary:    Plot barchart of QC history and save file
    
    @param input_list:    File name & (intensity or ratio) list
    @type input_list:    ListType (List of lists)
    @param list_type:    Intensity or ratio values
    @type list_type:    StringType
    @param plot_title:    Ronseal
    @type plot_title:    StringType
    @param filename:    Name of .CDF file being processed
    @type filename:    StringType
    @param cutoffs:    Cut-off values for test
    @type cutoffs:    TupleType (2-tuple)
    '''
    n_list = len(input_list)
    fnames = []
    raw_vals = []
    for ln in input_list:
        fnames.append(ln[0])
        raw_vals.append(ln[1])
    
    # Replace any None or NA values with zero
    vals = []
    for val in raw_vals:
        try:
            vals.append(float(val))
        except TypeError:
            vals.append(0.0)
    
    x_ind = arange(n_list)
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_xticklabels(fnames,rotation='vertical')
    ax.set_xticks(x_ind)
#    ax.set_xlim(-0.75, n_list-0.25)  # see below
    ax.axhline(cutoffs[0],color='b')
    ax.axhline(cutoffs[1],color='b')
    ax.bar(x_ind,vals,color='r',align='center')
    # ax.set_xlim needs to be called after ax.bar to plot zero values
    ax.set_xlim(-0.75, n_list-0.25)
    #ax.xticks(range(len(dates))+1,)
    #ax.yticks(arange(0,max(ints)+1,1000)
    plt.xlabel('File name')
    plt.ylabel((list_type).title())
    plt.title((plot_title).title())
    
    # Bottom shifted down by 0.4 (empirical)
    fig.subplots_adjust(bottom=0.4)
#    plt.show()
    # Create a new directory for the images
    # Cut the '.cdf' from the filename (case-insensitive)
    #fname = re.sub('(.*[^\.])\.[Cc][Dd][Ff]', '\\1', filename)
    fname = os.path.splitext(filename)[0]
    dir_path = ''.join([report_folder, 'images/', fname,'/'])
    # Ensure the existence of 'images' directory
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
    
    img_name = ''.join([dir_path, str2fname(plot_title)])
    # Output image as file
    fig.savefig(img_name,dpi=300,transparent=True)


def qc_db_update_badfile(db_file, expr_name, instr_name, expr_date):
    '''
    @summary:    Adds database info for files that had no peaks
    
    @param db_file:    Database filename
    @type db_file:    StringType
    @param expr_name:    QC file name
    @type expr_name:    StringType
    @param instr_name:    Name of instrument file was run on
    @type instr_name:    StringType
    @param expr_date:    Date of experiment from .cdf file
    @type expr_date:    StringType
    '''
    fname = os.path.basename(expr_name)
    #curr_date = time.strftime('%Y%m%d_%H%M')
    # Get experiment date for file as '%Y%m%d_%H%M'
    m_file = CDF(expr_name).attributes()['external_file_ref_0']
    
    # Create db file if it doesn't exist
    if not os.path.isfile(db_file):
        f = open(db_file,'wb')
        f.close()
    
    # Open connection to database file
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    
    #
    #    File processing metadata
    #
    # Append to table
    try:
        cur.execute(''.join([ \
            'insert into qc_files (file_name, instr_name, '
            'method_file, expr_date) ', \
            'values (?,?,?,?)']), (fname, instr_name, m_file, expr_date)
        )
    except:
        # Create table if it doesn't exist
        cur.execute(''.join([ \
            'create table if not exists qc_files (', \
            'fkey integer primary key autoincrement unique not null, ', \
            'file_name varchar, ', \
            'method_file varchar, ', \
            'instr_name varchar, ', \
            'expr_date datetime)']))
        # Then append to table
        cur.execute(''.join([ \
            'insert into qc_files (file_name, instr_name, '
            'method_file, expr_date) ', \
            'values (?,?,?,?)']), (fname, instr_name, m_file, expr_date)
        )
    
    conn.commit()
    cur.close()
    conn.close()

def qc_db_update(cpd_list, db_file, expr_name, \
    instr_name, expr_date, report_labels, test_results):
    '''
    @summary:    Creates a database table of intensities and areas, and
                 records what date files were processed on
    @param cpd_list:    List of compounds
    @type cpd_list:    ListType
    @param db_file:    Database filename
    @type db_file:    StringType
    @param expr_name:    QC file name
    @type expr_name:    StringType
    @param instr_name:    Name of instrument file was run on
    @type instr_name:    StringType
    @param expr_date:    Date of experiment from .cdf file
    @type expr_date:    StringType
    @param report_labels:    Names of tests used in report
    @type report_labels:    DictionaryType
    @param test_results:    Numeric and Pass/Watch/Fail results for all tests
    @type test_results:    ListType
    '''
    fname = os.path.basename(expr_name)
    #curr_date = time.strftime('%Y%m%d_%H%M')
    # Get experiment date for file as '%Y%m%d_%H%M'
    m_file = CDF(expr_name).attributes()['external_file_ref_0']
    
    # Create db file if it doesn't exist
    if not os.path.isfile(db_file):
        f = open(db_file,'wb')
        f.close()
    
    # Open connection to database file
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    
    #
    #    File processing metadata
    #
    # Append to table
    try:
        cur.execute(''.join([ \
            'insert into qc_files (file_name, instr_name, '
            'method_file, expr_date) ', \
            'values (?,?,?,?)']), (fname, instr_name, m_file, expr_date)
        )
    except:
        # Create table if it doesn't exist
        cur.execute(''.join([ \
            'create table if not exists qc_files (', \
            'fkey integer primary key autoincrement unique not null, ', \
            'file_name varchar, ', \
            'method_file varchar, ', \
            'instr_name varchar, ', \
            'expr_date datetime)']))
        # Then append to table
        cur.execute(''.join([ \
            'insert into qc_files (file_name, instr_name, '
            'method_file, expr_date) ', \
            'values (?,?,?,?)']), (fname, instr_name, m_file, expr_date)
        )
    
    #
    #    Compound areas and intensities
    #
    try:
        for cpd in cpd_list:
            cur.execute(''.join([ \
                'insert into qc_data (', \
                'file_name, expr_date, compound, area, ', \
                'intensity, RT, RTdelta, symmetry', \
                ') values (?, ?, ?, ?, ?, ?, ?, ?)']), \
                (fname, expr_date, cpd.get_name(), cpd.get_area(), \
                cpd.get_intensity(), cpd.get_rt(), \
                cpd.get_delta(), cpd.get_symmetry())
            )
    except:
        # Create if table doesn't exist
        cur.execute(''.join([ \
            'create table if not exists qc_data (', \
            'file_name varchar, expr_date varchar, ' \
            'compound varchar, area real, intensity real, ' \
            'RT real, RTdelta real, symmetry real)'])
        )
        
        for cpd in cpd_list:
            cur.execute(''.join([
                'insert into qc_data ('\
                'file_name, expr_date, compound, area, ' \
                'intensity, RT, RTdelta, symmetry'
                ') values (?, ?, ?, ?, ?, ?, ?, ?)']), \
                (fname, expr_date, cpd.get_name(), cpd.get_area(), \
                cpd.get_intensity(), cpd.get_rt(), \
                cpd.get_delta(), cpd.get_symmetry()) \
            )
    #
    #    Test results
    #
    try:
        cur.execute(''.join([\
            'insert into qc_results (file_name, ', \
                'test1r1, test1r2, ', \
               'test2r1, test2r2, ', \
                'test3r1, test3r2, ', \
                'test4r1, test4r2, ', \
                'test5r1, test5r2, ', \
                'test6r1, test6r2, ', \
                'test7r1, test7r2, ', \
                'test8r1, test8r2, ', \
                'test9r1, test9r2', \
           ') values (', \
                '?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?', \
            ')']), \
            (fname, test_results[0], test_results[1], test_results[2], \
                test_results[3], test_results[4], test_results[5], \
                test_results[6], test_results[7], test_results[8], \
                test_results[9], test_results[10], test_results[11], \
                test_results[12], test_results[13], test_results[14], \
                test_results[15], test_results[16], test_results[17] \
            )
        )
    
    # If table doesn't exist, nor does the test ID table
    except sqlite3.OperationalError:
        # Create testID table
        cur.execute(''.join([ \
            'create table if not exists qc_testID (test_number varchar, ', \
            'test_label varchar)'])
        )
        # Populate testID table
        ordered_keys = ['rt_lock', 'rbtl_i', 'pr_r', 'rbose_i', 'gr_r', \
            'ar_r', 'ala_i', 'av_r', 'mlt_i']
        for ii in range(len(ordered_keys)):
            cur.execute(''.join([ \
                'insert into qc_testID (test_number, test_label', \
                ') values (?, ?)']), \
                (''.join(['test', str(ii+1)]), \
                    report_labels[ordered_keys[ii]] \
                )
            )
        # Create qc_results table
        cur.execute(''.join([ \
            'create table if not exists qc_results (file_name varchar,', \
            'test1r1 real, test1r2 varchar, ', \
            'test2r1 real, test2r2 varchar, ', \
            'test3r1 real, test3r2 varchar, ', \
            'test4r1 real, test4r2 varchar, ', \
            'test5r1 real, test5r2 varchar, ', \
            'test6r1 real, test6r2 varchar, ', \
            'test7r1 real, test7r2 varchar, ', \
            'test8r1 real, test8r2 varchar, ', \
            'test9r1 real, test9r2 varchar)'])
        )
        # Populate qc_results table
        cur.execute(''.join([ \
            'insert into qc_results (file_name, ', \
                'test1r1, test1r2, ', \
                'test2r1, test2r2, ', \
                'test3r1, test3r2, ', \
                'test4r1, test4r2, ', \
                'test5r1, test5r2, ', \
                'test6r1, test6r2, ', \
                'test7r1, test7r2, ', \
                'test8r1, test8r2, ', \
                'test9r1, test9r2', \
            ') values (', \
                '?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?', \
            ')']), \
            (fname, test_results[0], test_results[1], test_results[2], \
                test_results[3], test_results[4], test_results[5], \
                test_results[6], test_results[7], test_results[8], \
                test_results[9], test_results[10], test_results[11], \
                test_results[12], test_results[13], test_results[14], \
                test_results[15], test_results[16], test_results[17] \
            )
        )
    

    
    conn.commit()
    cur.close()
    conn.close()

def db_plotter_tbs(db_file, report_vals, filename, instr_name, report_folder):
    '''
    @summary:    Plot histories of each test from a database file
    @param db_file:    Database filename
    @type db_file:    StringType
    @param filename:    Name of .CDF file being processed
    @type filename:    StringType
    @param report_vals:    Cut-off values for tests
    @type report_vals:    DictionaryType 
    '''
    # TODO: Not pass. Placeholder function.
    pass

# Deprecated: Now have a new table for qc_results. Reduces process time.
def db_plotter_old(db_file, report_vals, filename, instr_name, report_folder):
    '''
    @summary:    Plot histories of each test from a database file
    @param db_file:    Database filename
    @type db_file:    StringType
    @param filename:    Name of .CDF file being processed
    @type filename:    StringType
    @param report_vals:    Cut-off values for tests
    @type report_vals:    DictionaryType 
    '''
    
    # Open connection to database file
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    
    # Get last (8) results from db file 
    raw_fetch = cur.execute(''.join([ \
        'select file_name, compound, area, intensity from qc_data where ', \
        'file_name in (', \
            'select file_name from qc_files where instr_name = "', \
            instr_name, \
            '" order by expr_date desc limit 8', \
        ') order by expr_date asc'])).fetchall()
    
    # Tidy up the output
    clean_fetch = []
    for ii in range(len(raw_fetch)):
        tmp = [
            str(raw_fetch[ii][0]),              # File name
            str(raw_fetch[ii][1]),              # Compound name
            raw_fetch[ii][2],                   # Area
            raw_fetch[ii][3]                    # Intensity
        ]
        clean_fetch.append(tmp)
    
    cur.close()
    conn.close()
    
    # Get lists of vals
    # Ribitol intensity test
    rbtl_ilist = get_cpd_data('ribitol', 'intensity', clean_fetch)
    rbtl_cutoffs = (report_vals['rbtl_i_pass'], report_vals['rbtl_i_fail'])
    # Generate plot
    plot_history(rbtl_ilist, 'intensity','ribitol intensity', filename, \
        rbtl_cutoffs, report_folder)
    
    # Putrescine/Ribitol ratio test
    putr_ilist = get_cpd_data('putrescine 4tms', 'intensity', clean_fetch)
    pr_ratio = []
    for ii in range(len(rbtl_ilist)):
        try:
            pr_ratio.append(\
                [rbtl_ilist[ii][0],putr_ilist[ii][1]/rbtl_ilist[ii][1]]\
            )
        except ZeroDivisionError:
            pr_ratio.append([rbtl_ilist[ii][0], 0])
    
    pr_cutoffs = (report_vals['pr_r_pass'], report_vals['pr_r_fail'])
    # Generate plot
    plot_history(pr_ratio, 'ratio','Putrescine/Ribitol ratio (i/i)', \
        filename, pr_cutoffs, report_folder)
    
    ## Ribose intensity test
    rbose_ilist = get_cpd_data('ribose', 'intensity', clean_fetch)
    rbose_cutoffs = (report_vals['rbose_i_pass'], report_vals['rbose_i_fail'])
    # Generate plot
    plot_history(rbose_ilist, 'intensity', 'ribose intensity', filename, \
        rbose_cutoffs, report_folder)
    
    ## Glutamate/Ribose ratio test
    glu_ilist = get_cpd_data('glutamate 3tms', 'intensity', clean_fetch)
    grs_ratio = []
    for ii in range(len(glu_ilist)):
        try:
            grs_ratio.append(\
                [glu_ilist[ii][0],glu_ilist[ii][1]/rbose_ilist[ii][1]]\
            )
        except ZeroDivisionError:
            grs_ratio.append([glu_ilist[ii][0], 0])
    
    grs_cutoffs = (report_vals['gr_r_pass'], report_vals['gr_r_fail'])
    # Generate plot
    plot_history(grs_ratio, 'ratio','Glutamate/Ribose ratio (i/i)', \
        filename, grs_cutoffs, report_folder)
    
    ## Asparagine/Ribose ratio
    asn_ilist = get_cpd_data('asparagine 3tms', 'intensity', clean_fetch)
    ars_ratio = []
    for ii in range(len(asn_ilist)):
        try:
            ars_ratio.append(\
                [asn_ilist[ii][0],asn_ilist[ii][1]/rbose_ilist[ii][1]]\
            )
        except ZeroDivisionError:
            ars_ratio.append([asn_ilist[ii][0], 0])
    
    ars_cutoffs = (report_vals['ar_r_pass'], report_vals['ar_r_fail'])
    # Generate plot
    plot_history(ars_ratio, 'ratio','Asparagine/Ribose ratio (i/i)', \
        filename, ars_cutoffs, report_folder)
    
    ## Column condition - start
    ala_ilist = get_cpd_data('alanine 2tms', 'intensity', clean_fetch)
    ala_cutoffs = (report_vals['ala_i_pass'], report_vals['ala_i_fail'])
    # Generate plot
    plot_history(ala_ilist, 'intensity', 'alanine intensity', filename, \
        ala_cutoffs, report_folder)
    
    val_ilist = get_cpd_data('valine 2tms', 'intensity', clean_fetch)
#    val_alist = get_cpd_data('valine', 'area', clean_fetch)
    av_ratio = []
    for ii in range(len(ala_ilist)):
        try:
            av_ratio.append(\
                [ala_ilist[ii][0],ala_ilist[ii][1]/val_ilist[ii][1]]\
            )
        except ZeroDivisionError:
            av_ratio.append([ala_ilist[ii][0], 0])
    
    av_cutoffs = (report_vals['av_r_pass'], report_vals['av_r_fail'])
    # Generate plot
    plot_history(av_ratio, 'ratio','Alanine/Valine ratio (i/i)', \
        filename, av_cutoffs, report_folder)
    
    ##Column condition - end
    mlt_ilist = get_cpd_data('maltotriose mx1', 'intensity', clean_fetch)
    mlt_cutoffs = (report_vals['mlt_i_pass'], report_vals['mlt_i_fail'])
    # Generate plot
    plot_history(mlt_ilist, 'intensity','maltotriose intensity', filename, \
        mlt_cutoffs, report_folder)

def db_plotter(db_file, report_vals, filename, method_name, instr_name, \
    report_folder):
    '''
    @summary:    Plot histories of each test from a database file
    @param db_file:    Database filename
    @type db_file:    StringType
    @param report_vals:    Cut-off values for tests
    @type report_vals:    DictionaryType
    @param filename:    Name of .CDF file being processed
    @type filename:    StringType
    @param method_name:    Name of method file used for file
    @type method_name:    StringType
    @param instr_name:    Name of the instrument used for file
    @type instr_name:    StringType
    @param report_folder:    Name of the report folder to store plots
    @type report_folder:    StringType
    '''
    
    # Open connection to database file
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    
    # Get last (8) results from db file 
    raw_fetch = cur.execute(''.join([ \
    ### Old method, prior to results table being in database
    #    'select file_name, compound, area, intensity from qc_data where ', \
    #    'file_name in (', \
    #        'select file_name from qc_files where instr_name = "', \
    #        instr_name, \
    #        '" order by expr_date desc limit 8', \
    #    ') order by expr_date asc'])).fetchall()
        'select * from qc_results where file_name in (', \
            'select file_name from qc_files where instr_name = "', \
            instr_name, '" and method_file = "', method_name, \
            '" order by expr_date desc limit 8', \
        ')'])).fetchall()
    
    cur.close()
    conn.close()
    
    # This now has 8 entries, the most recent being first. Reverse this:
    raw_fetch.reverse()
    
    # Tidy up the output
    clean_fetch = []
    for ii in range(len(raw_fetch)):
        tmp = [
            str(raw_fetch[ii][0]),              # File name
            raw_fetch[ii][3],                   # Ribitol intensity
            raw_fetch[ii][5],                   # Putrescine/Ribitol ratio
            raw_fetch[ii][7],                   # Ribose intensity
            raw_fetch[ii][9],                   # Glutamate/Ribose ratio
            raw_fetch[ii][11],                  # Asparagine/Ribose ratio
            raw_fetch[ii][13],                  # Alanine intensity
            raw_fetch[ii][15],                  # Alanine/Valine ratio
            raw_fetch[ii][17]                   # Maltotriose intensity
        ]
        clean_fetch.append(tmp)
    
    # Collect cutoff values
    rbtl_cutoffs = (report_vals['rbtl_i_pass'], report_vals['rbtl_i_fail'])
    pr_cutoffs = (report_vals['pr_r_pass'], report_vals['pr_r_fail'])
    rbose_cutoffs = (report_vals['rbose_i_pass'], report_vals['rbose_i_fail'])
    grs_cutoffs = (report_vals['gr_r_pass'], report_vals['gr_r_fail'])
    ars_cutoffs = (report_vals['ar_r_pass'], report_vals['ar_r_fail'])
    ala_cutoffs = (report_vals['ala_i_pass'], report_vals['ala_i_fail'])
    av_cutoffs = (report_vals['av_r_pass'], report_vals['av_r_fail'])
    mlt_cutoffs = (report_vals['mlt_i_pass'], report_vals['mlt_i_fail'])
    
    # Prepare result lists
    rbtl_ilist = []
    pr_ratio = []
    rbose_ilist = []
    grs_ratio = []
    ars_ratio = []
    ala_ilist = []
    av_ratio = []
    mlt_ilist = []
    
    # Plotting lists require file name, (intensity/ratio) pairs
    for rw in clean_fetch:
        rbtl_ilist.append([rw[0], rw[1]])
        pr_ratio.append([rw[0], rw[2]])
        rbose_ilist.append([rw[0], rw[3]])
        grs_ratio.append([rw[0], rw[4]])
        ars_ratio.append([rw[0], rw[5]])
        ala_ilist.append([rw[0], rw[6]])
        av_ratio.append([rw[0], rw[7]])
        mlt_ilist.append([rw[0], rw[8]])
    
    # Create plots
    plot_history(rbtl_ilist, 'intensity','ribitol intensity', filename, \
        rbtl_cutoffs, report_folder)
    plot_history(pr_ratio, 'ratio','Putrescine/Ribitol ratio (i/i)', \
        filename, pr_cutoffs, report_folder)
    plot_history(rbose_ilist, 'intensity', 'ribose intensity', filename, \
        rbose_cutoffs, report_folder)
    plot_history(grs_ratio, 'ratio','Glutamate/Ribose ratio (i/i)', \
        filename, grs_cutoffs, report_folder)
    plot_history(ars_ratio, 'ratio','Asparagine/Ribose ratio (i/i)', \
        filename, ars_cutoffs, report_folder)
    plot_history(ala_ilist, 'intensity', 'alanine intensity', filename, \
        ala_cutoffs, report_folder)
    plot_history(av_ratio, 'ratio','Alanine/Valine ratio (i/i)', \
        filename, av_cutoffs, report_folder)
    plot_history(mlt_ilist, 'intensity','maltotriose intensity', filename, \
        mlt_cutoffs, report_folder)

def report_gen_tbs(cpd_list, config_vals, report_vals, report_folder):
    '''
    @summary:    Creates test results, watch items and maintenance items
                 for TBS method
    
    @param cpd_list:    List of compounds to check
    @type cpd_list:    ListType (Compound)
    @param report_vals:    Cut-off values for tests
    @type report_vals:    DictionaryType 
    
    @return:    test_report, watch_points, maint_reqs
    @rtype:    TupleType
    '''
    pass

def report_gen(cpd_list, config_vals, report_vals, report_folder):
    '''
    @summary:    Creates test results, watch items and maintenance items
    
    @param cpd_list:    List of compounds to check
    @type cpd_list:    ListType (Compound)
    @param report_vals:    Cut-off values for tests
    @type report_vals:    DictionaryType 
    
    @return:    test_report, watch_points, maint_reqs
    @rtype:    TupleType
    '''
    # Create individual compounds from cpd_list for testing
    #    Mannitol, ribitol and ribose compound abbreviations altered
    #    so there is no conflict in namespace
    for cpd in cpd_list:
        if cpd.get_name() == 'Alanine 2TMS':
            ala = cpd
        elif cpd.get_name() == 'Asparagine 3TMS':
            asn = cpd
        elif cpd.get_name() == 'Glutamate 3TMS':
            glu = cpd
        elif cpd.get_name() == 'Maltotriose MX1':
            mlt = cpd
        elif cpd.get_name() == 'Mannitol':
            mann_test = cpd
        elif cpd.get_name() == 'Putrescine 4TMS':
            putr = cpd
        elif cpd.get_name() == 'Ribitol':
            rbtl_test = cpd
        elif cpd.get_name() == 'Ribose':
            rbose_test = cpd
        elif cpd.get_name() == 'Valine 2TMS':
            val = cpd
    
    # Prepare report containers
    test_results = []
    test_report = []
    watch_points = []
    maint_reqs = []
    
    #=========================================================================
    #    Perform the tests
    #=========================================================================
    #
    #    Retention Time Lock
    #
    mann_t = mann_test.get_rt()
    lock_d = abs(config_vals['rt_lock_t'] - mann_t)
    rtl_pass = report_vals['rt_lock_pass']
    rtl_fail = report_vals['rt_lock_fail']
    
    if lock_d < rtl_pass:
        rt_lock_test = 'PASS'
    elif rtl_pass <= lock_d <= rtl_fail:
        rt_lock_test = 'WATCH'
        watch_points.append('Retention Time Locking suboptimal')
    else:
        rt_lock_test = '\\textbf{FAIL}'
        maint_reqs.append('Re-lock Mannitol')
    
    # Add result to report
    test_results.append(lock_d)
    if rt_lock_test in ['PASS', 'WATCH']:
        test_results.append(rt_lock_test)
    else:
        test_results.append('FAIL')
    
    test_report.append(('Retention Time Lock', rt_lock_test))
    
    #
    #    Inlet conditions
    #
    # Ribitol intensity test
    rbtl_i = rbtl_test.get_intensity()
    
    if rbtl_i > report_vals['rbtl_i_pass']:
        rbtl_i_test = 'PASS'
    elif report_vals['rbtl_i_fail'] <= rbtl_i <= report_vals['rbtl_i_pass']:
        rbtl_i_test = 'WATCH'
        watch_points.append('Possible air leak or MS tuning suboptimal')
    else:
        rbtl_i_test = '\\textbf{FAIL}'
        maint_reqs.append('Check for air leaks. Retune MS')
    
    # Add result to report
    test_results.append(rbtl_i)
    if rbtl_i_test in ['PASS', 'WATCH']:
        test_results.append(rbtl_i_test)
    else:
        test_results.append('FAIL')
    
    test_report.append(('Ribitol intensity test', rbtl_i_test))
    
    # Putrescine/Ribitol ratio test
    putr_i = putr.get_intensity()
    
    try:
        pr_val = putr_i/rbtl_i
        if pr_val > report_vals['pr_r_pass']:
            pr_r = 'PASS'
        elif report_vals['pr_r_fail'] <= pr_val <= report_vals['pr_r_pass']:
            pr_r = 'WATCH'
            watch_points.append('Inlet liner may need replacing soon')
        else:
            pr_r = '\\textbf{FAIL}'
            maint_reqs.append('Replace inlet liner')
    except ZeroDivisionError:
        pr_val = 0.0
        pr_r = '\\textbf{FAIL}'
        maint_reqs.append('Replace inlet liner')
    
    test_results.append(pr_val)
    if pr_r in ['PASS', 'WATCH']:
        test_results.append(pr_r)
    else:
        test_results.append('FAIL')
    
    test_report.append(('Putrescine/Ribitol ratio test', pr_r))
    
    # Ribose intensity test
    rbose_i = rbose_test.get_intensity()
    
    if rbose_i > report_vals['rbose_i_pass']:
        rbose_i_test = 'PASS'
    elif report_vals['rbose_i_fail'] <= rbose_i <= report_vals['rbose_i_pass']:
        rbose_i_test = 'WATCH'
        watch_points.append('Syringe and inlet liner may need\
             maintenance soon (1)')
    else:
        rbose_i_test = '\\textbf{FAIL}'
        maint_reqs.append(\
            'Clean syringe, replace inlet liner and bake oven (1)')
    
    test_results.append(rbose_i)
    if rbose_i_test in ['PASS', 'WATCH']:
        test_results.append(rbose_i_test)
    else:
        test_results.append('FAIL')
    
    test_report.append(('Ribose intensity test', rbose_i_test))
    
    # Glutamate/Ribose ratio test
    glu_i = glu.get_intensity()
    
    
    try:
        gr_val = glu_i/rbose_i
        if gr_val > report_vals['gr_r_pass']:
            gr_r = 'PASS'
        elif report_vals['gr_r_fail'] <= gr_val <= report_vals['gr_r_pass']:
            gr_r = 'WATCH'
            watch_points.append(\
                'Syringe and inlet liner may need maintenance soon (2)')
        else:
            gr_val = 0.0
            gr_r = '\\textbf{FAIL}'
            maint_reqs.append(\
                'Clean syringe, replace inlet liner and bake oven (2)')
    except ZeroDivisionError:
        gr_r = '\\textbf{FAIL}'
        maint_reqs.append(\
            'Clean syringe, replace inlet liner and bake oven (2)')
    
    test_results.append(gr_val)
    if gr_r in ['PASS', 'WATCH']:
        test_results.append(gr_r)
    else:
        test_results.append('FAIL')
    
    test_report.append(('Glutamate/Ribose ratio test', gr_r))
    
    # Asparagine/Ribose ratio
    asn_i = asn.get_intensity()
    
    try:
        ar_val = asn_i/rbose_i
        if ar_val > report_vals['ar_r_pass']:
            ar_r = 'PASS'
        elif report_vals['ar_r_fail'] <= ar_val <= report_vals['ar_r_pass']:
            ar_r = 'WATCH'
            watch_points.append(\
                'Inlet liner and seal may need replacing soon')
        else:
            ar_r = '\\textbf{FAIL}'
            maint_reqs.append(\
                'Clean inlet port, replace liner and seal, bake oven')
    except ZeroDivisionError:
        ar_val = 0.0
        ar_r = '\\textbf{FAIL}'
        maint_reqs.append(\
            'Clean inlet port, replace liner and seal, bake oven')
    
    test_results.append(ar_val)
    if ar_r in ['PASS', 'WATCH']:
        test_results.append(ar_r)
    else:
        test_results.append('FAIL')
    
    test_report.append(('Asparagine/Ribose ratio', ar_r))
    
    #
    #    Column Conditions
    #
    # Alanine intensity test - Column condition - start (1)
    ala_i = ala.get_intensity()
    
    if ala_i > report_vals['ala_i_pass']:
        ala_i_test = 'PASS'
    elif report_vals['ala_i_fail'] <= ala_i <= report_vals['ala_i_pass']:
        ala_i_test = 'WATCH'
        watch_points.append('Start of column may need cutting soon (1)')
    else:
        ala_i_test = '\\textbf{FAIL}'
        maint_reqs.append('Cut start of column (1)')
    
    test_results.append(ala_i)
    if ala_i_test in ['PASS', 'WATCH']:
        test_results.append(ala_i_test)
    else:
        test_results.append('FAIL')
    
    test_report.append(('Alanine intensity test', ala_i_test))
    
    # Column condition - start (2)
    val_i = val.get_intensity()
    
    try:
        av_val = ala_i/val_i
        if av_val > report_vals['av_r_pass']:
            av_r = 'PASS'
        elif report_vals['av_r_fail'] <= av_val <= report_vals['av_r_pass']:
            av_r = 'WATCH'
            watch_points.append('Start of column may need cutting soon (2)')
        else:
            av_r = '\\textbf{FAIL}'
            maint_reqs.append('Cut start of column (2)')
    except ZeroDivisionError:
        av_val = 0.0
        av_r = '\\textbf{FAIL}'
        maint_reqs.append('Cut start of column (2)')
    
    test_results.append(av_val)
    if av_r in ['PASS', 'WATCH']:
        test_results.append(av_r)
    else:
        test_results.append('FAIL')
    
    test_report.append(('Alanine/Valine ratio test', av_r))
    
    #Column condition - end
    mlt_i = mlt.get_intensity()
    
    if mlt_i > report_vals['mlt_i_pass']:
        mlt_i_test = 'PASS'
    elif report_vals['mlt_i_fail'] <= mlt_i <= report_vals['mlt_i_pass']:
        mlt_i_test = 'WATCH'
        watch_points.append('End of column may need cutting soon')
    else:
        mlt_i_test = '\\textbf{FAIL}'
        maint_reqs.append('Cut end of column')
    
    test_results.append(mlt_i)
    if mlt_i_test in ['PASS', 'WATCH']:
        test_results.append(mlt_i_test)
    else:
        test_results.append('FAIL')
    
    test_report.append(('Maltotriose test', mlt_i_test))
    
    return (test_results, test_report, watch_points, maint_reqs)

def symm_test(all_sym, cpd_sym):
    '''
    @summary:    Compare peak symmetry for all peaks and for the compound list
    
    @param all_sym:    Symmetry ratios for all peaks
    @type all_sym:    ListType
    @param cpd_sym:    Symmetry ratios for compound list items
    @type cpd_sym:    ListType
    '''
    # Remove NA from symmetry vals
    sub_all_sym = []
    sub_cpd_sym = []
    for val in all_sym:
        if type(val) == float:
            sub_all_sym.append(val)
    
    for val in cpd_sym:
        if type(val) == float:
            sub_cpd_sym.append(val)
    
    # All compounds
    a_mean = mean(sub_all_sym)
    a_median = median(sub_all_sym)
    a_std = std(sub_all_sym)
    c_mean = mean(sub_cpd_sym)
    c_median = median(sub_cpd_sym)
    c_std = std(sub_cpd_sym)
    symm = [['', 'Mean', 'Median', 'Std Dev', 'CV']]
    symm.append(['All peaks', a_mean, a_median, a_std, (a_std/a_mean)])
    symm.append(['QC5 components', c_mean, c_median, c_std, (c_std/c_mean)])
    
    return symm

def qc_report_badfile(andi_file, report_folder, \
        timecode=time.strftime('%Y%m%d_%H%M')\
    ):
    '''
    @summary:    Creates a .tex report for files with no peaks
    
    @param andi_file:    Name of file being processed
    @type andi_file:    StringType
    @param report_folder:    Where reports will be stored
    @type report_folder:    StringType
    @param timecode:    Time stamp 
    @type timecode:    StringType
    '''
    # Use timecode provided
    #report_date = time.strftime('%Y%m%d_%H%M')
    report_date = timecode
    andi_base = os.path.basename(andi_file)[:-4]
    # Add in document headers
    tex_doc = [ \
        ''.join(['% report_time = ',report_date]), \
        '\documentclass[10pt,a4paper]{report}', \
        '%%% Put title and MA logo in document header', \
        '\usepackage{graphicx}', \
        # The fancy header package is required. If not installed, report will
        # not be generated. This is part of the texlive-full package on most
        # linux distros. If not available, or not installed properly, the 
        # package can be downloaded from the ctan site at:  
        # http://www.ctan.org/tex-archive/macros/latex/contrib/fancyhdr/
        '\usepackage{fancyhdr}', \
        '\usepackage[%', \
        '    top=1.5cm,bottom=3.5cm,left=1.5cm,right=1.5cm,%', \
        '    paper=a4paper,includeheadfoot]{geometry}', \
        '% Allow multiple columns', \
        '\usepackage{multicol}', \
        '\headheight=72pt', \
        '\headsep=10pt', \
        '\pagestyle{fancy}', \
        '\\renewcommand{\headrulewidth}{0pt}', \
        '\\fancyhead{}', \
        '\\fancyhead[L]{\\textbf{\Large{GC-MS Quality Control Report}}\\\\}', \
        # This is where pdflatex is looking for the logo (used by the
        # report_mailer function)
        ''.join(['\\fancyhead[R]{\includegraphics[height=48pt]', \
            '{images/MA_logo_a3.png}%']), \
        '    \\\\ \\today}', \
        '% Compact lists', \
        '\\newenvironment{itemize*}%', \
        '    {\\begin{itemize}%', \
        '        \setlength{\itemsep}{2pt}%', \
        '        \setlength{\parskip}{0pt}}%', \
        '    {\end{itemize}}', \
        '% Compact headings', \
        '\setlength{\parskip}{0pt}', \
        '\setlength{\headsep}{0pt}', \
        '\setlength{\partopsep}{0pt}', \
        '\setlength{\\topskip}{0pt}', \
        '\usepackage[compact]{titlesec}', \
        '    \\titlespacing{\section}{2pt}{*1}{*1}', \
        '    \\titlespacing{\subsection}{2pt}{*0}{*0}', \
        '', \
        '\\begin{document}', \
        '', \
        '\section*{Processing Error}', \
        ''.join([ \
            'The file \\verb|', os.path.basename(andi_file), \
            '| appears to have some rather serious issues.', \
        ]), \
        '\\'*4, \
        'Please ensure everything is as it should be prior to submitting', \
        'an error report.', \
        '\\'*4, \
        '\\begin{figure}[h]', \
        ''.join([' '*4, '\centering']), \
        ''.join([' '*8, '\includegraphics[width=0.9\\textwidth]{%']), \
        ''.join([' '*12, \
            'images/', andi_base, '/', andi_base, '_tic}']), \
        '\end{figure}', \
        '\\\\', '\end{document}']
    
    # Use .cdf file name in report filename
    andi_file_noext = os.path.splitext(os.path.basename(andi_file))[0]
    fname_tex = ''.join(\
        [report_folder, 'GCQC_Report_', andi_file_noext,'.tex']\
    )
    
    with open(fname_tex, 'wb') as f:
        for line in tex_doc:
            f.write(line)
            f.write('\n')
    
    f.close()

def qc_report_tbs(andi_file, report_folder, \
    timecode=time.strftime('%Y%m%d_%H%M')):
    '''
    @summary:    Creates a .tex report for TBS files
    
    @param andi_file:    Name of file being processed
    @type andi_file:    StringType
    @param report_folder:    Where reports will be stored
    @type report_folder:    StringType
    @param timecode:    Time stamp 
    @type timecode:    StringType
    
    '''
    report_date = timecode
    andi_base = os.path.basename(andi_file)[:-4]
    # Add in document headers
    tex_doc = [ \
        ''.join(['% report_time = ',report_date]), \
        '\documentclass[10pt,a4paper]{report}', \
        '%%% Put title and MA logo in document header', \
        '\usepackage{graphicx}', \
        # The fancy header package is required. If not installed, report will
        # not be generated. This is part of the texlive-full package on most
        # linux distros. If not available, or not installed properly, the 
        # package can be downloaded from the ctan site at:  
        # http://www.ctan.org/tex-archive/macros/latex/contrib/fancyhdr/
        '\usepackage{fancyhdr}', \
        '\usepackage[%', \
        '    top=1.5cm,bottom=3.5cm,left=1.5cm,right=1.5cm,%', \
        '    paper=a4paper,includeheadfoot]{geometry}', \
        '% Allow multiple columns', \
        '\usepackage{multicol}', \
        '\headheight=72pt', \
        '\headsep=10pt', \
        '\pagestyle{fancy}', \
        '\\renewcommand{\headrulewidth}{0pt}', \
        '\\fancyhead{}', \
        '\\fancyhead[L]{\\textbf{\Large{GC-MS Quality Control Report}}\\\\}', \
        # This is where pdflatex is looking for the logo (used by the
        # report_mailer function)
        ''.join(['\\fancyhead[R]{\includegraphics[height=48pt]', \
            '{images/MA_logo_a3.png}%']), \
        '    \\\\ \\today}', \
        '% Compact lists', \
        '\\newenvironment{itemize*}%', \
        '    {\\begin{itemize}%', \
        '        \setlength{\itemsep}{2pt}%', \
        '        \setlength{\parskip}{0pt}}%', \
        '    {\end{itemize}}', \
        '% Compact headings', \
        '\setlength{\parskip}{0pt}', \
        '\setlength{\headsep}{0pt}', \
        '\setlength{\partopsep}{0pt}', \
        '\setlength{\\topskip}{0pt}', \
        '\usepackage[compact]{titlesec}', \
        '    \\titlespacing{\section}{2pt}{*1}{*1}', \
        '    \\titlespacing{\subsection}{2pt}{*0}{*0}', \
        '', \
        '\\begin{document}', \
        '\section*{TBS Derivatised QC5 Mix}', \
        ''.join(['Total Ion Chromatogram for \\verb|', andi_base, \
            '|']), \
        '\\begin{figure}[h]', \
        ''.join([' '*4, '\centering']), \
        ''.join([' '*8, '\includegraphics[width=0.9\\textwidth]{%']), \
        ''.join([' '*12, \
            'images/', andi_base, '/', andi_base, '_tic}']), \
        '\end{figure}', \
        '\\\\', \
        ''.join(['\\footnotesize Please review and compare the TIC ', \
            'visually, as no automated test ', \
            'has yet been developed.']), \
        '\end{document}' \
    ]
    
    # Use .cdf file name in report filename
    andi_file_noext = os.path.splitext(os.path.basename(andi_file))[0]
    fname_tex = ''.join(\
        [report_folder, 'GCQC_Report_', andi_file_noext,'.tex']\
    )
    
    # Write the .tex file    
    with open(fname_tex, 'wb') as f:
        for line in tex_doc:
            f.write(line)
            f.write('\n')
    
    f.close()

def qc_report(test_report, watch_points, maint_reqs, symm, params, \
        andi_file, report_folder, timecode=time.strftime('%Y%m%d_%H%M')\
    ):
    '''
    @summary:    Creates a .tex report
    
    @param test_report:    Pass/Watch/Fail results for each test
    @type test_report:    ListType (string)
    @param watch_points:    QC items that are close to fail values
    @type watch_points:    ListType (string)
    @param maint_reqs:    QC items that have failed
    @type maint_reqs:    ListType (string)
    @param symm:    Results of peak symmetry testing
    @type symm:    ListType (string, float)
    @param params:    Current parameter settings for processing (see ./proc.py)
    @type params:    ListType
    @param andi_file:    Name of file being processed
    @type andi_file:    StringType
    @param report_folder:    Where reports will be stored
    @type report_folder:    StringType
    @param timecode:    Time stamp 
    @type timecode:    StringType
    '''
    # Use timecode provided
    #report_date = time.strftime('%Y%m%d_%H%M')
    report_date = timecode
    andi_base = os.path.basename(andi_file)[:-4]
    # Add in document headers
    tex_head = [ \
        ''.join(['% report_time = ',report_date]), \
        '\documentclass[10pt,a4paper]{report}', \
        '%%% Put title and MA logo in document header', \
        '\usepackage{graphicx}', \
        # The fancy header package is required. If not installed, report will
        # not be generated. This is part of the texlive-full package on most
        # linux distros. If not available, or not installed properly, the 
        # package can be downloaded from the ctan site at:  
        # http://www.ctan.org/tex-archive/macros/latex/contrib/fancyhdr/
        '\usepackage{fancyhdr}', \
        '\usepackage[%', \
        '    top=1.5cm,bottom=3.5cm,left=1.5cm,right=1.5cm,%', \
        '    paper=a4paper,includeheadfoot]{geometry}', \
        '% Allow multiple columns', \
        '\usepackage{multicol}', \
        '\headheight=72pt', \
        '\headsep=10pt', \
        '\pagestyle{fancy}', \
        '\\renewcommand{\headrulewidth}{0pt}', \
        '\\fancyhead{}', \
        '\\fancyhead[L]{\\textbf{\Large{GC-MS Quality Control Report}}\\\\}',\
        # This is where pdflatex is looking for the logo (used by the
        # report_mailer function)
        ''.join(['\\fancyhead[R]{\includegraphics[height=48pt]', \
            '{images/MA_logo_a3.png}%']), \
        '    \\\\ \\today}', \
        '% Compact lists', \
        '\\newenvironment{itemize*}%', \
        '    {\\begin{itemize}%', \
        '        \setlength{\itemsep}{2pt}%', \
        '        \setlength{\parskip}{0pt}}%', \
        '    {\end{itemize}}', \
        '% Compact headings', \
        '\setlength{\parskip}{0pt}', \
        '\setlength{\headsep}{0pt}', \
        '\setlength{\partopsep}{0pt}', \
        '\setlength{\\topskip}{0pt}', \
        '\usepackage[compact]{titlesec}', \
        '    \\titlespacing{\section}{2pt}{*1}{*1}', \
        '    \\titlespacing{\subsection}{2pt}{*0}{*0}', \
        '', \
        '\\begin{document}', \
        '', \
        '\section*{Test results}', \
        '\\begin{table}[!h]', \
        '    \\begin{tabular*}{0.9\\textwidth}{lccccr}' \
    ]
    
    # Standard tests
    for line in test_report:
        tex_head.append(''.join([ \
            ' '*8, line[0], ' & ', '\ldots', ' & ', '\ldots', ' & ', \
            '\ldots', ' & ', '\ldots', ' & ', line[1], ' \\\\']) \
        )

    
    tex_head.append(''.join([' '*8, '\\\\ \n',' '*8, \
        '\\textbf{Peak Symmetry} \\\\ \n']))
    
    # Peak symmetry tests
    for line in symm:
        try:
            tex_head.append(''.join([' '*8, line[0], \
                ' &  & ', '%0.2f' % (line[1]), \
                ' & ', '%0.2f' % (line[2]), \
                ' & ', '%0.2f' % (line[3]), \
                ' & ', '%0.2f' % (line[4]), \
                ' \\\\'])
            )
        except:
            tex_head.append(''.join([' '*8, line[0], \
                ' & & ', line[1], \
                ' & ', line[2], \
                ' & ', line[3], \
                ' & ', line[4], \
                ' \\\\'])
            )
    
    tex_mid = [ \
        '    \end{tabular*}', \
        '\end{table}', \
        '',
        '\section*{Watch items}' \
    ]
    
    if len(watch_points) != 0:
        tex_mid.append('\\begin{itemize*}')
        for line in watch_points:
            tex_mid.append(''.join([' '*4, '\item ', line]))
        tex_mid.append('\end{itemize*}')
    else:
        tex_mid.append('None')
    
    tex_hindmid = [ \
        '', \
        '    \section*{Maintenance required}' \
    ]
    
    if len(maint_reqs) != 0:
        tex_hindmid.append('\\begin{itemize*}')
        for line in maint_reqs:
            tex_hindmid.append(''.join([' '*4, '\item ', line]))
        tex_hindmid.append('\end{itemize*}')
    else:
        tex_hindmid.append('None')
    
    tex_param = ['', \
        '\\vskip 10pt \\hrule \\vskip 5pt', \
        '', \
        '{\\footnotesize\subsection*{Parameters used:}}', \
        
        '\\begin{table}[!h]', \
        '    \\footnotesize'
        '    \\begin{tabular*}{0.9\\textwidth}{lcr}']
    for line in params:
        tex_param.append(''.join([' '*8, line]))
    
    for line in ['    \end{tabular*}', '\end{table}']:
        tex_param.append(line)
    
#    for cpd in cpd_list:
#        tex_sym.append(''.join([''*8, '&', cpd.get_name(), '&', \
#            cpd.get_symmetry, '\\\\']))
#    for line in ['    \end{tabular*}', '\end{table}']:
#        tex_sym.append(line)
    
    # Graphs of test values
    tex_figs = [ \
#        '\end{multicols}', \
        '', \
        '\\pagebreak', '', \
        '\\begin{figure}[p]',\
        ''.join([' '*4, '\centering']),\
        ''.join([' '*4, '\\begin{tabular}{cc}']),\
        ''.join([' '*8, '\includegraphics[width=0.4\\textwidth]{%']),\
        ''.join([' '*12, \
#            'images/', report_date, '/ribitol_intensity} &']),\
            'images/', andi_base, '/ribitol_intensity} &']),\
        ''.join([' '*8, '\includegraphics[width=0.4\\textwidth]{%']),\
        ''.join([' '*12, \
#            'images/', report_date, '/putrescine_ribitol_ratio_(i_i)} \\\\']),\
            'images/', andi_base, '/putrescine_ribitol_ratio_(i_i)} \\\\']),\
        ''.join([' '*8, '\includegraphics[width=0.4\\textwidth]{%']),\
        ''.join([' '*12, \
#            'images/', report_date, '/ribose_intensity} &']),\
            'images/', andi_base, '/ribose_intensity} &']),\
        ''.join([' '*8, '\includegraphics[width=0.4\\textwidth]{%']),\
        ''.join([' '*12, \
#            'images/', report_date, '/glutamate_ribose_ratio_(i_i)} \\\\']),\
            'images/', andi_base, '/glutamate_ribose_ratio_(i_i)} \\\\']),\
        ''.join([' '*8, '\includegraphics[width=0.4\\textwidth]{%']),\
        ''.join([' '*12, \
#            'images/', report_date, '/asparagine_ribose_ratio_(i_i)} &']),\
            'images/', andi_base, '/asparagine_ribose_ratio_(i_i)} &']),\
        ''.join([' '*8, '\includegraphics[width=0.4\\textwidth]{%']),\
        ''.join([' '*12, \
#            'images/', report_date, '/alanine_valine_ratio_(i_i)} \\\\']),\
            'images/', andi_base, '/alanine_valine_ratio_(i_i)} \\\\']),\
        ''.join([' '*8, '\includegraphics[width=0.4\\textwidth]{%']),\
        ''.join([' '*12, \
#            'images/', report_date, '/alanine_intensity} &']),\
            'images/', andi_base, '/alanine_intensity} &']),\
#        ''.join([' '*8, '\includegraphics[width=0.4\\textwidth]{%']),\
#        ''.join([' '*12, \
#            'images/',report_date,'/valine_area_intensity_ratio} \\\\']),\
        ''.join([' '*8, '\includegraphics[width=0.4\\textwidth]{%']),\
        ''.join([' '*12, \
#            'images/', report_date, '/maltotriose_intensity} \\\\']),\
            'images/', andi_base, '/maltotriose_intensity} \\\\']),\
        ''.join([' '*4, '\end{tabular}']),\
        '\end{figure}', \
        '', '\clearpage', \
#        '\\pagebreak' \
        ]
    
    # End document
    tex_tail = ['', \
        #'}', \
        '\end{document}']
    
    tex_doc = []
    for section in [tex_head, tex_mid, tex_hindmid, \
#        tex_sym, tex_endsym, \
        tex_param, tex_figs, tex_tail]:
        for line in section:
            tex_doc.append(line)
    
    # Use .cdf file name in report filename
    andi_file_noext = os.path.splitext(os.path.basename(andi_file))[0]
    fname_tex = ''.join(\
        [report_folder, 'GCQC_Report_', andi_file_noext,'.tex']\
    )
    
    #fname_tex = ''.join(['GCQC_Report_',time.strftime('%Y%m%d'),'.tex'])
    with open(fname_tex, 'wb') as f:
        for line in tex_doc:
            f.write(line)
            f.write('\n')
    
    f.close()
# EOF
