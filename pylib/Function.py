'''metaqc.pylib.Function.py
'''

import ConfigParser
import cPickle
import csv
import sqlite3
import string

def config_reader(config_file):
    '''
    @summary:    Parse a config file to get settings
    
    @param config_file:    Name of settings file to use
    @type config_file:    StringType 
    
    @return:    Settings grouped by section
    @rtype:    TupleType (An (n+1)-tuple of DictTypes). Length is equal to the 
        number of sections in config file plus the names of each section
    '''
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    sects = config.sections()
    
    ii = 1
    sect_headers = {}
    for sect in sects:
        sect_index = ii
        sect_name = sect
        sect_headers[sect_name]=sect_index
        ii += 1
    
    dict_list = [sect_headers]        
    for sect in sects:
        sec_dict = {}
        for opt in config.options(sect):
            sec_dict[opt] = config.get(sect, opt)
        dict_list.append(sec_dict)
    
#    print tuple(dict_list)
    return tuple(dict_list)



def writecsvout(table, filename):
    '''
    @summary:    Ronseal. writes a csv file.
    
    @param table:    Input table
    @type table:    ListType (List of lists)
    @param filename:    Output filename
    @type filename:    StringType
    
    @return:    None
    @rtype:    NoneType
    '''
    
    output = csv.writer(open(filename, "wb"))
    output.writerows(table)

def str2fname(input_string):
    '''
    @summary:    Convert a string into a string that is suitable for 
                 use as a filename
    @param input_string:    String to convert
    @type input_string:    StringType
    '''
    res_char = '</ \\>?%$*:;|"\''
    out_str = []
    for char in input_string:
        if char in res_char:
            out_str.append(''.join('_'))
        else:
            out_str.append(''.join(char))
    
    ret_str = string.lower(''.join(out_str))
    
    return ret_str

def in_pickle(filename):
    """
    @summary:    Reads a pickle object from file
    
    @param filename:    Filename of pickle object
    @type filename:    StringType
    
    @return:    Object that was pickled in file filename
    @rtype:    Arbitrary (defined by pickle object)
    """
    
    fp = open(filename, "rb")
    obj = cPickle.load(fp)
    fp.close()
    return obj

def out_pickle(obj, filename):
    """
    @summary:    Dumps a pickle object to file
    
    @param table:    Name that will be assigned to object
    @type table:    Arbitrary (pickle object)
    @param pkl_in:    Name of pickle object
    @type pkl_in:    StringType
    """
    
    fp = open(filename, "wb")
    cPickle.dump(obj,fp)
    fp.close()

def writedb2csv(db_file, csv_basename):
    '''
    @summary:    Writes the qc_files and qc_data tables to .csv files.
        These are the two tables that are most useful for linking to Excel.
        The name of the output files will be '<csv_file>_files.csv' and 
        '<csv_basename>_data.csv'.
    
    @param db_file:    metaqc database filename
    @type db_file:    StringType
    
    @param csv_basename:    Output .csv file base name
    @type csv_basename:    StringType
    '''
    
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    
    # qc_files and qc_data written for completeness. Merged data is
    # the main file to be used for tracking data in Excel.
    raw_f = cur.execute('select * from qc_files').fetchall()
    raw_d = cur.execute('select * from qc_data').fetchall()
    raw_r = cur.execute( \
        'select * from qc_files natural join qc_results').fetchall()
    
    cur.close()
    conn.close()
    
    clean_f = [['fkey', 'Filename', 'Method', 'Instrument', 'DateRun']]
    for ii in range(len(raw_f)):
        tmp = [int(raw_f[ii][0]),# fkey
            str(raw_f[ii][1]),   # file name
            str(raw_f[ii][2]),   # method name
            str(raw_f[ii][3]),   # instrument name
            str(raw_f[ii][4])]   # experiment date
        clean_f.append(tmp)
    
    clean_d = [['Filename', 'DateRun', 'Compound', 'Area', 'Intensity', 'RT',\
        'RTdelta', 'Symmetry']]
    for ii in range(len(raw_d)):
        tmp = [str(raw_d[ii][0]),# file name
            str(raw_d[ii][1]),   # experiment date
            str(raw_d[ii][2]),   # compound name
            raw_d[ii][3],        # area
            raw_d[ii][4],        # intensity
            raw_d[ii][5],        # retention time
            raw_d[ii][6],        # shift from expected retention time
            raw_d[ii][7]]        # ratio of left:right peak boundary lengths
        clean_d.append(tmp)
    
    clean_r = [['Filename', 'MethodFile', 'Instrument', 'ExperimentDate', \
        'Test', 'Value', 'Pass/Watch/Fail']]
    
    # Report labels & keys outside scope of function; create.
    report_labels = {'ala_i': 'Alanine intensity', 'ar_r': 'Asparagine/Ribose ratio', \
        'av_r': 'Alanine/Valine ratio', 'gr_r': 'Glutamate/Ribose ratio', \
        'mlt_i': 'Maltotriose intensity', 'pr_r': 'Putrescine/Ribitol ratio', \
        'rbtl_i': 'Ribitol intensity', 'rbose_i': 'Ribose intensity', \
        'rt_lock': 'Retention Time Lock delta'}
    ordered_keys = ['rt_lock', 'rbtl_i', 'pr_r', 'rbose_i', 'gr_r', \
        'ar_r', 'ala_i', 'av_r', 'mlt_i']
    # Tidy up results
    for rec in raw_r:
        n = 5
        for ii in range(len(ordered_keys)):
            tmp = [str(rec[1]), str(rec[2]), str(rec[3]), str(rec[4]), \
                report_labels[ordered_keys[ii]], rec[ii+n], rec[ii+n+1]]
            n += 1
            clean_r.append(tmp)
    
    # Files
    csv_f = ''.join([csv_basename, '_files.csv'])
    with open(csv_f, "wb") as f:
        outf = csv.writer(f)
        outf.writerows(clean_f)
    # Data
    csv_d = ''.join([csv_basename, '_data.csv'])
    with open(csv_d, "wb") as d:
        outd = csv.writer(d)
        outd.writerows(clean_d)
    
    # Report Values
    csv_r = ''.join([csv_basename, '_report.csv'])
    with open(csv_r, "wb") as r:
        outd = csv.writer(r)
        outd.writerows(clean_r)

# EOF
