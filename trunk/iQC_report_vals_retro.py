# Collect vals for creating table
import os, sqlite3, sys
sys.path.append('/x/PyMS/')
sys.path.append('/x/iQC/')

#from operator import itemgetter ## 
from gcqc.Function import get_cpd_data
from pylib.Function import config_reader, writedb2csv, writecsvout

db_file = 'iQC_gcqc.db'
#db_file = 'iQC2.db'
config_file = 'default.cfg'
conn = sqlite3.connect(db_file)
cur = conn.cursor()
# Don't include TBS files
raw_flist = cur.execute(''.join([ \
    'select file_name from qc_files where method_file not like "%TBS%" ', \
    'order by expr_date asc' ])).fetchall()

flist = []
for ii in range(len(raw_flist)):
    flist.append(str(raw_flist[ii][0]))

fdic = {}
for f in flist:
    # Get all prev vals for calc tests
    raw_fetch = cur.execute(''.join([ \
        'select compound, area, intensity from qc_data where ', \
        'compound in (', \
            '"Ribitol", "Putrescine 4TMS", "Ribose", "Glutamate 3TMS", ', \
            '"Asparagine 3TMS", "Alanine 2TMS", "Valine 2TMS", ' \
            '"Maltotriose MX1"',\
        ') and file_name="', f,
    '"'])).fetchall()
    rt_lock = cur.execute(''.join([ \
        'select RT, RTdelta from qc_data where file_name="', f, \
    '" and compound="Mannitol"'])).fetchall()
    
    cdic = {}
    for ii in range(len(raw_fetch)):
        # Compound name
        nm = str(raw_fetch[ii][0])
        # Area
        try:
            ar = float(raw_fetch[ii][1])
        except:
            ar = 0.0
        # Intensity
        try:
            its = float(raw_fetch[ii][2])
        except:
            its = 0.0
        
        tmp = [nm, ar, its]
        # RTdelta is in seconds, while RT is in minutes
        try:
            cdic['mtl_delta'] = float(rt_lock[0][1])/60
        except ValueError:
            cdic['mtl_delta'] = rt_lock[0][1]
        try:
            cdic['mtl_rt'] = float(rt_lock[0][0])
        except ValueError:
            cdic['mtl_rt'] = rt_lock[0][0]
        
        if tmp[0] == 'Ribitol':
            cdic['rbtl'] = tmp[2]
        elif tmp[0] == 'Putrescine 4TMS':
            cdic['putr'] = tmp[2]
        elif tmp[0] == 'Ribose':
            cdic['rbos'] = tmp[2]
        elif tmp[0] == 'Glutamate 3TMS':
            cdic['glu'] = tmp[2]
        elif tmp[0] == 'Asparagine 3TMS':
            cdic['asn'] = tmp[2]
        elif tmp[0] == 'Alanine 2TMS':
            cdic['ala'] = tmp[2]
        elif tmp[0] == 'Valine 2TMS':
            cdic['val'] = tmp[2]
        elif tmp[0] == 'Maltotriose MX1':
            cdic['mlt'] = tmp[2]
    
    fdic[f] = cdic


cur.close()
conn.close()

config_set = config_reader(config_file)
c_ind = config_set[0]['general']
l_ind = config_set[0]['test_names']
v_ind = config_set[0]['pf_vals']

config_vals = config_set[c_ind]
report_labels = config_set[l_ind]
report_vals = config_set[v_ind]

# Reformat strings for report vals
for key in report_vals.keys():
    report_vals[key] = float(report_vals[key])

tdic = {}
for f in flist:
    #
    #    Retention Time Lock
    #
    # 25deg locks mannitol at 11 min
    ############################################# check - potentially borky time?
    try:
        rt_lock_val = fdic[f]['mtl_delta']
    except KeyError:
        rt_lock_val = None
        rt_lock_test = 'FAIL'
    
    if rt_lock_val != None:
        if 10.0 <= fdic[f]['mtl_rt'] <= 12.0:
            report_vals.update(rt_lock_pass=0.025)
            report_vals.update(rt_lock_fail=0.05)
        # 15deg locks mannitol at 15.5 min
        elif 14.0 <= fdic[f]['mtl_rt'] <= 17.0:
            report_vals.update(rt_lock_pass=0.0325)
            report_vals.update(rt_lock_fail=0.075)
        # 7deg locks mannitol at 21.5 min
        elif 19.5 <= fdic[f]['mtl_rt'] <= 23.5:
            report_vals.update(rt_lock_pass=0.05)
            report_vals.update(rt_lock_fail=0.10)
        else:
            rt_lock_test = 'FAIL'
        # Test
        rtl_pass = report_vals['rt_lock_pass']
        rtl_fail = report_vals['rt_lock_fail']
        
        if rt_lock_val < rtl_pass:
            rt_lock_test = 'PASS'
        elif rtl_pass <= rt_lock_val <= rtl_fail:
            rt_lock_test = 'WATCH'
        else:
            rt_lock_test = 'FAIL'
    #
    #    Inlet conditions
    #
    # Ribitol intensity test
    rbtl_i_pass = report_vals['rbtl_i_pass']
    rbtl_i_fail = report_vals['rbtl_i_fail']
    try:
        rbtl_i = fdic[f]['rbtl']
    except KeyError:
        rbtl_i = None
        rbtl_i_test= 'FAIL'
    
    if rbtl_i != None:
        if rbtl_i > rbtl_i_pass:
            rbtl_i_test = 'PASS'
        elif rbtl_i_fail <= rbtl_i <= rbtl_i_pass:
            rbtl_i_test = 'WATCH'
        else:
            rbtl_i_test = 'FAIL'
    
    # Putrescine/Ribitol ratio test
    pr_pass = report_vals['pr_r_pass']
    pr_fail = report_vals['pr_r_fail']
    try:
        pr_val = fdic[f]['putr']/cdic['rbtl']
    except KeyError, ZeroDivisionError:
        pr_val = None
        pr_r = 'FAIL'
    
    if pr_val != None:
        if pr_val > pr_pass:
            pr_r = 'PASS'
        elif pr_fail <= pr_val <= pr_pass:
            pr_r = 'WATCH'
        else:
            pr_r = 'FAIL'
    
    # Ribose intensity test
    rbose_pass = report_vals['rbose_i_pass']
    rbose_fail = report_vals['rbose_i_fail']
    try:
        rbose_i = fdic[f]['rbos']
    except KeyError:
        rbose_i = None
        rbose_i_test = 'FAIL'
    
    if rbose_i != None:
        if rbose_i > report_vals['rbose_i_pass']:
            rbose_i_test = 'PASS'
        elif rbose_fail <= rbose_i <= rbose_pass:
            rbose_i_test = 'WATCH'
        else:
            rbose_i_test = 'FAIL'
    
    # Glutamate/Ribose ratio test
    gr_pass = report_vals['gr_r_pass']
    gr_fail = report_vals['gr_r_fail']
    try:
        gr_val = fdic[f]['glu']/cdic['rbos']
    except KeyError, ZeroDivisionError:
        gr_val = None
        gr_r = 'FAIL'
    
    if gr_val != None:
        if gr_val > gr_pass:
            gr_r = 'PASS'
        elif gr_fail <= gr_val <= gr_pass:
            gr_r = 'WATCH'
        else:
            gr_r = 'FAIL'
    
    # Asparagine/Ribose ratio
    ar_pass = report_vals['ar_r_pass']
    ar_fail = report_vals['ar_r_fail']
    try:
        ar_val = fdic[f]['asn']/cdic['rbos']
    except KeyError, ZeroDivisionError:
        ar_val = None
        ar_r = 'FAIL'
    
    if ar_val != None:
        if ar_val > ar_pass:
            ar_r = 'PASS'
        elif ar_fail <= ar_val <= ar_pass:
            ar_r = 'WATCH'
        else:
            ar_r = 'FAIL'
    
    #
    #    Column Conditions
    #
    # Alanine intensity test - Column condition - start (1)
    ala_pass = report_vals['ala_i_pass']
    ala_fail = report_vals['ala_i_fail']
    try:
        ala_i = fdic[f]['ala']
    except KeyError:
        ala_i = None
        ala_i_test = 'FAIL'
    
    if ala_i != None:
        if ala_i > ala_pass:
            ala_i_test = 'PASS'
        elif ala_fail <= ala_i <= ala_pass:
            ala_i_test = 'WATCH'
        else:
            ala_i_test = 'FAIL'
    
    # Column condition - start (2)
    av_pass = report_vals['av_r_pass']
    av_fail = report_vals['av_r_fail']
    
    try:
        av_val = fdic[f]['ala']/cdic['val']
    except KeyError, ZeroDivisionError:
        av_val = None
        av_r = 'FAIL'
    
    if av_val != None:
        if av_val > av_pass:
            av_r = 'PASS'
        elif av_fail <= av_val <= av_pass:
            av_r = 'WATCH'
        else:
            av_r = 'FAIL'
    
    # Column condition - end
    try:
        mlt_i = fdic[f]['mlt']
    except KeyError:
        mlt_i = None
        mlt_i_test = 'FAIL'
    
    if mlt_i != None:
        if mlt_i > report_vals['mlt_i_pass']:
            mlt_i_test = 'PASS'
        elif report_vals['mlt_i_fail'] <= mlt_i <= report_vals['mlt_i_pass']:
            mlt_i_test = 'WATCH'
        else:
            mlt_i_test = 'FAIL'
    
    # Assemble outlist
    tdic[f] = [rt_lock_val, rt_lock_test, rbtl_i, rbtl_i_test, pr_val, pr_r, \
        rbose_i, rbose_i_test, gr_val, gr_r, ar_val, ar_r, ala_i, ala_i_test, \
        av_val, av_r, mlt_i, mlt_i_test]

ordered_keys = ['rt_lock', 'rbtl_i', 'pr_r', 'rbose_i', 'gr_r', \
    'ar_r', 'ala_i', 'av_r', 'mlt_i']

######################
# Update db
conn = sqlite3.connect(db_file)
cur = conn.cursor()

for f in flist:
    try:
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
            (f, tdic[f][0], tdic[f][1], tdic[f][2], tdic[f][3], tdic[f][4], \
                tdic[f][5], tdic[f][6], tdic[f][7], tdic[f][8], tdic[f][9], \
                tdic[f][10], tdic[f][11], tdic[f][12], tdic[f][13], \
                tdic[f][14], tdic[f][15], tdic[f][16], tdic[f][17] \
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
                    report_labels[ordered_keys[ii]]\
                )\
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
            (f, tdic[f][0], tdic[f][1], tdic[f][2], tdic[f][3], tdic[f][4], \
                tdic[f][5], tdic[f][6], tdic[f][7], tdic[f][8], tdic[f][9], \
                tdic[f][10], tdic[f][11], tdic[f][12], tdic[f][13], \
                tdic[f][14], tdic[f][15], tdic[f][16], tdic[f][17] \
            )
        )

conn.commit()
cur.close()
conn.close()

###################
#    Write output table to .csv
#
# Headers
outdata = [['Filename', 'MethodFile', 'Instrument', 'ExperimentDate', \
    'Test', 'Value', 'Pass/Watch/Fail']]
#    
#    
#    'Retention Time Lock delta', 'Ribitol intensity', \
#    'Putrescine/Ribitol ratio', 'Ribose intensity', 'Glutamate/Ribose ratio', \
#    'Asparagine/Ribose ratio', 'Alanine intensity', 'Alanine/Valine ratio', \
#    'Maltotriose intensity']]

# Get data
conn = sqlite3.connect(db_file)
cur = conn.cursor()
raw_fetch = cur.execute( \
    'select * from qc_files natural join qc_results').fetchall()
cur.close()
conn.close()

# Tidy up
for rec in raw_fetch:
#    val = [str(rec[1]), str(rec[2]), str(rec[3]), str(rec[4]), 'Values', \
#        rec[5], rec[7], rec[9], rec[11], rec[13], rec[15], rec[17], \
#        rec[19], rec[21]]
#    pwf = [str(rec[1]), str(rec[2]), str(rec[3]), str(rec[4]), 'PWF', \
#        str(rec[6]), str(rec[8]), str(rec[10]), str(rec[12]), \
#        str(rec[14]), str(rec[16]), str(rec[18]), str(rec[20]), str(rec[22])]
#    outdata.append(val)
#    outdata.append(pwf)
    n = 5
    for ii in range(len(ordered_keys)):
        tmp = [str(rec[1]), str(rec[2]), str(rec[3]), str(rec[4]), \
            report_labels[ordered_keys[ii]], rec[ii+n], rec[ii+n+1]]
        n += 1
        outdata.append(tmp)

# Write output file
writecsvout(outdata, 'iQC_test_db_write.csv')









