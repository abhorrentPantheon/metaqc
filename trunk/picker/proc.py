'''metaqc.picker.proc.py
'''
### Depending on sys, if libnetcdf error occurs, try this before running
# export LD_LIBRARY_PATH=/usr/local/lib:/usr/lib

import os
import shutil
import sqlite3
import subprocess
import sys
import time

sys.path.append('/x/PyMS/')
sys.path.append('/x/metaqc/')

from mailer.Function import create_mail, send_mail
from pylib.Function import config_reader, writedb2csv

#
#    Use the same config file for all
#
config_set = config_reader('/x/metaqc/default.cfg')
# Get position numbers for section values
g_ind = config_set[0]['general']
m_ind = config_set[0]['mail']
config_vals, mail_settings = config_set[g_ind], config_set[m_ind]
preidf = config_vals['instr_data_folders'].split('\n')
instr_data_folders = []
for val in preidf:
    if os.path.exists(val):
        instr_data_folders.append(val)

# Specify reports folder
report_folder = config_vals['report_folder']

# Database file
db_file = config_vals['db_file']

#
#    Email settings
#
# Package maintainer (for mailing error logs etc)
maint_addr = mail_settings['maintainer']
sender = mail_settings['sender']
recipients = mail_settings['recipients']
server = mail_settings['server']


# Other settings
curr_dir = os.getcwd()
mailer_dir = '/x/metaqc/mailer/'
qc_dir = '/x/metaqc/gcqc/'

##############################################################################
# pragma db_file:
## pragma table_info(qc_files)
#(0, u'fkey', u'integer', 1, None, 1)        # identifier
#(1, u'file_name', u'varchar', 0, None, 0)   # QC file name (.cdf)
#(2, u'method_file', u'varchar', 0, None, 0) # Method used for qc file
#(3, u'instr_name', u'varchar', 0, None, 0)  # Instrument name
#(4, u'expr_date', u'datetime', 0, None, 0)  # QC file time stamp (run date)
#
## pragma table_info(qc_data)
#(0, u'file_name', u'varchar', 0, None, 0)   # as for qc_files[1]
#(1, u'expr_date', u'varchar', 0, None, 0)   # as for qc_files[3]
#(2, u'compound', u'varchar', 0, None, 0)    # compound name
#(3, u'area', u'real', 0, None, 0)           # area value
#(4, u'intensity', u'real', 0, None, 0)      # intensity value
#(5, u'RT', u'real', 0, None, 0)             # retention time (found)
#(6, u'RTdelta', u'real', 0, None, 0)        # retention time shift vs expected
#(7, u'symmetry', u'real', 0, None, 0)       # symmetry score
#
## pragma table_info(qc_reports)
#(0, u'file_name', u'varchar', 0, None, 0)   # qc file name (.cdf)
#(1, u'data_folder', u'varchar', 0, None, 0) # where .cdf file located (dir)
#(2, u'rep_loc', u'varchar', 0, None, 0)     # full path to .tex
#
## pragma table_info(pdf_info)
#(0, u'file_name', u'varchar', 0, None, 0)   # qc file name (.cdf)
#(1, u'pdf_loc', u'varchar', 0, None, 0)     # full path to .pdf
#
##############################################################################

#
#    Functionality begins here
#
# If database doesn't exist, create it
if not os.path.isfile(db_file):
    f = open(db_file,'wb')
    f.close()

completed_dct = {}
long_folder_list = []

# Open and close to ensure db not locked for other scripts
conn = sqlite3.connect(db_file)
cur = conn.cursor()
report_pragma = cur.execute('pragma table_info(qc_reports)').fetchall()
if len(report_pragma) != 0:
    raw_comp = cur.execute( \
        'select file_name, rep_loc, data_folder from qc_reports').fetchall()
    for raw in raw_comp:
        completed_dct[str(raw[0])] = str(raw[1])
        long_folder_list.append(str(raw[2]))

cur.close()
conn.close()

# Generate entire list of netCDF data files (with full path data)
file_list = []
# Also create a qc_file:full_path dict
file_dct = {}
for instr in instr_data_folders:
    for root, dirs, files in os.walk(instr):
        for f in files:
            strf = str(os.path.join(root, f))
            # If it's not a cdf, ignore it
            if strf.lower()[-4:] == '.cdf':
                # If it has iQC/QC5 (bio21/botany), add it
                if strf.lower().find('iqc') != -1:
#                if strf.lower().find('qc5') != -1:
                    file_list.append(strf)
                    file_dct[f] = strf

# Use sets to determine those files not yet processed
fset = set(file_dct.keys())
cset = set(completed_dct.keys())
uncomp_set = fset.difference(cset)
uncomp_files = sorted(list(uncomp_set))

os.chdir(qc_dir)
if len(uncomp_set) != 0:
    for f in uncomp_files:
        # If file is > 10min old, process (to ensure no access issues)
        # stat[7] is atime (access) - we want to use ctime (create) [9]
        time_diff = time.time() - os.stat(file_dct[f])[9]
        if time_diff >= (10*60):
            run_call = [sys.executable, '/x/metaqc/gcqc/proc.py', \
                '-c', '/x/metaqc/default.cfg', file_dct[f]]
            subprocess.call(run_call, stdout=open('./output_log.txt', 'wb'), \
                stderr=open('./error_log.txt', 'wb'))
            
            # Update db if no errors
            #    Old method is now deprecated due to generic 
            #    pycdf warnings on some systems.  If it's a proper
            #    error, it will have 'Traceback' in the file
            #if os.stat('./error_log.txt')[6] == 0:
            checkstr = ''
            with open('./error_log.txt', 'rb') as errf:
                for line in errf:
                    checkstr = checkstr+line
            errf.close()
            
            if 'traceback' not in checkstr.lower():
                conn = sqlite3.connect(db_file)
                cur = conn.cursor()
                try:
                    cur.execute(''.join([ \
                        'insert into qc_reports ', \
                            '(file_name, data_folder, rep_loc) ', \
                            'values (?, ?, ?)']), \
                        (f, os.path.dirname(file_dct[f]), ''.join([ \
                            report_folder, \
                            'GCQC_Report_', \
                            os.path.splitext(os.path.basename(f))[0], \
                            '.tex']) \
                        ) \
                    )
                except:
                    cur.execute(''.join([ \
                         'create table if not exists qc_reports ', \
                         '(file_name varchar, data_folder varchar, ', \
                         'rep_loc varchar)']))
                    cur.execute(''.join([ \
                        'insert into qc_reports ', \
                            '(file_name, data_folder, rep_loc) ', \
                            'values (?, ?, ?)']), \
                        (f, os.path.dirname(file_dct[f]), ''.join([ \
                            report_folder, \
                            
                            'GCQC_Report_', \
                            os.path.splitext(os.path.basename(f))[0], \
                            '.tex']) \
                        ) \
                    )
                # Commit changes and close cursor & connection
                conn.commit()
                cur.close()
                conn.close()
                
            # If there was an error, mail it to the maintainer
            else:
                # Join output and error for more info
                outlogstr = ''
                with open('./output_log.txt', 'rb') as logf:
                    for line in logf:
                        outlogstr = outlogstr+line
                logf.close()
                
                fullstr = ''.join(['Output Log:\n', \
                    '='*11, \
                    '\n', \
                    outlogstr, \
                    '\nError Log:\n', \
                    '='*10, \
                    '\n', \
                    checkstr])
                with open('./log_file.txt', 'wb') as outf:
                    outf.write(fullstr)
                outf.close()
                
                # If the error was no peaks, check db for report
                if 'no peaks' in fullstr.lower():
                    conn = sqlite3.connect(db_file)
                    cur = conn.cursor()
                    # Check if a report has already been created
                    try:
                        rep_loc = cur.execute(''.join([ \
                            'select rep_loc from qc_reports where ', \
                            'file_name="', f, '"'])).fetchall()
                    # If db is empty, it will throw one of these
                    except sqlite3.OperationalError:
                        rep_loc = []
                    
                    # If no report is in database, update db with new files
                    if len(rep_loc) == 0:
                        try:
                            cur.execute(''.join([ \
                                'insert into qc_reports ', \
                                    '(file_name, data_folder, rep_loc) ', \
                                    'values (?, ?, ?)']), \
                                (f, os.path.dirname(file_dct[f]), ''.join([ \
                                    report_folder, \
                                    'GCQC_Report_', \
                                    os.path.splitext(os.path.basename(f))[0], \
                                    '.tex']) \
                                ) \
                            )
                        except:
                            cur.execute(''.join([ \
                                 'create table if not exists qc_reports ', \
                                 '(file_name varchar, data_folder varchar, ', \
                                 'rep_loc varchar)']))
                            cur.execute(''.join([ \
                                'insert into qc_reports ', \
                                    '(file_name, data_folder, rep_loc) ', \
                                    'values (?, ?, ?)']), \
                                (f, os.path.dirname(file_dct[f]), ''.join([ \
                                    report_folder, \
                                    'GCQC_Report_', \
                                    os.path.splitext(os.path.basename(f))[0], \
                                    '.tex']) \
                                ) \
                            )
                    # Commit changes and close cursor & connection
                    conn.commit()
                    cur.close()
                    conn.close()
                
                else:
                    # Create mail
                    msg = create_mail(sender, maint_addr, server, \
                        'iQC processing error - GCQC @ Bio21', \
                        '$$'.join(['GCQC Reporting Tool', f, 
                            'There has been an error processing this file.', \
                            'Log file attached.']), \
                        './log_file.txt')
                    send_mail(sender, maint_addr, server, msg)
            
            # Clean up log files
            rmf = ['./error_log.txt', 'output_log.txt']
            if os.path.isfile('./log_file.txt'):
                rmf.append('./log_file.txt')
            for frm in rmf:
                os.remove(frm)

# Return
os.chdir(curr_dir)

#
#    PDF Generator/Mailer
#
# To generate a report for the last iQC file in a directory:
# Get unique folder names
# folder_list_set = set(long_folder_list)
# new_folders = []
# for f in uncomp_files:
#     new_folders.append(os.path.dirname(file_dct[f]))
# 
# new_set = set(new_folders)
# folder_set = folder_list_set.union(new_set)
# folder_list = list(folder_set)
#
#pdf_req = []
#for dirt in folder_list:
#    mtimes = {}
#    f_list = os.listdir(dirt)
#    for f in f_list:
#        if f.lower()[-4:] == '.cdf':
#            if f.lower().find('iqc') == -1:
#                # stat[8] is modification time
#                mtimes[os.stat(file_dct[f])[8]] = os.path.basename(f)
#    # Append the file with the most recent (i.e. last)
#    # modification time to pdf_req
#    pdf_req.append(mtimes[max(mtimes.keys())])

# Otherwise, create a report for each:
pdf_req = fset

conn = sqlite3.connect(db_file)
cur = conn.cursor()
pdf_pragma = cur.execute('pragma table_info(pdf_info)').fetchall()
pdf_done = []
if len(pdf_pragma) != 0:
    pdfraw = cur.execute('select file_name from pdf_info').fetchall()
    for pr in pdfraw:
        pdf_done.append(pr[0]) 

cur.close()
conn.close()

done_set = set(pdf_done)
req_set = set(pdf_req)
need_pdf = sorted(list(req_set.difference(done_set)))

os.chdir(mailer_dir)
if len(need_pdf) != 0:
    for f in need_pdf:
        run_call = [sys.executable, '/x/metaqc/mailer/proc.py', \
            '-m', '-s', '-r', '-c', '/x/metaqc/default.cfg', file_dct[f]] 
        subprocess.call(run_call, stdout=open('./output_log.txt', 'wb'), \
            stderr=open('./error_log.txt', 'wb'))
        
        # Update database file if there are no errors
        checkstr = ''
        with open('./error_log.txt', 'rb') as errf:
            for line in errf:
                checkstr = checkstr+line
        errf.close()
        
        if 'traceback' not in checkstr.lower():
            conn = sqlite3.connect(db_file)
            cur = conn.cursor()
            try:
                cur.execute(''.join([ \
                    'insert into pdf_info (file_name, pdf_loc) ', \
                    'values (?, ?)']), \
                    (f, ''.join([ \
                        report_folder, \
                        'GCQC_Report_', \
                        os.path.splitext(os.path.basename(f))[0], '.pdf'])\
                    ) \
                )
            except:
                cur.execute(''.join(['create table if not exists pdf_info ', \
                    '(file_name varchar, pdf_loc varchar)']))
                cur.execute(''.join([ \
                    'insert into pdf_info (file_name, pdf_loc) ', \
                    'values (?, ?)']), \
                    (f, ''.join([ \
                        report_folder, \
                        'GCQC_Report_', \
                        os.path.splitext(os.path.basename(f))[0], '.pdf'])\
                    ) \
                )
            # Commit and close
            conn.commit()
            cur.close()
            conn.close()
        
        # Otherwise mail the error log
        else:
            # Join output and error for more info
            outlogstr = ''
            with open('./output_log.txt', 'rb') as logf:
                for line in logf:
                    outlogstr = outlogstr+line
            logf.close()
            
            fullstr = ''.join(['Output Log:\n', \
                '='*11, \
                '\n', \
                outlogstr, \
                '\nError Log:\n', \
                '='*10, \
                '\n', \
                checkstr])
            with open('./log_file.txt', 'wb') as outf:
                outf.write(fullstr)
            outf.close()
            
            # Create mail
            msg = create_mail(sender, maint_addr, server, \
                'iQC processing error - Report @ Bio21', \
                '$$'.join(['GCQC Reporting Tool', f, 
                    'There has been an error processing this file.', \
                    'Log file attached.']), \
                './log_file.txt')
            send_mail(sender, maint_addr, server, msg)
        
        # Clean up log files
        rmf = ['./error_log.txt', 'output_log.txt']
        if os.path.isfile('./log_file.txt'):
            rmf.append('./log_file.txt')
        for frm in rmf:
            os.remove(frm)

# Update the copy of the database
db_path = os.path.split(db_file)[0]
db_copy = os.path.join(db_path, 'metaqc_gcqc_copy.db')
shutil.copy2(db_file, db_copy)
# Read and execute only on file (except owner - need to be able to update!)
os.chmod(db_copy, 0755)
# Write out backup .csv files
writedb2csv(db_file, '../metaqc_dbexport')
# Bio21
shutil.copy2('../metaqc_dbexport_files.csv', '/mnt/ma_srv/Shared_Data/Projects/Metabolomics Australia/GC-MS QC Method development/metaqc_dbexport_files.csv')
shutil.copy2('../metaqc_dbexport_data.csv', '/mnt/ma_srv/Shared_Data/Projects/Metabolomics Australia/GC-MS QC Method development/metaqc_dbexport_data.csv')
shutil.copy2('../metaqc_dbexport_report.csv', '/mnt/ma_srv/Shared_Data/Projects/Metabolomics Australia/GC-MS QC Method development/metaqc_dbexport_report.csv')
# Botany
#shutil.copy2('../metaqc_dbexport_files.csv', '/media/sf_main/groups/Metabolomics/QC REPORTS/metaqc_dbexport_files.csv')
#shutil.copy2('../metaqc_dbexport_data.csv', '/media/sf_main/groups/Metabolomics/QC REPORTS/metaqc_dbexport_data.csv')
#shutil.copy2('../metaqc_dbexport_report.csv', '/media/sf_main/groups/Metabolomics/QC REPORTS/metaqc_dbexport_report.csv')


os.chdir(curr_dir)
# EOF
