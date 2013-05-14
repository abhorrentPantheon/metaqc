'''metaqc.mailer.proc.py
'''

import glob
import os
import re
import string
import subprocess
import sys

from time import strftime

sys.path.append('/x/metaqc/')
# From module containing generic stuff that can be recycled
from pylib.Function import config_reader
from mailer.Function import create_mail, send_mail

#==========================================================================
#    Parse arguments from command line 
#==========================================================================
## FUTUREPY: argparse is used in python3 
## Tested. Differences to note are that argparse defaults to type=string and
##    action='store'. nargs defaults to the next argument for -x or --xxx 
##    and the argument itself in the case of positional arguments. Specifying
##    nargs here will cause the parser to store the result as a list of 
##    strings, not as a string.
##    Further, the positional argument is different. It is automatically stored
##    as a string with the name specified in the 'name or flag' parameter. To
##    display this differently in help, use metavar to change this.
#import argparse
#parser = argparse.ArgumentParser(
#    description='Send an email with GCQC reporting information.')
#parser.add_argument('-s', '--scan', action='store_true', dest='scan_req', \
#    help='Scan for a fail value in the most recent report file')
#parser.add_argument('-m', '--mail', action='store_true', dest='mail_req', \
#    help='Create an email regardless of whether or not there is a failure')
#parser.add_argument('-r', '--report', action='store_true', dest='report_req', \
#    help='Attach a summary report')
#parser.add_argument('-c', '--config', dest='config_file', metavar='FILE.cfg', \
#    help='''Configuration file to use. Specifying other options
#    on the command line overrides the settings from this file''')
#parser.add_argument('qc_file',  metavar='QC_FILE.CDF', 
#    help='QC file used for report')
## To display help:
##args = parser.parse_args(['-h'])
#
## To run from cmd line:
#args = parser.parse_args()
#==========================================================================

from optparse import OptionParser

usage = 'usage: %prog [options] QC_FILE.CDF'
parser = OptionParser(usage, description=\
    'Send an email with GCQC reporting information.')
parser.add_option('-s', '--scan', action='store_true', dest='scan_req', \
    help='Scan for a fail value in the most recent report file')
parser.add_option('-m', '--mail', action='store_true', dest='mail_req', \
    help='Create an email regardless of whether or not there is a failure')
parser.add_option('-r', '--report', action='store_true', dest='report_req', \
    help='Attach a summary report (.pdf)')
parser.add_option('-c', '--config', nargs=1, action='store', \
    type='string', dest='config_file', metavar='FILE.cfg', \
    help=''.join(['Configuration file to use. Specifying other options', \
    ' on the command line overrides settings in this file']))
# NB: QC_FILE.CDF is used to get the name of the report file.

# TODO: Determine if this is necessary, and implement if so
#parser.add_option("-v", action="store_true", dest="verbose")

(options, args) = parser.parse_args()
#==========================================================================

# Verify there is a QC file specified
if len(args) == 0:
    print ' !! No .cdf file specified; exiting.'
    raise RuntimeError
elif len(args) > 1:
    print ' !! More than one argument (.cdf file) specified; exiting.'
    raise RuntimeError
else:
    qc_file = args[0]

# Read config options
if options.config_file != None:
    config_set = config_reader(options.config_file)
    # Get position numbers for section values
    g_ind = config_set[0]['general']
    m_ind = config_set[0]['mail']
    config_vals, mail_settings = config_set[g_ind], config_set[m_ind]
    
else:
# Specify options for default settings
    config_vals = {}
    config_vals['report_folder'] = '/x/metaqc/reports'
    mail_settings = {}
    mail_settings['sender'] = 'qc_reports@unimelb.edu.au'
    mail_settings['recipients'] = 'jairus.bowne@gmail.com'
    mail_settings['maintainer'] = 'jairus.bowne@gmail.com'
    mail_settings['server'] = 'smtp.unimelb.edu.au'
    mail_settings['subject'] = 'GCQC Report - Test - Please ignore'
    # Message text. Use '$' to specify line breaks.
    mail_settings['msg_body'] = ''.join(['GCQC Report - ', \
        strftime('%d %b %Y'), '$$This is an automatically generated email.'])

# Check that the necessary images and styles are present before continuing
# These are defined in gcqc.Function.py
if not os.path.isfile(os.path.join(config_vals['report_folder'], \
        'images/MA_logo_a3.png')):
    print ' !! The logo that you are planning to use is not present.'
    raise RuntimeError

# TODO: check for other important latex stuff?
if not os.path.isfile(os.path.join(config_vals['report_folder'], \
        'styles/fancyhdr.sty')):
    print ' !! The fancy header style is not installed.'
    print '    Please install this style first, or the texlive-full package.'
    print '    If texlive-full has been installed, then you may wish to try'
    print '    http://www.ctan.org/tex-archive/macros/latex/contrib/fancyhdr/'
    raise RuntimeError

# Generate the email itself.  Insert date if necessary. Swap \n for '$' so
# that the html version is processed correctly in create_mail function.
body = string.replace(mail_settings['msg_body'], '\n', '$')
body = string.replace(body, 'DATE', strftime('%d %b %Y'))
body = ''.join([body, 'iQC file - ', os.path.basename(qc_file), '$'])

# If the mail_on_fail entry exists, scan_req must be true
if config_vals['mail_on_fail'] != '':
    mof_list = string.split(config_vals['mail_on_fail'],'\n')
    try:
        mof_list.remove('')
    except:
        pass
    if len(mof_list) != 0:
        options.scan_req = True

orig_dir = os.getcwd()
report_folder = config_vals['report_folder']

# Scan for fail values if requested
if options.scan_req is True:
    # If not in report directory, go there
    os.chdir(report_folder)
    # Get file name of report file
    qc_bfn = os.path.splitext(os.path.basename(qc_file))[0]
    
    try:
        report_file_tex = glob.glob(''.join(['*', qc_bfn, '.tex']))[0]
    except IndexError:
        print ' !! No .tex report found. Are you certain there is one?'
        raise RuntimeError
    # Collect fail lines, and trim them to only the test name
    fail_vals = []
    for ln in open(report_file_tex):
        if 'FAIL' in ln:
            fail_vals.append(re.sub('\s+(.+)\s&.+&.+&.+&.+&.+\n', '\\1', ln))
    
    if len(fail_vals) != 0:
        body = ''.join([body, '$The following tests have failed:$'])
        
        fail_text = ''
        for ln in fail_vals:
            if ln in mof_list:
                options.mail_req = True
            fail_text = ''.join([fail_text, ' '*4, ln,'$'])
        
        body = ''.join([body, fail_text])

# Return to previous dir
if os.getcwd() != orig_dir:
    os.chdir(orig_dir)

if options.mail_req is True:
    if  options.report_req is True:
        os.chdir(report_folder)
        # Create the pdf
        subprocess.call(['pdflatex', \
            os.path.join(report_folder, report_file_tex)])
        report_file_pdf = string.replace(report_file_tex, '.tex', '.pdf')
        
        # Clean up extra files from pdflatex
        os.remove(string.replace(report_file_tex, '.tex', '.log'))
        os.remove(string.replace(report_file_tex, '.tex', '.aux'))

        # Create the message
        msg = create_mail(mail_settings['sender'], \
            mail_settings['recipients'], \
            mail_settings['server'], \
            mail_settings['subject'], \
            body, \
            report_file_pdf)
    else:
        msg = create_mail(mail_settings['sender'], \
            mail_settings['recipients'], \
            mail_settings['server'], \
            mail_settings['subject'], \
            body)
    
    # Send message
    send_mail(mail_settings['sender'], \
        mail_settings['recipients'], \
        mail_settings['server'], \
        msg)
    
    os.chdir(orig_dir)
    
    print 'Mail sent'

# EOF
