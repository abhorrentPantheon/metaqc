'''iQC.mailer.Function.py
'''

import os
import smtplib
import string

from email import encoders
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.application import MIMEApplication
from email.Utils import formatdate

def create_mail(sender, recipients, server, subject, body, *attachment):
    '''
    @summary:    Creates an email message
    
    @param sender:    Address mail is sent from
    @type sender:    StringType
    @param recipients:    To whom the mail is sent
    @type recipients:    StringType
    @param server:    Mail server that message is sent through
    @type server:    StringType
    @param subject:    Contents of Subject: field for the message
    @type subject:    StringType
    @param body:    Message text to send. Line breaks in final message body 
        are indicated by '$'
    @type body:    StringType
    @param attachment:    Filename of the file to attach to the message
    @type attachment:    StringType
    
    @return:    A message object
    @rtype:    MIMEMultipartType
    '''
    
    # Create the base of the message
    msg = MIMEMultipart()
    msg['From'] = sender
    msg['To'] = recipients
    msg['Subject'] = subject
    msg['Date'] = formatdate(localtime=True)
    
    # Create the message body in both plain text and html
    sub_msgt = ['']
    sub_msgh = [''.join(['<html>\n', ' '*4, '<head></head>\n', ' '*4, \
        '<body>\n', ' '*8, '<p>\n'])]
    
    for msgl in string.split(body, '$'):
        sub_msgt.append(string.join([msgl, '\n'], ''))
        sub_msgh.append(string.join([' '*12, msgl, '<br>\n'], ''))
    
    sub_msgh.append(''.join([' '*8, '</p>\n', ' '*4, '</body>\n</html>']))
    
    msgt = string.join(sub_msgt, '')
    msgh = string.join(sub_msgh, '')
    
    # Add the text portion of the message in both html and plain text
    msg_alt = MIMEMultipart('alternative')
    
    mtxt = MIMEText(msgt, 'plain')
    mhtm = MIMEText(msgh, 'html')
    
    msg_alt.attach(mtxt)
    msg_alt.attach(mhtm)
    
    if len(attachment) != 0:
        # Check file type
        if attachment[0][-3:].lower() == 'txt':
            log = open(attachment[0], 'rb')
            mlog = MIMEBase('text', 'plain')
            mlog.set_payload(log.read())
            log.close()
            encoders.encode_base64(mlog)
            mlog.add_header('Content-Disposition', 'attachment', \
                filename = os.path.basename(attachment[0]))
            mlog.add_header('Content-type', 'text/plain', \
                filename = os.path.basename(attachment[0]))
            msg.attach(msg_alt)
            msg.attach(mlog)
        elif attachment[0][-3:].lower() == 'pdf':
            # Create the attachment
            rep = open(attachment[0], 'rb')
            #mrep = MIMEApplication(rep.read(), 'pdf')
            #rep.close()
            mrep = MIMEBase('application', 'pdf')
            mrep.set_payload(rep.read())
            rep.close()
            encoders.encode_base64(mrep)
            mrep.add_header('Content-Disposition', 'attachment', \
                filename = os.path.basename(attachment[0]))
            mrep.add_header('Content-type', 'application/pdf', \
                filename = os.path.basename(attachment[0]))
            
            # Build message structure
            msg.attach(msg_alt)
            msg.attach(mrep)
    else:
        msg.attach(msg_alt)
    
    return msg

def send_mail(sender, recipients, server, msg):
    '''
    @summary:    Sends an email message
    
    @param sender:    Address mail is sent from
    @type sender:    StringType
    @param recipients:    To whom the mail is sent
    @type recipients:    StringType
    @param server:    Mail server that message is sent through
    @type server:    StringType
    @param msg:   Message to be sent
    @type msg:   MIMEMultipartType  
    '''
    rcp_list = string.split(recipients, ',')
    
    pyserver = smtplib.SMTP(server)
    pyserver.sendmail(sender, \
        rcp_list, \
        msg.as_string())
    pyserver.quit()

# EOF
