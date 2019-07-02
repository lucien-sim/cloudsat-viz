#!/usr/bin/python3

import os

def download_ftp_cloudsat(ftp_url, usr, pwd, dest_path):
    """Downloads the specified files via FTP.
    
    PARAMETERS: 
    -----------

    ftp_url:   FTP address for directory that contains all the 
               CloudSat files
    usr:       Username for FTP login
    pwd:       Password for FTP login
    dest_path: Path to directory to which the CloudSat files will be saved. 

    OUTPUTS: 
    --------
    None
    """

    wget_cmd = 'wget -nc -nd -P '+dest_path + \
        ' --user='+usr+' --password='+pwd+' '+ftp_url
    os.system(wget_cmd)
    return None


if __name__ == '__main__': 

    message1 = """
Before the download begins, you must order the CloudSat files from the CloudSat 
Data Processing Center: http://www.cloudsat.cira.colostate.edu/

Instructions: 
1. Set up an account (provides username and password)
2. Use the following page to order the desired data: 
    http://www.cloudsat.cira.colostate.edu/order-data. 

    This visualization requires the following data products: 
    2B-CLDCLASS.P1_RO5; 2B-CWC-RO.P1_R05; 2B-GEOPROF.P1_R05. 

When the ordered data becomes available, you will receive an email and 
you're ready to proceed. When you receive an email with the download
details...
"""

    message2 = "Great! Now you're ready to download the data."

    message3 = """
The email should have included a link to the FTP directory that's holding 
the data you ordered. 

Example: ftp://ftp1.cloudsat.cira.colostate.edu/.orders/.5d164bf23187fd351f5e3e10
"""

    from global_variables import data_path

    print(message1)
    input1 = input("Enter any letter and hit 'ENTER': ")

    print(message2)
    usr = input("Enter your Cloudsat Data Processing Center username: ")
    pwd = input("Enter your Cloudsat Data Processing Center password: ")

    print(message3) 
    ftp_url = input("Enter that link here: ")+'/*'

    dest_path = os.path.join(data_path, 'cloudsat')

    print("""
Downloading CloudSat data! 
Data will be saved to the following local directory:
""")
    print(dest_path)
    download_ftp_cloudsat(ftp_url, usr, pwd, dest_path)
