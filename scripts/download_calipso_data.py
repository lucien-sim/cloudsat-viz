#!/usr/bin/python3

import os

def download_ftp_calipso(ftp_url, dest_path):
    """Downloads the specified files via FTP.
    
    PARAMETERS: 
    -----------

    ftp_url:   FTP address for directory that contains all the 
               CALIPSO files
    dest_path: Path to directory to which the CALIPSO files will be saved. 

    OUTPUTS: 
    --------
    None
    """
    wget_cmd = 'wget -nc -nd -P '+dest_path+' '+ftp_url
    os.system(wget_cmd)
    return None


if __name__ == '__main__':

    message1 = """
Before the download begins, you must order the CALIPSO files from 
the CALIPSO Search and Subsetting web application: 
https://subset.larc.nasa.gov/calipso/login.php

Instructions: 
1. Set up an account (provides username and password)
2. Navigate to the following webpage: 
   https://subset.larc.nasa.gov/calipso/login.php
3. Log in from that page, then use the web application to 
   download the data. 

   You'll need to complete two orders: 

   One for the level one products: 
   Atmospheric, 
   1064 nm calibration/backscatter

   And one for the level two products: 
   Vertical Feature Mask: Layer Information

When the ordered data becomes available, you will receive an 
email for each order and will be ready to proceed. When 
When you receive the emails with the download
details...
"""

    message2 = "Great! Now you're ready to download the CALIPSO data."

    message3 = """
The emails should have included links to the FTP directories that are 
holding the data you ordered. 

Example: ftp://xfr139.larc.nasa.gov/sflops/Distribution/2019183112143_43051
"""

    from global_variables import data_path

    print(message1)
    input1 = input("Enter any letter and hit 'ENTER': ")

    print(message2)

    print(message3)
    ftp_url1 = input("Enter the first link here: ")+'/*'
    ftp_url2 = input("Enter the second link here: ")+'/*'

    dest_path = os.path.join(data_path, 'calipso')

    print("""
Downloading CALIPSO data from the first FTP directory! 
Data will be saved to the following local directory:
""")
    print(dest_path)
    download_ftp_calipso(ftp_url1, dest_path)

    print("""
Downloading CALIPSO data from the second FTP directory! 
Data will be saved to the following local directory:
""")
    print(dest_path)
    download_ftp_calipso(ftp_url2, dest_path)
