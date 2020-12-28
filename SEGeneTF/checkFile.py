# -*- coding: utf-8 -*-
"""
Created on 2018-12-04 14:09:10
Last Modified on 2018-12-04 14:09:10

Functions to check file or directories.

@Author: Ying Huang
"""
import os
import re

def check_if_file_dup(ifile, suffix='', if_rmfile=False):
    """Check if file existed. If ture, create a new file, 
    and rename it by adding a number to the end of the 
    origin file name.
    
    ifile: <Path> path of a file that will be checked.
    suffix: <str> If provied, a number will be add to end of
        filename before the suffix, when renaming file.
        
    Example:

        # If files named "/home/test.csv, /home/test_1.csv" are
        # already exist, check_if_file_dup() will work as below:

        >>> check_if_file_dup("/home/test.csv", suffix='.csv')
        >>> '/home/test_2.csv'
    """
    assert isinstance(suffix, str)

    if if_rmfile is True:
        if os.path.exists(ifile):
            os.remove(ifile)
        return ifile

    regexp = r'(.+)' + suffix

    filename = ifile
    name_no_suffix = re.search(regexp, ifile).group(1)
    num = 1

    while os.path.exists(filename):

        filename = "{}_{}{}".format(name_no_suffix, str(num), suffix)

        num += 1

    return filename


