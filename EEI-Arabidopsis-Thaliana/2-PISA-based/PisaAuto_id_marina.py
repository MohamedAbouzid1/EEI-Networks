#!/usr/bin/python3
"""
Code to automatically run PDBePISA web server on the given pdb id and downloading the
generated xml files.

  How to use
  ----------
First you need to have the python packages selenium, halo and argparse installed.

Then you can run the script with the following command :

    python PisaAuto_id.py pdb_id store_dir

Note that right now it's made for firefox browser but adding other browsers 
isn't hard to implement (ex: driver = webdriver.Chrome() for chrome).
Also note that only the interface table, the residues interaction and interfacing 
residues xml files are downloaded.

  Author
  ------
    Hocine Meraouna

"""

import time
import sys
import os
import argparse
import re
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.firefox.service import Service
from selenium.common.exceptions import NoSuchElementException
from halo import Halo


def check_exists_by_name(name, driver):
    """
    The function to check if an element is present on the webdriver.

    Parameters
    ----------
    driver : selenium webdriver
    name : string
        the name of the element

    Returns
    -------
    boolean
    """
    try:
        driver.find_element('name', name)
    except NoSuchElementException:
        return False
    return True


def start():
    """
    The function to access to the pisa web server.
    I'm using firefox but it can be changed for other browsers.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    selenium webdriver
    """
    # Remove or comment out print statements unless they're critical
    # print("1- Accessing to PISA website :")
    
    # Set environment variables for Firefox to find the GTK library
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        # Save the original LD_LIBRARY_PATH to restore it later if needed
        old_ld_library_path = os.environ.get('LD_LIBRARY_PATH', '')
        
        # Set LD_LIBRARY_PATH to include conda environment libraries
        if old_ld_library_path:
            os.environ['LD_LIBRARY_PATH'] = f"{conda_prefix}/lib:{old_ld_library_path}"
        else:
            os.environ['LD_LIBRARY_PATH'] = f"{conda_prefix}/lib"
    
    # Set temporary directory
    os.environ["TMPDIR"] = "/home/bbf3630/tmp"

    # Setup Firefox options
    options = webdriver.FirefoxOptions()
    options.set_preference("browser.cache.disk.parent_directory", "/home/bbf3630/tmp")
    options.add_argument("--headless")
    
    # Setup the service with the geckodriver path
    service = Service("/home/bbf3630/miniconda3/envs/selenium_env/bin/geckodriver")
    
    # Create the driver
    driver = webdriver.Firefox(service=service, options=options)
    driver.get("https://www.ebi.ac.uk/pdbe/pisa/")
    launch = driver.find_element("name", "start_server")
    launch.click()
    return driver


def launch_pdb_id(driver, pdb_id):
    """
    The function to run pisa web service on the pdb id.
    
    Parameters
    ----------
    driver : selenium webdriver
        given by the function start()
    pdb_id : string
        pdb id given by the user
    
    Returns
    -------
    selenium webdriver
    """
    # Only keep essential status messages
    print(f"Processing {pdb_id}", file=sys.stderr)
    
    pdb_entry = driver.find_element("name", "edt_pdbcode")
    pdb_entry.clear()
    pdb_entry.send_keys(pdb_id)
    time.sleep(10)
    
    interface = driver.find_element("name", "btn_submit_interfaces")
    interface.click()
    
    while(not check_exists_by_name('downloadXML', driver)):
        pass
    
    time.sleep(2)
    return driver


def download_xmls(driver, pdb_id, store_dir):
    """
    The function to download the xml files.

    Parameters
    ----------
    driver : selenium webdriver
    pdb_id : string
    store_dir : string
        directory to store the results
    
    Returns
    -------
    Nothing
    """
    # Remove spinner and most print statements
    driver.find_element('name', 'downloadXML').click()
    time.sleep(5)

    driver.switch_to.window(driver.window_handles[1])
    xml = driver.current_url

    # Only log errors to stderr
    if not os.path.exists(store_dir):
        print(f"Creating directory {store_dir}", file=sys.stderr)
        os.makedirs(store_dir)

    if not os.path.exists(store_dir + '/' + pdb_id):
        os.makedirs(store_dir + '/' + pdb_id)

    with open(store_dir + '/' + pdb_id + '/' + xml.split('/')[-1], 'w') as f:
        f.write(driver.page_source)

    time.sleep(3)

    driver.close()

    time.sleep(2)

    driver.switch_to.window(driver.window_handles[0])

    inter_lst = []

    # Fix for interfaces > 99
    with open(store_dir + '/' + pdb_id + '/' + xml.split('/')[-1], "r") as f_xml:
        for line in f_xml:
            if line.strip().startswith("<INTERFACENO>"):
                 # Use regex to extract the interface number properly
                match = re.search(r"<INTERFACENO>(\d+)</INTERFACENO>", line)
                if match:
                    inter_lst.append(match.group(1))


    # Remove spinner and use minimal progress indication
    for i in inter_lst:
        print(f"Interface {i}/{len(inter_lst)}", file=sys.stderr)

        time.sleep(2)

        driver.find_element('link text', i).click()

        time.sleep(5)

        xmls = driver.find_elements('name', 'downloadXML')

        for j in range(1, len(xmls)):
            xmls[j].click()

            time.sleep(3)

            driver.switch_to.window(driver.window_handles[1])
            xml = driver.current_url

            with open(store_dir + '/' + pdb_id + '/' + xml.split('/')[-1], 'w') as f:
                f.write(driver.page_source)

            time.sleep(3)

            driver.close()

            time.sleep(2)

            driver.switch_to.window(driver.window_handles[0])

        driver.back()

        time.sleep(3)


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser()

    PARSER.add_argument("pdb_id", help="the id of the pdb you want to run pisa on", type=str)
    PARSER.add_argument("store_dir", help="the folder to store results", type=str)

    ARGS = PARSER.parse_args()

    PDB_ID = ARGS.pdb_id
    STR_DIR = ARGS.store_dir

    download_xmls(launch_pdb_id(start(), PDB_ID), PDB_ID, STR_DIR)