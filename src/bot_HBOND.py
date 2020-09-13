#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 16:20:13 2020

@author: Debbah Nagi
"""

import os
import pandas
from selenium import webdriver
from selenium.webdriver.chrome.options import Options  

url = "http://pic.mbu.iisc.ernet.in/job.html"
list_elems = ["hbond3", "hbond4", "hbond5", "Submit"]


def launching_HBONDS(url):

    chrome_options = Options()  
    chrome_options.add_argument("--headless")  
    driver = webdriver.Chrome(options=chrome_options)
    #driver = webdriver.Chrome()
    
    driver.get(url)
    
    upload = driver.find_element_by_name("pdbname")

    #driver.get("http://caps.ncbs.res.in/iws/hbond.html")
    #upload = driver.find_element_by_name("file")
    #upload.send_keys("/home/kambei/Bureau/M2_Bioinfo/Projet-Court/data/1BTA.pdb")
    #submit = driver.find_element_by_xpath("//input[@value='HBOND']")
    return driver, upload

def find_and_click(list_elements,file, upload, driver):
    absolute_path = os.path.abspath(file)
    upload.send_keys(absolute_path)

    for elem in list_elements:
        click_elem = driver.find_element_by_name(elem)
        click_elem.click()
    return driver 

def body_to_list(driver):
    count = 0
    main_ch_main = []
    main_ch_side = []
    side_ch_side = []
    body = (driver.find_element_by_xpath("//body").text).split("\n")
    driver.close()
    
    for line in body:
        if "Main Chain-Main" in line:
            main_main = True
        elif "Main Chain-Side" in line:
            main_side = True
        elif "Side Chain-Side" in line:
            side_side = True
        if len(line) > 9:
            if line[:1].isdigit() and main_main == True:
                main_ch_main.append(line)
            elif line.startswith("Dd-a") and count ==0:
                main_main = False
                count += 1
            elif line[:1].isdigit() and main_side == True:
                main_ch_side.append(line)
            elif line.startswith("Dd-a") and count == 1:
                main_side = False
                count += 1
            elif line[:1].isdigit() and side_side == True:
                side_ch_side.append(line)
            elif line.startswith("Dd-a") and count == 2:
                side_side=False
    return main_ch_main, main_ch_side, side_ch_side    

def list_to_df(list_interactions):
    df = pandas.DataFrame(list_interactions)
    df = df[0].str.split(pat=" ", expand=True)
    return df

# Main program
if __name__ == "__main__":
    file = "1BTA.pdb"

    page_hbonds, upload = launching_HBONDS(url)

    page_results = find_and_click(list_elems,file, upload, page_hbonds)
    main_ch_main, main_ch_side, side_ch_side = body_to_list(page_results)
    df_main_main, df_main_side, df_side_side = list_to_df(main_ch_main), list_to_df(main_ch_side), list_to_df(side_ch_side)
