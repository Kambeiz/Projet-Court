#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 16:20:13 2020

@author: kambei
"""

from selenium import webdriver
from selenium.webdriver.chrome.options import Options  

main_ch_main = []
main_ch_side = []
side_ch_side = []


url = "http://pic.mbu.iisc.ernet.in/job.html"
'''hbd_url = "http://pic.mbu.iisc.ernet.in/TEMP/"
pdb_name = "2dsq"

files= {'pdb_name': open("2dsq.pdb", "rb")}
values = {"hbond3":"func3", "hbond4":"func4", "hbond5":"func5"}
values_b= {'pdb_name': open("2dsq.pdb", "rb"), "hbond3":"func3", "hbond4":"func4", "hbond5":"func5"}

r2 = requests.post(url_2, files=files, data=values, allow_redirects=True)
get= requests.get("http://pic.mbu.iisc.ernet.in/TEMP/{}.hbd".format(pdb_name))

t = get.text
print(type(t))
v = t.split("\n")
v = [x for x in v if not x.startswith('#')]

'''
'''
session = HTMLSession()
r = session.get(url_2)
r.html.render()
#s = session.post(url_2, data=values_b, allow_redirects=True)
print(r.content)
'''

'''web = Browser()
web.go_to('http://caps.ncbs.res.in/iws/hbond.html')
up = web.find_elements(tag="input TYPE ='FILE'", number=1)
#up.send_keys("1BTA.pdb")

web.click('HBOND')
web.quit()
'''

chrome_options = Options()  
chrome_options.add_argument("--headless")  
#driver = webdriver.Chrome(options=chrome_options)
driver = webdriver.Chrome()

driver.get(url)

upload = driver.find_element_by_name("pdbname")
upload.send_keys("/home/kambei/Bureau/M2_Bioinfo/Projet-Court/data/1BTA.pdb")
main_chain_main = driver.find_element_by_name("hbond3")
main_chain_side = driver.find_element_by_name("hbond4")
side_chain_side = driver.find_element_by_name("hbond5")

#driver.get("http://caps.ncbs.res.in/iws/hbond.html")
#upload = driver.find_element_by_name("file")
#upload.send_keys("/home/kambei/Bureau/M2_Bioinfo/Projet-Court/data/1BTA.pdb")
#submit = driver.find_element_by_xpath("//input[@value='HBOND']")
main_chain_main.click()
main_chain_side.click()
side_chain_side.click()
submit = driver.find_element_by_name("Submit")
submit.click()
body = (driver.find_element_by_xpath("//body").text).split("\n")

count = 0
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
