import requests 
import re

url = "https://www.ebi.ac.uk/chembl/api/data/molecule/search?q=pain killers"

r = requests.get(url, headers={ "Content-Type" : "application/json"})


results = re.findall(r'<molfile>([^<]*)?</molfile>', r.text)


with open("results.txt", "w") as file: 
    for result in results: 
        file.write(result)
        file.write("$$$$\n")