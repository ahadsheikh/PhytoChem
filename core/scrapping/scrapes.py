from bs4 import BeautifulSoup
import requests


def single_smiles_scrape(smiles: str) -> str:
    # Give SwissADME data from  a smiles

    url = "http://www.swissadme.ch/index.php"
    s = requests.post(url, data={'smiles': smiles})

    soup = BeautifulSoup(s.content, 'html.parser')
    csv_div = soup.find('div', attrs={
        "style": "float: right; width: 250px; height: 30px; text-align: left; font: 13pt  Helvetica,Verdana, "
                 "sans-serif;"})
    csv_link = "http://www.swissadme.ch/" + csv_div.a['href']

    csv_data = requests.get(csv_link)

    return csv_data.text