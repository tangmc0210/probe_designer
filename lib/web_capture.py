import requests
from bs4 import BeautifulSoup


url = 'https://mouse.brain-map.org/experiment/show/75081210'
response = requests.get(url)
html = response.text

soup = BeautifulSoup(html, 'html.parser')
probe_sequence_element = soup.find(text='Sequence:').parent
probe_sequence = probe_sequence_element.next_sibling.strip()

print(probe_sequence)