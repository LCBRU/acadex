#/usr/bin/python

from pymed import PubMed
from pprint import pprint as pp
import lxml.etree as etree
from bs4 import BeautifulSoup

# Query Help
#
# https://pubmed.ncbi.nlm.nih.gov/advanced/
#
pubmed = PubMed(tool="Acadex", email="richard.a.bramley@uhl-tr.nhs.uk")
query = '(David Adlam[Author])'

results = pubmed.query(query, max_results=1)

print('Hello')
for r in results:
    # bs = BeautifulSoup(r.xml, 'xml')
    # print(bs.prettify())
    # pp(etree.tostring(r.xml.getroot(), pretty_print=True))
    # pp(r.toDict())
    pass
