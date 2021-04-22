#/usr/bin/python

from Bio import Entrez
from pprint import pprint as pp

Entrez.email = "richard.a.bramley@uhl-tr.nhs.uk"

def iterdict(d, lvl=0):
    for k,v in d.items():        
        print (' ' * lvl * 4, k)
        if isinstance(v, dict):
            iterdict(v, lvl + 1)


# handle = Entrez.esearch(db="pubmed", retstart=5, retmax=10, term="(N Samani[Author])")
handle = Entrez.esearch(db="pubmed", retstart=5, retmax=10, term="The Monomorphic Clear Cell Tumor: A Report Of Two Cases")
records = Entrez.read(handle)
print(records['Count'])
print(records['RetMax'])
print(records)
handle.close()
# print(records['IdList'])

handle = Entrez.efetch(db="pubmed", id=records['IdList'], retmode='xml')

from dataclasses import dataclass

@dataclass
class Abstract:
    label: str
    text: str

    def summary(self):
        if len(self.label) > 0:
            return f'{self.label}: {self.text}'
        else:
            return self.text
    
    
print(handle)

records = Entrez.read(handle)
for r in records['PubmedArticle']:
    art = r['MedlineCitation']['Article']
    print(art)

    # for au in art['AuthorList']:
    #     print(' '.join([aff.get('Affiliation', None) for aff in au.get("AffiliationInfo", [])]))

    # print(art['ArticleDate'][0]['Year'])

    # abstracts = []

    # if 'Abstract' in art:
    #     for at in [v for k, v in art['Abstract'].items() if k == 'AbstractText']:
    #         for s in at:
    #             abstracts.append(Abstract(label=s.attributes.get("Label", ""), text=s))

    # print('\n'.join(a.summary() for a in abstracts))

handle.close()
