#/usr/bin/python
import json, sys
from scholarly import scholarly
from pprint import pprint as pp

search_query = scholarly.search_author('adlam, leicester')
author = scholarly.fill(next(search_query), sections=['indices', 'counts', 'publications'])
p = author['publications'][0]
pp(p)
p = scholarly.fill(p)
pp(p)

# print(scholarly.bibtex(p))
# author = scholarly.search_author_id('the id')
# author = author.fill(sections=['indices'])

# if author:
#     _dict = {
#         'id': author.id, 
#         'name': author.name,
#         'affiliation': author.affiliation,
#         'hindex': author.hindex,
#         'hindex5y': author.hindex5y,
#         'i10index': author.i10index,
#         'i10index5y': author.i10index5y,
#     }
#     print(json.dumps(_dict))
