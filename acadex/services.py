from datetime import date, datetime
from acadex.model import AbstractSection, Academic, Author, Publication
from flask_security import current_user
from lbrc_flask.celery import celery
from lbrc_flask.database import db
from lbrc_flask.validators import parse_date_or_none
from lbrc_flask.security.model import User
from flask import current_app
from scholarly import scholarly
from Bio import Entrez


def add_or_update_academic(google_scholar_id):
    _add_or_update_academic.delay(google_scholar_id, current_user.id)


@celery.task()
def _add_or_update_academic(google_scholar_id, user_id):
    user = User.query.get(user_id)

    current_app.logger.info(f'Adding or Updating Academic: {google_scholar_id} as user {user}')

    Entrez.email = user.email

    resp = scholarly.fill(scholarly.search_author_id(google_scholar_id), sections=['indices'])

    if resp:
        a = Academic.query.filter(Academic.google_scholar_id == google_scholar_id).one_or_none()

        if a is None:
            a = Academic(google_scholar_id=google_scholar_id)
        
        a.name=resp['name']
        a.affiliation=resp['affiliation']
        a.cited_by=resp['citedby']
        a.h_index=resp['hindex']
        a.i10_index=resp['i10index']
        a.last_update_date=datetime.utcnow()

        db.session.add(a)
        db.session.commit()

        _update_publications(a)


    current_app.logger.info(f'Adding or Updating Academic Completed: {google_scholar_id}')


def _update_publications(academic):
    current_app.logger.info('Updating Publications')

    retstart = 0
    count = 1
    batch_size = 100

    while retstart < count:
        handle = Entrez.esearch(db="pubmed", retstart=retstart, retmax=batch_size, term=f"({academic.pubmed_name}[Author]) ")
        pubmed_records = Entrez.read(handle)
        handle.close()

        count = int(pubmed_records['Count'])
        retstart += batch_size

        handle = Entrez.efetch(db="pubmed", id=pubmed_records['IdList'], retmode='xml')

        pubmed_records = Entrez.read(handle)
        for pubmed_record in pubmed_records['PubmedArticle']:
            pm_id = int(pubmed_record['MedlineCitation']['PMID'])

            publication = Publication.query.filter(Publication.pm_id == pm_id).one_or_none()

            if publication is None:
                publication = Publication(pm_id=pm_id)

            _update_publication(pubmed_record, publication)

            publication.academics.add(academic)

            db.session.add(publication)
            db.session.commit()

        handle.close()

    current_app.logger.info('Updating Publications Completed')


def _update_publication(pubmed_record, publication):
    art = pubmed_record['MedlineCitation']['Article']
    publication.journal = art['Journal']['Title']

    art['Journal']['Title']

    if len(art['ArticleDate']) > 0:
        artdate = art['ArticleDate'][0]
        publication.published_date = date(
                    int(artdate['Year']),
                    int(artdate['Month']),
                    int(artdate['Day']),
                )
    else:
        pub_date = art['Journal']['JournalIssue']['PubDate']
        if 'MedlineDate' in pub_date.keys():
            publication.published_date = parse_date_or_none(pub_date["MedlineDate"])
        else:
            publication.published_date = parse_date_or_none(
                f'{pub_date.get("Day", 1)} {pub_date.get("Month", 1)} {pub_date["Year"]}'
            )

    publication.title = art['ArticleTitle']

    if publication.id:
        AbstractSection.query.filter(AbstractSection.publication_id == publication.id).delete()
        Author.query.filter(Author.publication_id == publication.id).delete()

    if 'Abstract' in art:
        for at in [v for k, v in art['Abstract'].items() if k == 'AbstractText']:
            for s in at:
                db.session.add(AbstractSection(
                    publication=publication,
                    label=s.attributes.get("Label", None),
                    text=s,
                ))

    for au in art['AuthorList']:
        db.session.add(Author(
            publication=publication,
            last_name=au.get("LastName", None),
            fore_name=au.get("ForeName", None),
            initials=au.get("Initials", None),
            affiliation=' '.join([aff.get('Affiliation', None) for aff in au.get("AffiliationInfo", [])]),
        ))
