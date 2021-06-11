from lbrc_flask.database import db
from lbrc_flask.security import AuditMixin
from lbrc_flask.model import CommonMixin


class Academic(AuditMixin, CommonMixin, db.Model):

    id = db.Column(db.Integer(), primary_key=True)
    google_scholar_id = db.Column(db.String(255))
    name = db.Column(db.String(500))
    affiliation = db.Column(db.String(500))
    cited_by = db.Column(db.Integer)
    h_index = db.Column(db.Integer)
    i10_index = db.Column(db.Integer)
    is_updating = db.Column(db.Boolean)

    @property
    def pubmed_name(self):
        firstname, *_, lastname = self.name.split()
        return f'{lastname} {firstname[0]}'


academics_publications = db.Table(
    'academics_publications',
    db.Column(
        'academic_id',
        db.Integer(),
        db.ForeignKey('academic.id'),
        primary_key=True,
    ),
    db.Column(
        'publication_id',
        db.Integer(),
        db.ForeignKey('publication.id'),
        primary_key=True,
    ),
)


class Publication(AuditMixin, CommonMixin, db.Model):

    id = db.Column(db.Integer(), primary_key=True)
    pm_id = db.Column(db.Integer())
    journal = db.Column(db.String(200))
    published_date = db.Column(db.Date)
    title = db.Column(db.UnicodeText())

    academics = db.relationship(
        "Academic", secondary=academics_publications, collection_class=set, backref=db.backref("publications", lazy="joined")
    )


class AbstractSection(AuditMixin, CommonMixin, db.Model):

    id = db.Column(db.Integer(), primary_key=True)
    publication_id = db.Column(db.Integer(), db.ForeignKey(Publication.id))
    publication = db.relationship(Publication, lazy="joined", backref='abstracts')
    label = db.Column(db.String(200))
    text = db.Column(db.UnicodeText())


class Author(AuditMixin, CommonMixin, db.Model):

    id = db.Column(db.Integer(), primary_key=True)
    publication_id = db.Column(db.Integer(), db.ForeignKey(Publication.id))
    publication = db.relationship(Publication, lazy="joined", backref='authors')
    last_name = db.Column(db.String(100))
    fore_name = db.Column(db.String(100))
    initials = db.Column(db.String(100))
    affiliation = db.Column(db.UnicodeText())

    @property
    def full_name(self):
        return f'{self.fore_name} {self.last_name}'

