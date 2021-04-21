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

    @property
    def pubmed_name(self):
        firstname, *_, lastname = self.name.split()
        return f'{lastname} {firstname[0]}'


academics_articles = db.Table(
    'academics_articles',
    db.Column(
        'academic_id',
        db.Integer(),
        db.ForeignKey('academic.id'),
        primary_key=True,
    ),
    db.Column(
        'article_id',
        db.Integer(),
        db.ForeignKey('article.id'),
        primary_key=True,
    ),
)


class Article(AuditMixin, CommonMixin, db.Model):

    id = db.Column(db.Integer(), primary_key=True)
    pm_id = db.Column(db.Integer())
    journal = db.Column(db.String(200))
    published_date = db.Column(db.Date)
    title = db.Column(db.UnicodeText())
    abstract = db.Column(db.UnicodeText())
    abstract_funding = db.Column(db.UnicodeText())

    authors = db.relationship(
        "Academic", secondary=academics_articles, collection_class=set, backref=db.backref("articles", lazy="joined")
    )
