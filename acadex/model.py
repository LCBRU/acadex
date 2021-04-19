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
