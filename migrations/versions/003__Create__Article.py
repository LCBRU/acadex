from sqlalchemy import MetaData, Table, Column, Integer, NVARCHAR, UnicodeText, Date
from sqlalchemy.sql.schema import ForeignKey
from lbrc_flask.security.migrations import get_audit_mixin_columns

meta = MetaData()

def upgrade(migrate_engine):
    meta.bind = migrate_engine

    a = Table("academic", meta, autoload=True)

    t = Table(
        "article",
        meta,
        Column("id", Integer, primary_key=True),
        Column("pm_id", Integer, unique=True),
        Column("journal", UnicodeText),
        Column("published_date", Date),
        Column("title", UnicodeText),
        Column("abstract", UnicodeText),
        Column("abstract_funding", UnicodeText),
        *get_audit_mixin_columns(),
    )

    t.create()


def downgrade(migrate_engine):
    meta.bind = migrate_engine
    t = Table("article", meta, autoload=True)
    t.drop()