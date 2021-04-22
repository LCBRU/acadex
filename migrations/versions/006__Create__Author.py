from sqlalchemy import MetaData, Table, Column, Integer, UnicodeText, Date, NVARCHAR
from sqlalchemy.sql.schema import ForeignKey
from lbrc_flask.security.migrations import get_audit_mixin_columns

meta = MetaData()

def upgrade(migrate_engine):
    meta.bind = migrate_engine

    p = Table("publication", meta, autoload=True)

    t = Table(
        "author",
        meta,
        Column("id", Integer, primary_key=True),
        Column("publication_id", Integer, ForeignKey(p.c.id), index=True, nullable=False),
        Column("last_name", NVARCHAR(100)),
        Column("fore_name", NVARCHAR(100)),
        Column("initials", NVARCHAR(100)),
        Column("affiliation", UnicodeText),
        *get_audit_mixin_columns(),
    )

    t.create()


def downgrade(migrate_engine):
    meta.bind = migrate_engine
    t = Table("author", meta, autoload=True)
    t.drop()