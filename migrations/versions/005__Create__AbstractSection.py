from sqlalchemy import MetaData, Table, Column, Integer, UnicodeText, Date, NVARCHAR
from sqlalchemy.sql.schema import ForeignKey
from lbrc_flask.security.migrations import get_audit_mixin_columns

meta = MetaData()

def upgrade(migrate_engine):
    meta.bind = migrate_engine

    p = Table("publication", meta, autoload=True)

    t = Table(
        "abstract_section",
        meta,
        Column("id", Integer, primary_key=True),
        Column("publication_id", Integer, ForeignKey(p.c.id), index=True, nullable=False),
        Column("label", NVARCHAR(200)),
        Column("text", UnicodeText),
        *get_audit_mixin_columns(),
    )

    t.create()


def downgrade(migrate_engine):
    meta.bind = migrate_engine
    t = Table("abstract_section", meta, autoload=True)
    t.drop()