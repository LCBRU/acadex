from sqlalchemy import MetaData, Table, Column, Integer, NVARCHAR
from sqlalchemy.sql.schema import ForeignKey
from lbrc_flask.security.migrations import get_audit_mixin_columns

meta = MetaData()

def upgrade(migrate_engine):
    meta.bind = migrate_engine

    t = Table(
        "academic",
        meta,
        Column("id", Integer, primary_key=True),
        Column("google_scholar_id", NVARCHAR(255)),
        Column("name", NVARCHAR(500)),
        Column("affiliation", NVARCHAR(500)),
        Column("cited_by", Integer),
        Column("h_index", Integer),
        Column("i10_index", Integer),
        *get_audit_mixin_columns(),
    )

    t.create()


def downgrade(migrate_engine):
    meta.bind = migrate_engine
    t = Table("academic", meta, autoload=True)
    t.drop()