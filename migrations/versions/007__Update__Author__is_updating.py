from sqlalchemy import MetaData, Table, Column, Boolean
from sqlalchemy.sql.sqltypes import Boolean

meta = MetaData()

def upgrade(migrate_engine):
    meta.bind = migrate_engine

    t = Table("academic", meta, autoload=True)
    column = Column("is_updating", Boolean)
    column.create(t)


def downgrade(migrate_engine):
    meta.bind = migrate_engine
    t = Table("academic", meta, autoload=True)
    t.c.is_updating.drop()
