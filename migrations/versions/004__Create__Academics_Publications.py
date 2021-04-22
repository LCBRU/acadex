from sqlalchemy import (
    MetaData,
    Table,
    Column,
    Integer,
    ForeignKey,
)

meta = MetaData()


def upgrade(migrate_engine):
    meta.bind = migrate_engine

    academic = Table("academic", meta, autoload=True)
    publication = Table("publication", meta, autoload=True)

    t = Table(
        "academics_publications",
        meta,
        Column("id", Integer, primary_key=True),
        Column("academic_id", Integer, ForeignKey(academic.c.id), index=True, nullable=False),
        Column("publication_id", Integer, ForeignKey(publication.c.id), index=True, nullable=False),
    )
    t.create()


def downgrade(migrate_engine):
    meta.bind = migrate_engine
    t = Table("academics_publications", meta, autoload=True)
    t.drop()
