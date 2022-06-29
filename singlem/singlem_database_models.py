"""
models.py
- SQLAlchemy Data classes for SingleM databases
"""

from sqlalchemy import Column, Integer, String, ForeignKey, Float
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

# declarative base class
Base = declarative_base()

class Otu(Base):
    '''
    sqlite> select * from otus limit 3;
    id|sample_name|num_hits|coverage|taxonomy_id|marker_id|sequence_id|sequence
    1|GB_GCA_000309865.1_protein|1|1.03|0|16|351|GCCGACCCCAATATCATCGCTGATCTGGACTCCCATCATCTACTATTCAAAGAAGGCATC
    2|GB_GCA_000309865.1_protein|1|1.03|7|46|1087|ACCAGTAAGAACTGGGTGATCTGGGCAGCTGACTTTATGGAGAAATTTGATGCGGATCTG
    3|GB_GCA_000309865.1_protein|1|1.03|8|50|1140|CGCTGGGAAGCTGGTGGAGCC------------AAAGGCCTGGATCGCGTGCATGAATTC
    '''
    __tablename__ = 'otus'
    id = Column(Integer, primary_key=True)
    # sample_name|num_hits|coverage|taxonomy|marker_id|sequence_id
    sample_name = Column(String, nullable=False, index=True)
    num_hits = Column(Integer, nullable=False)
    coverage = Column(Float, nullable=False)
    taxonomy_id = Column(Integer, ForeignKey('taxonomy.id'), nullable=False, index=True)
    marker_id = Column(Integer, ForeignKey('markers.id'), nullable=False, index=True)
    sequence_id = Column(Integer, ForeignKey('nucleotides.id'), nullable=False, index=True)
    # We include the sequence itself so dumping the database is faster
    sequence = Column(String, nullable=False)

class Taxonomy(Base):
    __tablename__ = 'taxonomy'
    id = Column(Integer, primary_key=True)
    taxonomy = Column(String, nullable=False, index=True)

class NucleotideSequence(Base):
    '''
    id|marker_id|sequence|marker_wise_id
    1|1|GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC|0
    2|2|GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT|0
    3|3|CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC|0
    '''
    __tablename__ = 'nucleotides'
    id = Column(Integer, primary_key=True)
    marker_id = Column(Integer, ForeignKey('markers.id'), nullable=False, index=True)
    sequence = Column(String, nullable=False, index=True)
    marker_wise_id = Column(Integer, nullable=False, index=True)

    # nucleotide_proteins = relationship('NucleotidesProteins', back_populates='nucleotide_sequence.nucleotide_id')

class Marker(Base):
    __tablename__ = 'markers'
    id = Column(Integer, primary_key=True)
    marker = Column(String, nullable=False, index=True)
    otus = relationship('Otu', backref='marker')
    nucleotides = relationship('NucleotideSequence', backref='marker')

class ProteinSequence(Base):
    '''
    id|marker_wise_id|protein_sequence
    1|0|GKANPAPPVGPALGQAGVNI
    2|0|AKLGDIVKIQETRPLSATKR
    3|0|RRWNPKMKKYIFTERNGIYI
    '''
    __tablename__ = 'proteins'
    id = Column(Integer, primary_key=True)
    marker_id = Column(Integer, ForeignKey('markers.id'), nullable=False, index=True)
    protein_sequence = Column(String, nullable=False, index=True)

    # nucleotide_proteins = relationship('NucleotidesProteins', back_populates='nucleotide_sequence.protein_id')

class NucleotidesProteins(Base):
    __tablename__ = 'nucleotides_proteins'
    id = Column(Integer, primary_key=True)
    nucleotide_id = Column(Integer, ForeignKey('nucleotides.id'), nullable=False, index=True)
    protein_id = Column(Integer, ForeignKey('proteins.id'), nullable=False, index=True)