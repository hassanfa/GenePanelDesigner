#!/usr/bin/env python
"""
# Copyright (c) Hassan Foroughi <hassan.foroughi@scilifelab.s>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
#
# A simple script to extract coordinates for input genes, transcriptts, exons, etc.
# This is to facilitate gene panel design.
"""
import logging
import os
import click
import coloredlogs
import pandas
import copy
import pybedtools

from genepanel_utils import read_json

LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
LOG = logging.getLogger(__name__)

__version__ = 0.01

@click.command()
@click.option(
    '--log-level',
    default='DEBUG',
    type=click.Choice(LOG_LEVELS),
    help='Set the level of log output.',
    show_default=True)
@click.option(
    '-r',
    '--reference',
    type=click.Path(),
    required=True,
    help='Reference refFlat file. See ../data/ncbiRefSeq.txt as an example.')
@click.option(
    '-i',
    '--input-json',
    required=True,
    help='''
              Input string json with the following keys:
              genename, transcript, exon_number, coordinate.
              ''')
@click.option('-o', '--output-bed', required=True, help='Output file name.')
@click.option(
    '--strand-match/--no-strand-match',
    default=False,
    show_default=True,
    help='''
              Consider strandness when counting exon number.
              e.g. if strand is negative, the exon with largest coordinates
              is exon 1.
              ''')
def genepanel(log_level, reference, input_json, output_bed, strand_match):
    """
    Read input json string of genename, transcript, exon_number, and coordinate
    and write a 4 column bedfile.
    """
    coloredlogs.install(level=log_level)

    LOG.info("Running version %s", __version__)
    LOG.debug("Debug logging enabled.")
    LOG.debug("Input string %s", input_json)
    LOG.debug("Reference file is %s", reference)
    LOG.debug("Strandness search is set to %s", strand_match)
    LOG.debug("Output file is %s", output_bed)

    try:
        input_json = read_json(input_json)
    except:
        LOG.error("Invalid input json string")
        raise click.Abort()

    try:
        LOG.debug("Loading reference file into a dataframe")
        reference_df = pandas.read_csv(reference,delimiter='\t',encoding='utf-8')
    except:
        LOG.error("Couldn't load reference file into dataframe")
        raise click.Abort()

    df = copy.deepcopy(reference_df.loc[reference_df['name2']==input_json['genename'],
        ['chrom','exonStarts', 'exonEnds','strand', 'name2','name','exonCount']])
    df['bed_4_name'] = df[['name2','name','exonCount']].apply(lambda x: '|'.join(x.values.astype(str)), axis=1)
    df = df.drop(columns=['name2', 'name', 'exonCount'])
    df_bed=pybedtools.BedTool.from_dataframe(df).sort()
    df_bed.merge().saveas(output_bed)


if __name__ == '__main__':
    genepanel()
