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
@click.option('-o', '--output-bed', help='Output file name. If nothing is specified, it will be genename.bed')
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

    try:
        input_json = read_json(input_json)
    except:
        LOG.error("Invalid input json string")
        raise click.Abort()

    if not output_bed:
        output_bed = input_json['genename']+".bed"
        LOG.info("Setting output file name as %s",output_bed) 
    else:
        LOG.debug("Output file is %s", output_bed)

    try:
        LOG.debug("Loading reference file into a dataframe")
        reference_df = pandas.read_csv(
            reference, delimiter='\t', encoding='utf-8')
    except:
        LOG.error("Couldn't load reference file into dataframe")
        raise click.Abort()

    try:
        LOG.debug("Extracting regions based on genename")
        df = copy.deepcopy(
            reference_df.loc[reference_df['name2'] == input_json['genename'], [
                'chrom', 'exonStarts', 'exonEnds', 'strand', 'name2', 'name',
                'exonCount'
            ]])
    except:
        LOG.error("Extracting genename from reference file failed")
        raise click.Abort()

    if df.shape[0] == 0:
        LOG.warning(
            "No enteries for input string was found in reference file.")

    try:
        LOG.debug("Converting dataframe to pybedtools object and expanding exonStarts and exonEnds")
        df_bed = pybedtools.BedTool.from_dataframe(df).expand(c="2,3").sort()
    except:
        LOG.error("Bed objet creation failed")
        raise click.Abort()

    try:
        LOG.debug("Merging bed regions and collapsing disctint strand, genename, transcript")
        df_bed.merge(c='4,5,6', o="distinct").saveas(output_bed)
    except:
        LOG.error("Merge and collapse of BED object failed")
        raise click.Abort()

if __name__ == '__main__':
    genepanel()
