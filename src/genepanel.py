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

def parseIntSet(nputstr=""):
    '''
    Expand range of numbers.
    https://stackoverflow.com/questions/712460/interpreting-number-ranges-in-python
    '''

    selection = set()
    invalid = set()
    # tokens are comma seperated values
    tokens = [x.strip() for x in nputstr.split(',')]
    for i in tokens:
        if len(i) > 0:
            if i[:1] == "<":
                i = "1-%s"%(i[1:])
        try:
            # typically tokens are plain old integers
            selection.add(int(i))
        except:
            # if not, then it might be a range
            try:
                token = [int(k.strip()) for k in i.split('-')]
                if len(token) > 1:
                    token.sort()
                    # we have items seperated by a dash
                    # try to build a valid range
                    first = token[0]
                    last = token[len(token)-1]
                    for x in range(first, last+1):
                        selection.add(x)
            except:
                # not an int and not a range...
                invalid.add(i)
    # Report invalid tokens before returning valid selection
    if len(invalid) > 0:
        LOG.error("Invalid set: " + str(invalid))
        return None
    else:
        return selection
# end parseIntSet


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

    if not input_json['genename']:
        LOG.error('Input string must have genename entry')
        raise click.Abort()

    # if there is exons entry, expand the range
    exon_range = False
    if 'exons' in input_json.keys():
        if input_json['exons']:
            exon_range = list(parseIntSet(input_json['exons']))
            if exon_range is None:
                LOG.error('Incorrect exon range. Example: <3,4-7,10-12,24')
                raise click.Abort()
            input_json['exons'] = ",".join(str(x) for x in exon_range)
            LOG.info('Exon ranges converted to %s', exon_range)

    # if there is no output file specified
    if not output_bed:
        output_bed = '_'.join(str(x) for x in input_json.values())
        for ch in ",:-":
            output_bed = output_bed.replace(ch,'_')
        output_bed += ".bed"
        LOG.info("Setting output file name as %s",output_bed) 
    else:
        LOG.debug("Output file is %s", output_bed)

    # read input reference file
    try:
        LOG.debug("Loading reference file into a dataframe")
        reference_df = pandas.read_csv(
            reference, delimiter='\t', encoding='utf-8')
        # create a new df by removing all characters after "." in transcript field
        # i.e. "name" column
        transcript_df = reference_df["name"].str.split(".", n = 1, expand = True)
        # replace transcript names with first value
        reference_df['name'] = transcript_df[0]
    except:
        LOG.error("Couldn't load reference file into dataframe")
        raise click.Abort()

    # filter reference file based on input json
    try:
        LOG.debug("Extracting regions based on genename")
        df = copy.deepcopy(
            reference_df.loc[reference_df['name2'] == input_json['genename'], [
                'chrom', 'exonStarts', 'exonEnds', 'strand', 'name2', 'name',
                'exonCount'
            ]])
        if 'transcript' in input_json.keys():
            if input_json['transcript']:
                df = df.loc[df['name'] == input_json['transcript']]
    except:
        LOG.error("Extracting genename from reference file failed")
        raise click.Abort()

    # if entries found in the reference, then exit
    if df.shape[0] == 0:
        LOG.warning(
            "No enteries for input string was found in reference file.")
        raise click.Abort()

    if 'exons' in input_json.keys():
        if df.name.count() > 1 or df.name2.count() > 1 or df.strand.count() > 1:
            LOG.error('More than one entry name (transcript) or name2 (gene) or strand column')
            LOG.error(df.head())
            raise click.Abort()

    # prepare filtered dataframe as bed file to write to output
    try:
        LOG.debug("Converting dataframe to pybedtools object and expanding exonStarts and exonEnds")
        df_bed = pybedtools.BedTool.from_dataframe(df).expand(c="2,3").sort().to_dataframe()
        df_out = copy.deepcopy(df_bed)
        if exon_range:
            exonNum = list(range(1, len(df_out)+1))

            if df.strand.unique()[0] == "-":
                exonNum.sort(reverse=True)

            df_out['thickStart'] = exonNum
            df_out = df_out[df_out['thickStart'].isin(exon_range)]
            df_out['thickStart'] = 'exon_num_' + df_out['thickStart'].astype(str)
        else:
            df_out['thickStart'] = 'total_exon_' + df_out['thickStart'].astype(str)
        
        df_bed = pybedtools.BedTool.from_dataframe(df_out).sort()
    except:
        LOG.error("Bed object creation failed")
        raise click.Abort()

    # collapse columns 
    try:
        LOG.debug("Merging bed regions and collapsing disctint strand, genename, transcript")
        df_bed.merge(c='4,5,6,7', o="distinct").saveas(output_bed)
    except:
        LOG.error("Merge and collapse of BED object failed")
        raise click.Abort()

if __name__ == '__main__':
    genepanel()
