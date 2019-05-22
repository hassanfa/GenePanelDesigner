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
import re
import pybedtools

from genepanel_utils import read_json
from genepanel_utils import get_exon_range
from genepanel_utils import build_reference
from genepanel_utils import build_output_file
from genepanel_utils import filter_gene
from genepanel_utils import filter_exon
from genepanel_utils import filter_transcript

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
    help='''
              Input string json with the following keys:
              genename, transcript, exon_number, coordinate.
              ''')
@click.option(
    '-o',
    '--output-bed',
    help='Output file name. If nothing is specified, it will be genename.bed')
@click.option(
    '--strand-match/--no-strand-match',
    default=False,
    show_default=True,
    help='''
              Consider strandness when counting exon number.
              e.g. if strand is negative, the exon with largest coordinates
              is exon 1.
              ''')
@click.option(
    '-g',
    '--gene',
    help="Genename/gene symbol, Accepted examples: TP53 (only one gene)")
@click.option(
    '-t',
    '--transcript',
    help="Transcript. Accepted examples: NM_123545 (only one transcript)")
@click.option(
    '-e',
    '--exon',
    help="Exon numbers. Accepted examples: '<3', '1-6', '1,4,10'")
@click.option(
    '-c',
    '--coordinate',
    help="Coordinate. Accepted examples: '16:123456-123457,18:123456-123457'")
@click.option(
    '-d', '--directory', help="Output directory/path for the bedfiles")
@click.option(
    '-p', '--padding', help="Padding in bp.", default=30, show_default=True)
def genepanel(log_level, reference, input_json, output_bed, strand_match, gene,
              transcript, exon, coordinate, directory, padding):
    """
    Read input json string of genename, transcript, exon_number, and coordinate
    and write a 4 column bedfile.
    """
    coloredlogs.install(level=log_level)

    LOG.info("Running version %s", __version__)
    LOG.debug("Debug logging enabled.")
    LOG.debug("Reference file is %s", reference)
    LOG.debug("Strandness search is set to %s", strand_match)

    if (input_json in locals()) ^ (gene in locals()):
        LOG.error("Only one or either of input_json or gene is required")
        raise click.Abort()

    if input_json:
        input_json = read_json(input_json)

    if gene:
        input_json = dict()
        input_json['genename'] = gene

    if transcript:
        input_json['transcript'] = transcript

    if exon:
        input_json['exons'] = exon

    if coordinate:
        input_json['coordinate'] = coordinate

    for key in list(input_json.keys()):
        for ch in '" ':
            input_json[key] = input_json[key].replace(ch, '')
        if input_json[key] == "":
            input_json.pop(key)

    LOG.info("Input string %s", input_json)

    if not input_json['genename']:
        LOG.error('Input string must have genename entry')
        raise click.Abort()

    exon_range = get_exon_range(input_json)

    output_bed = build_output_file(
        output_bed=output_bed, input_json=input_json, directory=directory)

    reference_df = build_reference(reference)

    df = filter_gene(reference_df, input_json)
    if not 'coordinate' in input_json.keys():
        df = filter_transcript(df, input_json)

        if exon_range:
            if df.name.count() > 1 or df.name2.count() > 1 or df.strand.count(
            ) > 1:
                LOG.error(
                    'More than one entry name (transcript) or name2 (gene) or strand column'
                )
                LOG.error(df.head())
                raise click.Abort()

        # prepare filtered dataframe as bed file to write to output
#        try:
        if True:
            LOG.debug(
                "Converting dataframe to pybedtools object and expanding exonStarts and exonEnds"
            )
            df = filter_exon(df, exon_range)
            df_bed = pybedtools.BedTool.from_dataframe(df).sort()


#        except:
#            LOG.error("Bed object creation failed")
#            raise click.Abort()
    else:
        bed_ranges = input_json['coordinate'].split(',')
        bed_start = ""
        bed_end = ""
        df_bed = copy.deepcopy(df[0:0])
        df_bed['name2'] = df.name2.unique()
        df_bed['chrom'] = df.chrom.unique()
        df_bed['strand'] = df.strand.unique()
        for bed in bed_ranges:
            bed = re.split(':|-', bed)
            chrom = 'chr' + str(bed[0])
            if not chrom == df_bed['chrom'][0]:
                LOG.error('chormosome not matching refernce')
                LOG.error(df)
                LOG.error(chrom)
                raise click.Abort()
            bed_start = bed[1] + "," + bed_start
            bed_end = bed[2] + "," + bed_end
        df_bed['exonStarts'] = bed_start
        df_bed['exonEnds'] = bed_end
        df_bed['exonCount'] = "1"
        df_bed['name'] = "."

        df_bed = df_bed[[
            'chrom', 'exonStarts', 'exonEnds', 'name2', 'exonCount', 'strand'
        ]]
        df_bed = pybedtools.BedTool.from_dataframe(df_bed).expand(
            c="2,3").sort()

    # collapse columns
    try:
        LOG.debug(
            "Merging bed regions and collapsing disctint strand, genename, transcript"
        )

        df_bed.slop(
            b=padding, g='hg19.genome').merge(
                c='4,5,6', o="distinct").saveas(output_bed)
    except:
        LOG.error("Merge and collapse of BED object failed")
        LOG.error(df_bed)
        LOG.error(output_bed)
        raise click.Abort()

if __name__ == '__main__':
    genepanel()
