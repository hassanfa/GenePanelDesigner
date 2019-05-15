import json
import copy
import logging
import pandas
import click
import pybedtools

LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
LOG = logging.getLogger(__name__)


def read_json(input_string):
    '''
    Reads a json string and outputs a dictionary
    '''

    try:
        input_json = json.loads(input_string)
    except:
        LOG.error("Invalid input json string")
        raise click.Abort()

    # remove empty keys
    for key in list(input_json.keys()):
        if not input_json[key]:
            LOG.warning("Removed %s with empty value", key)
            input_json.pop(key)

    return input_json


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
                i = "1-%s" % (i[1:])
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
                    last = token[len(token) - 1]
                    for x in range(first, last + 1):
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


def get_exon_range(input_json):
    '''
    prepare an exon range list
    '''
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

    return exon_range


def build_output_file(output_bed, input_json):
    '''
    Build output filename
    '''
    # if there is no output file specified
    if not output_bed:
        output_bed = '_'.join(str(x) for x in input_json.values())
        for ch in ",:-":
            output_bed = output_bed.replace(ch, '_')
        output_bed += ".bed"
        LOG.info("Setting output file name as %s", output_bed)
    else:
        LOG.debug("Output file is %s", output_bed)

    return output_bed


def build_reference(reference):
    '''
    Build reference dataframe from reference file
    '''

    # read input reference file
    try:
        LOG.debug("Loading reference file into a dataframe")
        reference_df = pandas.read_csv(
            reference, delimiter='\t', encoding='utf-8')
        # create a new df by removing all characters after "." in transcript field
        # i.e. "name" column
        transcript_df = reference_df["name"].str.split(".", n=1, expand=True)
        # replace transcript names with first value
        reference_df['name'] = transcript_df[0]
    except:
        LOG.error("Couldn't load reference file into dataframe")
        raise click.Abort()

    return reference_df


def filter_gene(reference_df, input_json):
    '''
    Filter reference dataframe based on genename
    '''
    # filter reference file based on input json
    try:
        LOG.debug("Extracting regions based on genename")
        df = copy.deepcopy(
            reference_df.loc[reference_df['name2'] == input_json['genename'], [
                'chrom', 'exonStarts', 'exonEnds', 'strand', 'name2', 'name',
                'exonCount'
            ]])
        # if entries found in the reference, then exit
        if df.shape[0] == 0:
            LOG.warning(
                "No enteries for input string was found in reference file.")
            raise click.Abort()
    except:
        LOG.error("Extracting genename from reference file failed")
        raise click.Abort()

    return df


def filter_transcript(df, input_json):
    '''
    Filter reference dataframe based on genename
    '''
    # filter reference file based on input json
    try:
        LOG.debug("Extracting regions based on transcript.")
        if 'transcript' in input_json.keys():
            if input_json['transcript']:
                df = df.loc[df['name'] == input_json['transcript']]
                # if entries found in the reference, then exit
                if df.shape[0] == 0:
                    LOG.warning(
                        "No enteries for input string was found in reference file."
                    )
                    raise click.Abort()
    except:
        LOG.error("Extracting transcript from reference file failed")
        raise click.Abort()

    return df


def filter_exon(df, exon_range):
    '''
    Filter bed file dataframe based on exon number and strand information
    '''

    df_bed = pybedtools.BedTool.from_dataframe(df).expand(
        c="2,3").sort().to_dataframe()
    df_out = copy.deepcopy(df_bed)
    if exon_range:
        exonNum = list(range(1, len(df_out) + 1))

        if df.strand.unique()[0] == "-":
            exonNum.sort(reverse=True)

        df_out['thickStart'] = exonNum
        df_out = df_out[df_out['thickStart'].isin(exon_range)]
        df_out['thickStart'] = 'exon_num_' + df_out['thickStart'].astype(str)
    else:
        df_out['thickStart'] = 'total_exon_' + df_out['thickStart'].astype(str)

    return df_out
