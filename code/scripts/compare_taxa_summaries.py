#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

from os.path import basename, join
from cogent.util.misc import create_dir
from qiime.parse import parse_taxa_summary_table
from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option)

from taxcompare.compare_taxa_summaries import (add_filename_suffix,
        comparison_modes, compare_taxa_summaries, correlation_types,
        parse_sample_id_map)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Compares taxa summary files"
script_info['script_description'] = """
This script compares two taxa summary files by computing the correlation
coefficient between pairs of samples. This is useful, for example, if you want
to compare the taxonomic composition of mock communities that were assigned
using different taxonomy assigners in order to see if they are correlated or
not. Another example use-case is to compare the taxonomic composition of
several mock community replicate samples to a single expected, or known, sample
community.

This script is also useful for sorting and filling taxa summary files so that
each sample has the same taxa listed in the same order (with missing taxa
reporting an abundance of zero). The sorted and filled taxa summary files can
then be passed to a script, such as plot_taxa_summary.py, to visually compare
the differences using the same taxa coloring scheme.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("Paired sample comparison",
"Compare all samples that have matching sample IDs between the two input taxa "
"summary files using the pearson correlation coefficient.",
"%prog -i ts_rdp_0.60.txt,ts_rdp_0.80.txt -m paired -o taxa_comp"))
script_info['script_usage'].append(("Paired sample comparison with sample ID "
"map", "Compare all samples that have sample ID mappings in the sample ID "
"map.", "%prog -i ts_rdp_0.60.txt,ts_rdp_0.80.txt -m paired -o taxa_comp"))
script_info['script_usage'].append(("Expected sample comparison",
"Compare all samples in the first taxa summary file to a taxa summary file "
"that contains a single expected, or known, sample using the spearman "
"correlation coefficient.",
"%prog -i ts_rdp_0.60.txt,expected.txt -m expected -o taxa_comp -c spearman"))
script_info['output_description'] = """
The script will always output three files to the specified output directory.
Two files will be the sorted and filled versions of the input taxa summary
files, which can then be used in plot_taxa_summary.py to visualize the
differences in taxonomic composition. The third output file is a tab-separated
file containing the correlation coefficients that were computed between each of
the samples. Each line will contain the sample IDs of the samples that were
compared, followed by the correlation coefficient that was computed.
"""

script_info['required_options'] = [
    make_option('-i', '--taxa_summary_fps', type='existing_filepaths',
        help='the two input taxa summary filepaths, comma-separated. These '
        'will usually be the files that are output by summarize_taxa.py. '
        'These taxa summary files do not need to have the same taxa in the '
        'same order, as the script will make them compatible before comparing '
        'them'),
    options_lookup['output_dir'],
    make_option('-m', '--comparison_mode', type='choice',
        choices=comparison_modes, help='the type of comparison to '
        'perform. Valid choices: ' + ' or '.join(comparison_modes) +
        '. "paired" will compare each sample in the taxa summary '
        'files that match based on sample ID, or that match given a sample ID '
        'map (see the --sample_id_map_fp for more information). "expected" '
        'will compare each sample in the first taxa summary file to an '
        'expected sample (contained in the second taxa summary file). If '
        '"expected", the second taxa summary file must contain only a single '
        'sample that all other samples will be compared to. Note that the '
        '"expected" mode is provided purely for convenience, as this type of '
        'comparison can be accomplished by providing a sample ID map that '
        'maps each sample in the first taxa summary file to a single sample '
        'in the second taxa summary file')
]
script_info['optional_options'] = [
    make_option('-c', '--correlation_type', type='choice',
        choices=correlation_types, help='the type of correlation coefficient '
        'to compute. Valid choices: ' + ' or '.join(correlation_types) +
        ' [default: %default]', default='pearson'),
     make_option('-s','--sample_id_map_fp', type='existing_filepath',
        help='Map of sample IDs from the first taxa summary file to compare '
        'against the sample IDs from the second taxa summary file. Each line '
        'should contain two sample IDs separated by a tab. More than one '
        'sample ID from the first taxa summary file may map to the same '
        'sample ID in the second taxa summary file if desired. This option '
        'only applies if the comparison mode is "paired". If not provided, '
        'only sample IDs that exist in both taxa summary files will be '
        'compared [default: %default]', default=None)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if len(opts.taxa_summary_fps) != 2:
        option_parser.error("Exactly two taxa summary files are required. You "
                            "provided %d." % len(opts.taxa_summary_fps))

    # Create the output dir if it doesn't already exist.
    try:
        create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create or access output directory "
                            "specified with the -o option.")

    sample_id_map = None
    if opts.sample_id_map_fp:
        sample_id_map = parse_sample_id_map(open(opts.sample_id_map_fp, 'U'))

    results = compare_taxa_summaries(
            parse_taxa_summary_table(open(opts.taxa_summary_fps[0], 'U')),
            parse_taxa_summary_table(open(opts.taxa_summary_fps[1], 'U')),
            opts.comparison_mode, correlation_type=opts.correlation_type,
            sample_id_map=sample_id_map)

    # Write out the sorted and filled taxa summaries, basing their
    # filenames on the original input filenames. If the filenames are the same,
    # append a number to each filename.
    same_filenames = False
    if basename(opts.taxa_summary_fps[0]) == \
       basename(opts.taxa_summary_fps[1]):
        same_filenames = True

    for orig_ts_fp, filled_ts_lines, file_num in zip(opts.taxa_summary_fps,
                                                     results[:2], range(0, 2)):
        filename_suffix = '_sorted_and_filled'
        if same_filenames:
            filename_suffix += '_%d' % file_num
        filled_ts_fp = add_filename_suffix(orig_ts_fp, filename_suffix)
        filled_ts_f = open(join(opts.output_dir, filled_ts_fp), 'w')
        filled_ts_f.write(filled_ts_lines)
        filled_ts_f.close()

    # Write the correlation vector.
    corr_vec_f = open(join(opts.output_dir,
                           'taxa_summary_comparison.txt'), 'w')
    corr_vec_f.write(results[2])
    corr_vec_f.close()


if __name__ == "__main__":
    main()
