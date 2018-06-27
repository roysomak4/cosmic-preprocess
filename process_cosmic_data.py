from __future__ import division 
from collections import Counter
import sys
import os
import gzip
import sh
import progressbar


def process_cosmic_db(args):
    if len(args) < 2:
        print ('Insufficient arguments.')
        print('Expected: python process_cosmic_data.py <cosmic_coding.vcf.gz> <cosmic_mutant_export.tsv.gz>')
        sys.exit(1)

    # COSMIC coding mut VCF file
    cosmic_coding_vcf = args[0]

    # COSMIC mut export file
    cosmic_export = args[1]

    # cosmic processed file
    cosmic_output_file = "cosmic_processed.txt"

    # process vcf file first
    unique_variants = {}
    process_vcf(cosmic_coding_vcf, unique_variants)

    # process mut export file
    sites = {}
    process_mut_export(cosmic_export, sites)

    # merge unique variants and sites count data and write to output file
    get_sites_per_variant(unique_variants, sites, cosmic_output_file)


def process_vcf(vcf_file, unique_variants):
    print('Calculating file size...')
    row_count = int(sh.bash("-c", "gunzip -c {0} | wc -l".format(vcf_file)))
    print('Parsing coding cosmic VCF file...')
    counter = 0
    with progressbar.ProgressBar(max_value=row_count) as pbar:
        for vcf_entry in read_vcf(vcf_file):
            temp_arr = vcf_entry.split('\t')
            variant = temp_arr[0] + ':' + temp_arr[1] + ':' + temp_arr[3] + ':' + temp_arr[4]
            cosmic_id = temp_arr[2]
            unique_variants.setdefault(variant, []).append(cosmic_id)
            counter += 1
            pbar.update(counter)
    print('Parsing complete.')
    

def process_mut_export(filename, sites):
    # run a shell command to create intermediate file
    sites_file = 'sites.tmp'
    print('Preprocessing cosmic mutant export file...')
    # check operating system
    gz_decomp = ""
    if sys.platform.startswith("linux"):
        gz_decomp = "zcat"
    elif sys.platform.startswith("darwin"):
        gz_decomp = "gunzip -c"
    else:
        pass

    sh.bash("-c", "{app} {filename} | cut -f 7,8,17 >{sites_file}"
                  .format(app=gz_decomp,
                          filename=filename, 
                          sites_file=sites_file))
    print('Calculating file size...')
    row_count = int(sh.bash("-c", "cat {0} | wc -l".format(sites_file)))
    
    # parse the intermediate file
    print ('Parsing intermediate file...')
    counter = 0
    with progressbar.ProgressBar(max_value=row_count) as pbar:
        for sites_item in read_sites_file(sites_file):
            temp_arr = sites_item.strip('\n').split('\t')
            cosmic_id = temp_arr[2]
            site = temp_arr[1]
            tumor_id = temp_arr[0]
            if sites.has_key(cosmic_id):
                if sites[cosmic_id].has_key(site):
                    sites[cosmic_id][site].add(tumor_id)
                else:
                    sites[cosmic_id][site] = {tumor_id}
            else:
                temp_dict = {}
                temp_dict[site] = {tumor_id}
                sites[cosmic_id] = temp_dict
            counter += 1
            pbar.update(counter)
    print('Parsing complete.')
    print('Removing temporary file...')
    sh.rm(sites_file)
    

def get_sites_per_variant(unique_variants, sites, output_file):
    print ('Processing unique variants and associated sites...')
    with open(output_file, 'w') as outf:
        # output buffer before writing to file
        outbuffer = []
        buffer_cap = 10000
        row_count = 0
        # parse the unique variant list
        counter = 0
        var_count = len(unique_variants)
        with progressbar.ProgressBar(max_value=var_count) as pbar:
            for variant, ids in unique_variants.iteritems():
                # hold column element of each row in a list
                var_line_ele = variant.split(':')
                var_line_ele.append(';'.join(ids))
                # parse each of the cosmic ids to extract the site counts
                temp_site_info = {}
                for cosmic_id in ids:
                    for site, tumor_ids in sites[cosmic_id].iteritems():
                        if temp_site_info.has_key(site):
                            temp_site_info[site] = temp_site_info[site] | tumor_ids
                        else:
                            temp_site_info[site] = tumor_ids
                # parse the temp_site_info
                site_list = []
                for site,tumor_ids in temp_site_info.iteritems():
                    site_list.append(site + '=' + str(len(tumor_ids)))
                var_line_ele.append(';'.join(site_list))
                # check if outbuffer is full
                if row_count > buffer_cap:
                    # write to file
                    outf.write('\n'.join(outbuffer) + '\n')
                    # clear the buffer
                    del outbuffer[:]
                    # reset row_count
                    row_count = 0
                # add current line to buffer
                outbuffer.append('\t'.join(var_line_ele))
                row_count += 1
                counter += 1
                pbar.update(counter)
                # print 'Progress: {0:.2f}%      \r'.format((counter/var_count) * 100),
        # write the last set of lines to the file
        outf.write('\n'.join(outbuffer) + '\n')
        # clear the buffer
        del outbuffer[:]
        print('Cosmic files processed. Cosmic mutation and site counts written to {0}.'.format(output_file)) 
        

def read_vcf(vcf_file):
    with gzip.open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                yield line


def read_sites_file(filename):
    with open(filename, 'r') as input_file:
        for line in input_file:
            if not 'Primary' in line:
                yield line


if __name__ == '__main__':
    process_cosmic_db(sys.argv[1:])
