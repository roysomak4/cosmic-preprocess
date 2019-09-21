#from collections import Counter
import sys
import os
import gzip
import sh
import progressbar
import multiprocessing as mp

# instantiate multiprocessing manager function
mp_manager = mp.Manager()
    
# global variables
max_allele_len = 1000
unique_variants = mp_manager.dict()

def process_cosmic_db(args):
    if len(args) < 2:
        print ('Insufficient arguments.')
        print('Expected: python process_cosmic_data.py <cosmic_coding.vcf.gz> <cosmic_mutant_export.tsv.gz><truncate large variant [true or false]>')
        print ('Truncate large variant is optional. Default = False. If true, removes variants with ref or alt fields >1000bp in size')
        sys.exit(1)

    # COSMIC coding mut VCF file
    cosmic_coding_vcf = args[0]

    # COSMIC mut export file
    cosmic_export = args[1]

    # Truncate large variants beyond a defined size (default = 1000bp)
    truncate_large_var = False
    if len(args) > 2:
        truncate_var = args[2].lower()
        if truncate_var == 'true':
            truncate_large_var = True
        elif truncate_var == 'false':
            truncate_large_var = False
        else:
            print ('Invalid truncate variant argument.')
            print('Expected: python process_cosmic_data.py <cosmic_coding.vcf.gz> <cosmic_mutant_export.tsv.gz><truncate large variant [true or false]>')
            print ('Truncate large variant is optional. Default = False. If true, removes variants with ref or alt fields >1000bp in size')
            sys.exit(1)


    # cosmic processed file
    cosmic_output_file = "cosmic_processed.txt.gz"

    process_cosmic_vcf(cosmic_coding_vcf, truncate_large_var)

    print (unique_variants)
    # process mut export file
    # sites = {}
    # process_mut_export(cosmic_export, sites)

    # # merge unique variants and sites count data and write to output file
    # get_sites_per_variant(unique_variants, sites, cosmic_output_file)


def process_cosmic_vcf(filename, truncate_large_var):
    # init parallel processing
    pool = mp.Pool(10)
    jobs = []
    for vcf_rec in read_vcf(filename):
        jobs.append(pool.apply_async(process_vcf_record, (vcf_rec, truncate_large_var)))
    
    for job in jobs:
        job.get()


def process_vcf_record(vcf_record, truncate_large_var):
    temp_arr = vcf_record.split('\t')
    variant = temp_arr[0] + ':' + temp_arr[1] + ':' + temp_arr[3] + ':' + temp_arr[4]
    cosmic_id = temp_arr[2]
    # if truncation is activated, filter variants with 
    # ref or alt fields > max allele length
    if truncate_large_var:
        if len(temp_arr[3]) <= max_allele_len and len(temp_arr[4]) <= max_allele_len:
            unique_variants.setdefault(variant, []).append(cosmic_id)
    else:
        unique_variants.setdefault(variant, []).append(cosmic_id)


def process_mut_export(cosmic_export, sites):
    # run a shell command to create intermediate file
    sites_file = 'sites.tmp'
    print('Calculating file size...')
    row_count = get_file_row_count(sites_file)
    generate_site_file(cosmic_export, sites_file)
    # parse the intermediate file
    print ('Parsing intermediate file...')
    counter = 0
    with progressbar.ProgressBar(max_value=row_count) as pbar:
        for sites_item in read_sites_file(sites_file):
            temp_arr = sites_item.strip('\n').split('\t')
            cosmic_id = temp_arr[2]
            site = temp_arr[1]
            tumor_id = temp_arr[0]
            if cosmic_id in sites:
                if site in sites[cosmic_id]:
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
    with gzip.open(output_file, 'wb') as outf:
        # output buffer before writing to file
        outbuffer = []
        buffer_cap = 10000
        row_count = 0
        # parse the unique variant list
        counter = 0
        var_count = len(unique_variants)
        with progressbar.ProgressBar(max_value=var_count) as pbar:
            for variant, ids in unique_variants.items():
                # hold column element of each row in a list
                var_line_ele = variant.split(':')
                var_line_ele.append(';'.join(ids))
                # parse each of the cosmic ids to extract the site counts
                temp_site_info = {}
                for cosmic_id in ids:
                    if cosmic_id in sites:
                        for site, tumor_ids in sites[cosmic_id].items():
                            if site in temp_site_info:
                                temp_site_info[site] = temp_site_info[site] | tumor_ids
                            else:
                                temp_site_info[site] = tumor_ids
                    else:
                        print(f'{cosmic_id} not found in sites dictionary')
                # parse the temp_site_info
                site_list = []
                for site,tumor_ids in temp_site_info.items():
                    site_list.append(site + '=' + str(len(tumor_ids)))
                var_line_ele.append(';'.join(site_list))
                # check if outbuffer is full
                if row_count > buffer_cap:
                    # write to file
                    outf.write(('\n'.join(outbuffer) + '\n').encode('UTF-8'))
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
        outf.write(('\n'.join(outbuffer) + '\n').encode('UTF-8'))
        # clear the buffer
        del outbuffer[:]
        print(f'Cosmic files processed. Cosmic mutation and site counts written to {output_file}.')
        

def get_gzip_app():
    gz_app = ''
    if sys.platform.startswith("linux"):
        gz_app = "zcat"
    elif sys.platform.startswith("darwin"):
        gz_app = "gunzip -c"
    else:
        pass
    return gz_app


def get_file_row_count(filename):
    gz_decomp = get_gzip_app()
    print('Calculating number of variants...')
    return int(sh.bash("-c", f"{gz_decomp} {filename} | wc -l"))


def generate_site_file(cosmic_export, sites_file):
    print('Preprocessing cosmic mutant export file...')
    # check operating system
    gz_decomp = get_gzip_app()
    sh.bash("-c", f"{gz_decomp} {cosmic_export} | cut -f 7,8,17 >{sites_file}")


def read_vcf(vcf_file):
    with gzip.open(vcf_file, 'rb') as vcf:
        for line in vcf:
            linestr = line.decode('UTF-8')
            if not linestr.startswith('#'):
                yield linestr


def read_sites_file(filename):
    with open(filename, 'r') as input_file:
        for line in input_file:
            if not 'Primary' in line:
                yield line


def data_chunks(filename, size=1024):
    '''
    Considerations
    1. In python3 relative seeking in non-binary mode is not possible.
    Therefore we will be using gzipped files
    2. os.path.getsize() function on a gzipped file returns the byte size
    of the compressed file which does not match with the length of the
    uncompressed gzip file stream. 
    Therefore we will be using seek to get the size of the file and then 
    iterate over
    '''
    #file_end = os.path.getsize(filename)
    with gzip.open(filename, 'r') as f:
        chunk_end = f.tell()
        while True:
            chunk_start = chunk_end
            f.seek(size, 1)
            f.readline()
            chunk_end = f.tell()
            if chunk_end - chunk_start > 0:
                yield chunk_start, chunk_end - chunk_start
            else:
                break
            

if __name__ == '__main__':
    process_cosmic_db(sys.argv[1:])
