# -*- coding: utf-8 -*-
"""
Wrapped workflows.

Workflow 1:
    References:
        http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto;
        https://doi.org/10.1093/gigascience/giz159 (A genome alignment of 120 mammals highlights ultraconserved element variability and placenta-associated enhancers).
    Workflow:
        1. Masking genome: RepeatModeler --> RepeatMasker.
        2. Aligning: LASTZ.
        3. Chaining: axtChain --> RepeatFiller.py --> chainSort --> chainMergeSort --> chainCleaner.
        4. Netting: chainPreNet --> chainNet --> netSyntenic --> netClass --> NetFilterNonNested.perl.
        5. Maf'ing: netToAxt --> axtSort --> axtToMaf.
        6. PastCons.

@author: Rui Yang
"""

import os
from components import (split_fa, calc_chr_size, fa_to_2bit, par_lastz, par_axt_chain, par_repeat_filler,
                        par_chain_sort, chain_merge_sort, chain_cleaner, chain_pre_net, chain_net,
                        net_syntenic, net_class, net_filter_non_nested, net_to_axt, axt_to_maf)

def wga_workflow_1(input_dict, workflow_dict = {'block_before_netclass': True,
                                                'return_before_netclass': False,
                                                'skip_netclass': False,
                                                'block_after_netclass': True}):
    """Parameters:

           input_dict: a Python dictionary object containing all settings for this workflow.
           
           workflow_dict: to specify which steps to be run or not.
               
               block_before_netclass: run all steps before netClass.
               
               return_before_netclass: return a Python dictionary object containing all inputs and outputs before running netClass.
               
               skip_netclass: skip netClass (do not run netClass).
               
               block_after_netclass: run all steps after netClass.
               
       Workflow 1:
        
       References:
           
           http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto;
           
           https://doi.org/10.1093/gigascience/giz159 (A genome alignment of 120 mammals highlights ultraconserved element variability and placenta-associated enhancers).
       
       Workflow:
           
           1. Masking genome: RepeatModeler --> RepeatMasker (not included here).
           
           2. Aligning: LASTZ.
           
           3. Chaining: axtChain --> RepeatFiller.py --> chainSort --> chainMergeSort --> chainCleaner.
           
           4. Netting: chainPreNet --> chainNet --> netSyntenic --> netClass --> NetFilterNonNested.perl.
           
           5. Maf'ing: netToAxt --> axtSort --> axtToMaf.
           
           6. PastCons (not added yet).

       Note: UCSC database should be configured before running netClass (for your own genome, you may need to build your genome database by yourself; for genomes existing in UCSC database, you can use it via internet or deploy it on your local machine)."""
    if workflow_dict['block_before_netclass']:
        # preprocess reference and query genomes
        # mandatory inputs outside workflow: ref_fa, query_fa, root_dir
        input_dict['ref_fas'] = split_fa(fa = input_dict['ref_fa'],
                                         out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'ref_fas')
        input_dict['query_chrom_sizes'] = calc_chr_size(fa = input_dict['query_fa'],
                                                        out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'chrom_sizes')
        input_dict['ref_chrom_sizes'] = calc_chr_size(fa = input_dict['ref_fa'],
                                                      out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'chrom_sizes')
        input_dict['query_2bit'] = fa_to_2bit(fa = input_dict['query_fa'],
                                              out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'two_bits')
        input_dict['ref_2bit'] = fa_to_2bit(fa = input_dict['ref_fa'],
                                            out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'two_bits')
        
        # run LASTZ
        # mandatory inputs outside workflow: log_dir, lastz_params, lastz_type, concurrent_jobs
        input_dict['lastz.axts'] = par_lastz(ref_fas = input_dict['ref_fas'],
                                             query_fa = input_dict['query_fa'],
                                             out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'lastz_axts',
                                             log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'lastz',
                                             lastz_params = input_dict['lastz_params'],
                                             lastz_type = input_dict['lastz_type'],
                                             concurrent_jobs = input_dict['concurrent_jobs'])
        
        # run axtChain
        # mandatory inputs outside workflow: axt_chain_params
        input_dict['axt_chain.chains'] = par_axt_chain(axts = list(input_dict['lastz.axts'].values()),
                                                       ref_2bit = input_dict['ref_2bit'],
                                                       query_2bit = input_dict['query_2bit'],
                                                       axt_chain_params = input_dict['axt_chain_params'],
                                                       out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'axt_chain_chains',
                                                       log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'axt_chain',
                                                       concurrent_jobs = input_dict['concurrent_jobs'])
        
        # run RepeatFiller.py
        input_dict['repeat_filler.chains'] = par_repeat_filler(chains = list(input_dict['axt_chain.chains'].values()),
                                                               ref_2bit = input_dict['ref_2bit'],
                                                               query_2bit = input_dict['query_2bit'],
                                                               out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'repeat_filler_chains',
                                                               log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'repeat_filler',
                                                               concurrent_jobs = input_dict['concurrent_jobs'])
        
        # run chainSort
        input_dict['chain_sort.chains'] = par_chain_sort(chains = list(input_dict['repeat_filler.chains'].values()),
                                                         out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'chain_sort_chains',
                                                         log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'chain_sort',
                                                         concurrent_jobs = input_dict['concurrent_jobs'])
        
        # run chainMergeSort
        # mandatory inputs outside workflow: chain_merge_sort.prefix
        input_dict['chain_merge_sort.chain'] = chain_merge_sort(chains = list(input_dict['chain_sort.chains'].values()),
                                                                prefix = input_dict['chain_merge_sort.prefix'],
                                                                out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'chain_merge_sort_chain',
                                                                log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'chain_merge_sort')
        
        # run chainCleaner
        # mandatory inputs outside workflow: chain_cleaner_params
        input_dict['chain_cleaner.chain'] = chain_cleaner(chain = input_dict['chain_merge_sort.chain'],
                                                          ref_2bit = input_dict['ref_2bit'],
                                                          query_2bit = input_dict['query_2bit'],
                                                          ref_chr_size = input_dict['ref_chrom_sizes'],
                                                          query_chr_size = input_dict['query_chrom_sizes'],
                                                          chain_cleaner_params = input_dict['chain_cleaner_params'],
                                                          out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'chain_cleaner_chain',
                                                          log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'chain_cleaner')
        
        # run chainPreNet
        input_dict['chain_pre_net.chain'] = chain_pre_net(chain = input_dict['chain_cleaner.chain'][0],
                                                          ref_chr_size = input_dict['ref_chrom_sizes'],
                                                          query_chr_size = input_dict['query_chrom_sizes'],
                                                          out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'chain_pre_net_chain',
                                                          log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'chain_pre_net')
        
        # run chainNet
        # mandatory inputs outside workflow: chain_net_params, chain_net.ref_prefix, chain_net.query_prefix
        input_dict['chain_net.nets'] = chain_net(chain = input_dict['chain_pre_net.chain'],
                                                 ref_2bit = input_dict['ref_2bit'],
                                                 query_2bit = input_dict['query_2bit'],
                                                 ref_chr_size = input_dict['ref_chrom_sizes'],
                                                 query_chr_size = input_dict['query_chrom_sizes'],
                                                 chain_net_params = input_dict['chain_net_params'],
                                                 ref_prefix = input_dict['chain_net.ref_prefix'],
                                                 query_prefix = input_dict['chain_net.query_prefix'],
                                                 out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'chain_net_nets',
                                                 log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'chain_net')
        
        # run netSyntenic
        input_dict['net_syntenic.net'] = net_syntenic(net = input_dict['chain_net.nets'][0],
                                                      out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'net_syntenic_net',
                                                      log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'net_syntenic')
    
    if workflow_dict['return_before_netclass']:
        # do not continue and return here
        return input_dict
    
    if not workflow_dict['skip_netclass']:
        # run netClass
        # mandatory inputs outside workflow: net_class.ref_db, net_class.query_db, net_class_params
        input_dict['net_class.net'] = net_class(net = input_dict['net_syntenic.net'],
                                                ref_db = input_dict['net_class.ref_db'],
                                                query_db = input_dict['net_class.query_db'],
                                                net_class_params = input_dict['net_class_params'],
                                                out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'net_class_net',
                                                log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'net_class')
    
    if workflow_dict['block_after_netclass']:
        # run NetFilterNonNested.perl
        # mandatory inputs outside workflow: net_filter_non_nested_params
        input_dict['net_filter_non_nested.net'] = net_filter_non_nested(net = input_dict['net_class.net'],
                                                                        net_filter_non_nested_params = input_dict['net_filter_non_nested_params'],
                                                                        out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'net_filter_non_nested_net',
                                                                        log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'net_filter_non_nested')
        
        # run netToAxt & axtSort
        input_dict['net_to_axt.axt'] = net_to_axt(net = input_dict['net_filter_non_nested.net'],
                                                  chain = input_dict['chain_pre_net.chain'],
                                                  ref_2bit = input_dict['ref_2bit'],
                                                  query_2bit = input_dict['query_2bit'],
                                                  out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'net_to_axt_axt',
                                                  log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'net_to_axt')
        
        # run axtToMaf
        # mandatory inputs outside workflow: axt_to_maf.ref_prefix, axt_to_maf.query_prefix
        input_dict['axt_to_maf.maf'] = axt_to_maf(axt = input_dict['net_to_axt.axt'],
                                                  ref_chr_size = input_dict['ref_chrom_sizes'],
                                                  query_chr_size = input_dict['query_chrom_sizes'],
                                                  ref_prefix = input_dict['axt_to_maf.ref_prefix'],
                                                  query_prefix = input_dict['axt_to_maf.query_prefix'],
                                                  out_dir = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'axt_to_maf_maf',
                                                  log_dir = os.path.dirname(input_dict['log_dir'] + '/') + '/' + 'axt_to_maf')
    
    return input_dict
