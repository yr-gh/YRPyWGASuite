# -*- coding: utf-8 -*-
"""
Building blocks.

@author: Rui Yang
"""

import os, glob
from YRPyShellSuite.baseUtils import log, mkdir_p, run_shell_cmd, del_invalid_files
from YRPyShellSuite.coreUtils import par_thread

def mask_genome(fa, prefix, cpus = 1):
    """Masking genome using RepeatModeler + RepeatMasker.
    
    fa: genome FASTA file.
    
    prefix: containing both output path and prefix (e.g., ~/out/hs).
    
    If '~/out' does not exist, it will be created.
    
    cpus: the number of threads."""
    if not os.path.isfile(fa):
        log.error('{fa} does not exist! Program will exit ...'.format(fa = fa))
        exit(1)

    out_dir = os.path.expanduser(os.path.dirname(prefix))
    mkdir_p(out_dir)
    
    log.info('Starting to build database from {fa} ...'.format(fa = fa))
    builddatabase_log  = os.path.join(out_dir, 'run_builddatabase.log')
    builddatabase_cmd = 'BuildDatabase -name {prefix} {fa} &> {builddatabase_log}'.format(
        prefix = prefix,
        fa = fa,
        builddatabase_log = builddatabase_log)
    # log.info('CMD: {cmd}'.format(cmd = builddatabase_cmd))
    run_shell_cmd(builddatabase_cmd)
    log.info('Building database done!')
    
    log.info('Starting to detect repeat sequences de novo ...')
    repeatmodeler_log = os.path.join(out_dir, 'run_repeatmodeler.log')
    repeatmodeler_cmd = 'RepeatModeler -database {prefix} -pa {cpus} -LTRStruct &> {repeatmodeler_log}'.format(
        prefix = prefix,
        cpus = cpus,
        repeatmodeler_log = repeatmodeler_log)
    # log.info('CMD: {cmd}'.format(cmd = repeatmodeler_cmd))
    run_shell_cmd(repeatmodeler_cmd)
    log.info('Detecting repeat sequences de novo done!')
    
    log.info('Starting to mask genome ({fa}) with lower cases ...'.format(fa = fa))
    repeatmasker_log = os.path.join(out_dir, 'run_repeatmasker.log')
    repeatmasker_cmd = 'RepeatMasker {fa} -lib {lib} -xsmall -s -gff -pa {cpus} &> {repeatmasker_log}'.format(
        fa = fa,
        lib = prefix + '-families.fa',
        cpus = cpus,
        repeatmasker_log = repeatmasker_log)
    # log.info('CMD: {cmd}'.format(cmd = repeatmasker_cmd))
    run_shell_cmd(repeatmasker_cmd)
    log.info('Masking genome done!')

def split_fa(fa, out_dir):
    """Splitting FASTA file into chromosomes by using UCSC faSplit."""
    if not os.path.isfile(fa):
        log.error('{fa} does not exist! Program will exit ...'.format(fa = fa))
        exit(1)

    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/') + '/')
    mkdir_p(out_dir)
    
    log.info('Starting to split {fa} into chromosomes ...'.format(fa = fa))
    cmd = 'faSplit byname {fa} {out_dir}'.format(
        fa = fa,
        out_dir = out_dir)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Splitting {fa} into chromosomes done!'.format(fa = fa))
    
    log.info('Starting to retrieve chromosome-level FATSA files from {out_dir} ...'.format(out_dir = out_dir))
    fas = glob.glob(out_dir + '*.fa')
    if len(fas) > 0:
        log.info('Successfully retrieving {num} chromosome-level FASTA files from {out_dir}.'.format(num = len(fas), out_dir = out_dir))
        return fas
    else:
        log.error('No split FASTA files found!')
        exit(1)

def calc_chr_size(fa, out_dir):
    """Calculating each chromosome size by using UCSC faSize."""
    if not os.path.isfile(fa):
        log.error('{fa} does not exist! Program will exit ...'.format(fa = fa))
        exit(1)

    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    chrom_size = out_dir + '/' + os.path.splitext(os.path.basename(fa))[0] + '.chrom.sizes'
    
    mkdir_p(out_dir)
    
    log.info('Starting to calculate each chromosome size from {fa} ...'.format(fa = fa))
    cmd = 'faSize -detailed {fa} > {chrom_size}'.format(
        fa = fa,
        chrom_size = chrom_size)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Calculating each chromosome size from {fa} done!'.format(fa = fa))
    
    return chrom_size

def fa_to_2bit(fa, out_dir):
    """Converting FASTA file to 2bit file by using UCSC faToTwoBit."""
    if not os.path.isfile(fa):
        log.error('{fa} does not exist! Program will exit ...'.format(fa = fa))
        exit(1)
    
    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    twobit = out_dir + '/' + os.path.splitext(os.path.basename(fa))[0] + '.2bit'
    
    mkdir_p(out_dir)
    
    log.info('Starting to convert {fa} to {twobit} ...'.format(fa = fa, twobit = twobit))
    cmd = 'faToTwoBit {fa} {twobit}'.format(
        fa = fa,
        twobit = twobit)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Converting {fa} to {twobit} done!'.format(fa = fa, twobit = twobit))
    
    return twobit

def par_lastz(ref_fas, query_fa, out_dir, log_dir, lastz_params, lastz_type = 'lastz_32', concurrent_jobs = 1):
    """Performing the pairwise whole genome alignment using lastz/lastz_32 parallelly.
    
    ref_fas: an iterable object, containing a collection of the reference FASTA files.
    
    query_fa: the query FASTA file.
    
    out_dir: the output directory where to output the aligned files.
    
    log_dir: the log directory where to output the log files.
    
    lastz_params: a string of the lastz parameters.
    
    lastz_type: can be either 'lastz' or 'lastz_32'.
    
    concurrent_jobs: the number of jobs to run lastz parallelly."""
    def lastz(ref_fa, query_fa, out_dir, log_dir, lastz_params, lastz_type):
        axt = out_dir + '/' + os.path.splitext(os.path.basename(ref_fa))[0] + '.axt'
        out = log_dir + '/' + os.path.splitext(os.path.basename(ref_fa))[0] + '.out'
        err = log_dir + '/' + os.path.splitext(os.path.basename(ref_fa))[0] + '.err'
        
        log.info('Starting to align {query_fa} to {ref_fa} ...'.format(query_fa = query_fa, ref_fa = ref_fa))
        cmd = '{lastz_type} {ref_fa} {query_fa} {lastz_params} --output={axt} > {out} 2> {err}'.format(
            lastz_type = lastz_type,
            ref_fa = ref_fa, query_fa = query_fa,
            lastz_params = lastz_params,
            axt = axt, out = out, err = err)
        # log.info('CMD: {cmd}'.format(cmd = cmd))
        run_shell_cmd(cmd)
        log.info('Aligning {query_fa} to {ref_fa} done!'.format(query_fa = query_fa, ref_fa = ref_fa))
        
        return axt
    
    ref_fas = del_invalid_files(ref_fas)

    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    res = par_thread(ref_fas, lastz, return_type = 'dict', max_workers = concurrent_jobs, check_files = True, query_fa = query_fa, out_dir = out_dir, log_dir = log_dir, lastz_params = lastz_params, lastz_type = lastz_type)
    
    return res

def par_axt_chain(axts, ref_2bit, query_2bit, axt_chain_params, out_dir, log_dir, concurrent_jobs = 1):
    """Running UCSC axtChain parallelly.
    
    axts: an iterable object, containing a collection of AXT files.
    
    ref_2bit: the reference 2bit file.
    
    query_2bit: the query 2bit file.
    
    axt_chain_params: a string of the axtChain parameters.
    
    out_dir: the output directory where to output the chained files.
    
    log_dir: the log directory where to output the log files.
    
    concurrent_jobs: the number of jobs to run axtChain parallelly.
    
    Note: I strongly recommend using UCSC Kent tools contained in https://github.com/hillerlab/GenomeAlignmentTools if available!
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools."""
    def axt_chain(axt, ref_2bit, query_2bit, axt_chain_params, out_dir, log_dir):
        chain = out_dir + '/' + os.path.splitext(os.path.basename(axt))[0] + '.chain'
        out = log_dir + '/' + os.path.splitext(os.path.basename(axt))[0] + '.out'
        err = log_dir + '/' + os.path.splitext(os.path.basename(axt))[0] + '.err'
        
        log.info('Starting to run axtChain on {axt} ...'.format(axt = axt))
        cmd = 'axtChain {axt_chain_params} {axt} {ref_2bit} {query_2bit} {chain} > {out} 2> {err}'.format(
            axt_chain_params = axt_chain_params,
            axt = axt, ref_2bit = ref_2bit,
            query_2bit = query_2bit, chain = chain,
            out = out, err = err)
        # log.info('CMD: {cmd}'.format(cmd = cmd))
        run_shell_cmd(cmd)
        log.info('Running axtChain on {axt} done!'.format(axt = axt))
        
        return chain
    
    axts = del_invalid_files(axts)

    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    res = par_thread(axts, axt_chain, return_type = 'dict', max_workers = concurrent_jobs, check_files = True, ref_2bit = ref_2bit, query_2bit = query_2bit, axt_chain_params = axt_chain_params, out_dir = out_dir, log_dir = log_dir)
    
    return res

def par_repeat_filler(chains, ref_2bit, query_2bit, out_dir, log_dir, concurrent_jobs = 1):
    """Running RepeatFiller.py parallelly.
    
    chains: an iterable object, containing a collection of CHAIN files.
    
    ref_2bit: the reference 2bit file.
    
    query_2bit: the query 2bit file.
    
    out_dir: the output directory where to output the chained files.
    
    log_dir: the log directory where to output the log files.
    
    concurrent_jobs: the number of jobs to run axtChain parallelly.
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools."""
    def repeat_filler(chain, ref_2bit, query_2bit, out_dir, log_dir):
        filled_chain = out_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.rf.chain'
        out = log_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.out'
        err = log_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.err'
        
        log.info('Starting to run RepeatFiller.py on {chain} ...'.format(chain = chain))
        cmd = 'RepeatFiller.py -c {chain} -T2 {ref_2bit} -Q2 {query_2bit} -o {filled_chain} > {out} 2> {err}'.format(
            chain = chain, ref_2bit = ref_2bit,
            query_2bit = query_2bit, filled_chain = filled_chain,
            out = out, err = err)
        # log.info('CMD: {cmd}'.format(cmd = cmd))
        run_shell_cmd(cmd)
        log.info('Running RepeatFiller.py on {chain} done!'.format(chain = chain))
        
        return filled_chain
    
    chains = del_invalid_files(chains)

    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    res = par_thread(chains, repeat_filler, return_type = 'dict', max_workers = concurrent_jobs, check_files = True, ref_2bit = ref_2bit, query_2bit = query_2bit, out_dir = out_dir, log_dir = log_dir)
    
    return res

def par_chain_sort(chains, out_dir, log_dir, concurrent_jobs = 1):
    """Running UCSC chainSort parallelly.
    
    Note: I strongly recommend using UCSC Kent tools contained in https://github.com/hillerlab/GenomeAlignmentTools if available!
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools.
    
    chains: an iterable object, containing a collection of CHAIN files.
    
    out_dir: the output directory.
    
    log_dir: the log directory.
    
    concurrent_jobs: the number of jobs to run chainSort parallelly."""
    def chain_sort(chain, out_dir, log_dir):
        sorted_chain = out_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.srt.chain'
        out = log_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.out'
        err = log_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.err'
        
        log.info('Starting to run chainSort on {chain} ...'.format(chain = chain))
        cmd = 'chainSort {chain} {sorted_chain} > {out} 2> {err}'.format(
            chain = chain, sorted_chain = sorted_chain,
            out = out, err = err)
        # log.info('CMD: {cmd}'.format(cmd = cmd))
        run_shell_cmd(cmd)
        log.info('Running chainSort on {chain} done!'.format(chain = chain))
        
        return sorted_chain
    
    chains = del_invalid_files(chains)

    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    res = par_thread(chains, chain_sort, return_type = 'dict', max_workers = concurrent_jobs, check_files = True, out_dir = out_dir, log_dir = log_dir)
    
    return res

def chain_merge_sort(chains, out_dir, log_dir, prefix = 'merged'):
    """Running UCSC chainMergeSort.
    
    Note: I strongly recommend using UCSC Kent tools contained in https://github.com/hillerlab/GenomeAlignmentTools if available!
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools.
    
    chains: an iterable object, containing a collection of CHAIN files.
    
    out_dir: the output directory.
    
    log_dir: the log directory.
    
    prefix: the prefix of output file."""
    chains = del_invalid_files(chains)
    
    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    merged_chain = out_dir + '/' + prefix + '.chain'
    err = log_dir + '/' + os.path.splitext(os.path.basename(chains[0]))[0] + '.err'
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    log.info('Starting to merge sorted chain files into a larger chain file ...')
    cmd = 'chainMergeSort {chains} > {merged_chain} 2> {err}'.format(
        chains = ' '.join(chains),
        merged_chain = merged_chain,
        err = err)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Merging sorted chain files done!')
    
    return merged_chain

def chain_cleaner(chain, ref_2bit, query_2bit, ref_chr_size, query_chr_size, chain_cleaner_params, out_dir, log_dir):
    """Running UCSC chainCleaner.
    
    Note: I strongly recommend using UCSC Kent tools contained in https://github.com/hillerlab/GenomeAlignmentTools if available!
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools.
    
    chain: the CHAIN file.
    
    ref_2bit: the reference 2bit file.
    
    query_2bit: the query 2bit file.
    
    ref_chr_size: the reference chromosome sizes.
    
    query_chr_size: the query chromosome sizes.
    
    chain_cleaner_params: the parameters to be used.
    
    out_dir: the output directory.
    
    log_dir: the log directory."""
    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    cleaned_chain = out_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.cleaned.chain'
    removed_chain_bed = out_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.removed.bed'
    out = log_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.out'
    err = log_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.err'
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    log.info('Starting to clean chain file ...')
    cmd = 'chainCleaner {chain} {ref_2bit} {query_2bit} {cleaned_chain} {removed_chain_bed} -tSizes={ref_chr_size} -qSizes={query_chr_size} {chain_cleaner_params} > {out} 2> {err}'.format(
        chain = chain, ref_2bit = ref_2bit, query_2bit = query_2bit,
        cleaned_chain = cleaned_chain, removed_chain_bed = removed_chain_bed,
        ref_chr_size = ref_chr_size, query_chr_size = query_chr_size,
        chain_cleaner_params = chain_cleaner_params,
        out = out, err = err)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Cleaning chain file done!')
    
    return [cleaned_chain, removed_chain_bed]

def chain_pre_net(chain, ref_chr_size, query_chr_size, out_dir, log_dir):
    """Running UCSC chainPreNet.
    
    Note: I strongly recommend using UCSC Kent tools contained in https://github.com/hillerlab/GenomeAlignmentTools if available!
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools.
    
    chain: the CHAIN file.
    
    ref_chr_size: the reference chromosome sizes.
    
    query_chr_size: the query chromosome sizes.
    
    out_dir: the output directory.
    
    log_dir: the log directory."""
    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    prenetted_chain = out_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.prenetted.chain'
    out = log_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.out'
    err = log_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.err'
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    log.info('Starting to pre-net chain file ...')
    cmd = 'chainPreNet {chain} {ref_chr_size} {query_chr_size} {prenetted_chain} > {out} 2> {err}'.format(
        chain = chain, ref_chr_size = ref_chr_size,
        query_chr_size = query_chr_size, prenetted_chain = prenetted_chain,
        out = out, err = err)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Pre-netting chain file done!')
    
    return prenetted_chain

def chain_net(chain, ref_2bit, query_2bit, ref_chr_size, query_chr_size, chain_net_params, out_dir, log_dir, ref_prefix = "target", query_prefix = "query"):
    """Running UCSC chainNet, which is a modified version from https://github.com/hillerlab/GenomeAlignmentTools.
    
    Note: I strongly recommend using UCSC Kent tools contained in https://github.com/hillerlab/GenomeAlignmentTools if available!
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools.
    
    chain: the CHAIN file.
    
    ref_2bit: the reference 2bit file.
    
    query_2bit: the query 2bit file.
    
    ref_chr_size: the reference chromosome sizes.
    
    query_chr_size: the query chromosome sizes.
    
    chain_net_params: the parameters to be used.
    
    out_dir: the output directory.
    
    log_dir: the log directory.
    
    ref_prefix: the prefix for the reference.
    
    query_prefix: the prefix for the query."""
    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    # Reference genome as the template
    ref_net = out_dir + '/' + ref_prefix + '.net'
    # Query genome as the template
    query_net = out_dir + '/' + query_prefix + '.net'
    out = log_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.out'
    err = log_dir + '/' + os.path.splitext(os.path.basename(chain))[0] + '.err'
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    log.info('Starting to net chain file ...')
    cmd = 'chainNet {chain_net_params} -tNibDir={ref_2bit} -qNibDir={query_2bit} {chain} {ref_chr_size} {query_chr_size} {ref_net} {query_net} > {out} 2> {err}'.format(
        chain_net_params = chain_net_params, chain = chain,
        ref_2bit = ref_2bit, query_2bit = query_2bit,
        ref_chr_size = ref_chr_size, query_chr_size = query_chr_size,
        ref_net = ref_net, query_net = query_net,
        out = out, err = err)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Netting chain file done!')
    
    return [ref_net, query_net]

def net_syntenic(net, out_dir, log_dir):
    """Running UCSC netSyntenic.
    
    Note: I strongly recommend using UCSC Kent tools contained in https://github.com/hillerlab/GenomeAlignmentTools if available!
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools.
    
    net: the NET file.
    
    out_dir: the output directory.
    
    log_dir: the log directory."""
    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    syned_net = out_dir + '/' + os.path.splitext(os.path.basename(net))[0] + '.syned.net'
    out = log_dir + '/' + os.path.splitext(os.path.basename(net))[0] + '.out'
    err = log_dir + '/' + os.path.splitext(os.path.basename(net))[0] + '.err'
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    log.info('Starting to run netSyntenic ...')
    cmd = 'netSyntenic {net} {syned_net} > {out} 2> {err}'.format(
        net = net, syned_net = syned_net,
        out = out, err = err)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Running netSyntenic done!')
    
    return syned_net

def net_class(net, ref_db, query_db, net_class_params, out_dir, log_dir):
    """Running UCSC netClass.
    
    Note: I strongly recommend using UCSC Kent tools contained in https://github.com/hillerlab/GenomeAlignmentTools if available!
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools.
    
    net: the NET file.
    
    ref_db: the reference database name.
    
    query_db: the query database name.
    
    net_class_params: the parameters to be used.
    
    out_dir: the output directory.
    
    log_dir: the log directory."""
    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    classed_net = out_dir + '/' + os.path.splitext(os.path.basename(net))[0] + '.classed.net'
    out = log_dir + '/' + os.path.splitext(os.path.basename(net))[0] + '.out'
    err = log_dir + '/' + os.path.splitext(os.path.basename(net))[0] + '.err'
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    log.info('Starting to run netClass ...')
    cmd = 'netClass {net_class_params} {net} {ref_db} {query_db} {classed_net} > {out} 2> {err}'.format(
        net_class_params = net_class_params, net = net,
        ref_db = ref_db, query_db = query_db, classed_net = classed_net,
        out = out, err = err)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Running netClass done!')
    
    return classed_net

def net_filter_non_nested(net, net_filter_non_nested_params, out_dir, log_dir):
    """Running NetFilterNonNested.perl.
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools.
    
    net: the NET file.
    
    net_filter_non_nested_params: the parameters to be used.
    
    out_dir: the output directory.
    
    log_dir: the log directory."""
    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    filtered_net = out_dir + '/' + os.path.splitext(os.path.basename(net))[0] + '.filtered.net'
    err = log_dir + '/' + os.path.splitext(os.path.basename(net))[0] + '.err'
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    log.info('Starting to run NetFilterNonNested.perl ...')
    cmd = 'NetFilterNonNested.perl {net_filter_non_nested_params} {net} > {filtered_net} 2> {err}'.format(
        net_filter_non_nested_params = net_filter_non_nested_params,
        net = net, filtered_net = filtered_net, err = err)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Running NetFilterNonNested.perl done!')
    
    return filtered_net

def net_to_axt(net, chain, ref_2bit, query_2bit, out_dir, log_dir):
    """Running UCSC netToAxt, which converts the NET file to the AXT file.
    
    Note: I strongly recommend using UCSC Kent tools contained in https://github.com/hillerlab/GenomeAlignmentTools if available!
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools.
    
    net: the NET file.
    
    chain: the CHAIN file (I think that the CHAIN file generated by chainPreNet should be fine).
    
    ref_2bit: the reference 2bit file.
    
    query_2bit: the query 2bit file.
    
    out_dir: the output directory.
    
    log_dir: the log directory."""
    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    axt = out_dir + '/' + os.path.splitext(os.path.basename(net))[0] + '.axt'
    out = log_dir + '/' + os.path.splitext(os.path.basename(net))[0] + '.out'
    err = log_dir + '/' + os.path.splitext(os.path.basename(net))[0] + '.err'
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    log.info('Starting to convert NET to AXT ...')
    cmd = 'netToAxt {net} {chain} {ref_2bit} {query_2bit} stdout | axtSort stdin {axt} > {out} 2> {err}'.format(
        net = net, chain = chain,
        ref_2bit = ref_2bit, query_2bit = query_2bit,
        axt = axt, out = out, err = err)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Converting NET to AXT done!')
    
    return axt

def axt_to_maf(axt, ref_chr_size, query_chr_size, out_dir, log_dir, ref_prefix = 'target', query_prefix = 'query'):
    """Running UCSC axtToMaf, which converts the AXT file to the MAF file.
    
    Note: I strongly recommend using UCSC Kent tools contained in https://github.com/hillerlab/GenomeAlignmentTools if available!
    
    For more info, see https://github.com/hillerlab/GenomeAlignmentTools.
    
    axt: the AXT file.
    
    ref_chr_size: the reference chromosome sizes.
    
    query_chr_size: the query chromosome sizes.
    
    out_dir: the output directory.
    
    log_dir: the log directory.
    
    ref_prefix: the prefix for the reference.
    
    query_prefix: the prefix for the query."""
    out_dir = os.path.expanduser(os.path.dirname(out_dir + '/'))
    log_dir = os.path.expanduser(os.path.dirname(log_dir + '/'))
    
    maf = out_dir + '/' + os.path.splitext(os.path.basename(axt))[0] + '.maf'
    out = log_dir + '/' + os.path.splitext(os.path.basename(axt))[0] + '.out'
    err = log_dir + '/' + os.path.splitext(os.path.basename(axt))[0] + '.err'
    
    mkdir_p(out_dir)
    mkdir_p(log_dir)
    
    log.info('Starting to convert AXT to MAF ...')
    cmd = 'axtToMaf -tPrefix={ref_prefix} -qPrefix={query_prefix} {axt} {ref_chr_size} {query_chr_size} {maf} > {out} 2> {err}'.format(
        ref_prefix = ref_prefix, query_prefix = query_prefix,
        axt = axt, ref_chr_size = ref_chr_size,
        query_chr_size = query_chr_size, maf = maf,
        out = out, err = err)
    # log.info('CMD: {cmd}'.format(cmd = cmd))
    run_shell_cmd(cmd)
    log.info('Converting AXT to MAF done!')
    
    return maf
