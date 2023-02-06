#!/usr/bin/env python3

import time
import argparse
import os
import sys
import shutil
from pathlib import Path, PurePath


from biobb_guild import guild
from biobb_guild.guild import parsehippie
from biobb_guild.guild import nodes
from biobb_guild.guild import call_guild
from biobb_guild.guild import gda_disgenet
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu


def prep_output(destination, source):
    wdir = PurePath(source).parents[2]
    if not os.path.isdir(wdir):
        os.mkdir(wdir)
    print (wdir, source, destination)
    print (str(wdir)+'/'+destination)
    file = str(source) +"/output_netScore"
    prova = str(wdir)+'/'+destination
    print (file)
    shutil.copy(file, prova)
    if os.path.isfile(prova):
        print ("File copied.")
    else:
        print ("Some error.")


def main(args):
    start_time= time.time()
    conf = settings.ConfReader(args.config_path)
    print (conf)
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()
    
    global_log.info("step1_gda_retrieval: paths and prop {} ".format(global_prop["step1_gda_retrieval"]))
    
    global_log.info("step1_gda_retrieval: Retrieving file from DisGenet database API")
    gda_disgenet.gda_disgenet(**global_paths["step1_gda_retrieval"], properties=global_prop["step1_gda_retrieval"])

    global_log.info("step2_parse_hippie: Parsing Hippie Network")
    parsehippie.parsehippie(**global_paths["step2_parse_hippie"], properties=global_prop["step2_parse_hippie"])

    global_log.info("step3_making_nodes: Making Network Nodes based on DisGenet results")
    nodes.nodes(**global_paths["step3_making_nodes"], properties=global_prop["step3_making_nodes"])

    global_log.info("step4_call_guild: Calling Guild algo over Network")
    call_guild.call_guild(**global_paths["step4_call_guild"], properties=global_prop["step4_call_guild"])

    prep_output(args.output_netscore_path, global_paths["step4_call_guild"]["output_path"])
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % args.config_path)
    if args.system:
        global_log.info('  System: %s' % system)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Based on the official BioBB tutorial")
    parser.add_argument('--config', dest="config_path", required=True)
    parser.add_argument('--system', dest="system", required=False)
    parser.add_argument('--output_netscore_path', dest='output_netscore_path', required=False)
    args = parser.parse_args()
    main(args)
