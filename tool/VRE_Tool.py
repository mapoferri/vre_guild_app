#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2020-2022 Barcelona Supercomputing Center (BSC), Spain
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import subprocess
import time
from glob import glob

from basic_modules.tool import Tool
from utils import logger
import re
import yaml
from collections import defaultdict


class myTool(Tool):
    """
    This class define <myTool> Tool.
    """
    DEFAULT_KEYS = ['execution', 'project', 'description']
    """config.json default keys"""
    PYTHON_SCRIPT_PATH = "/workflow.py"
    """<myApplication>"""

    def __init__(self, configuration=None):
        """
        Init function.

        :param configuration: A dictionary containing parameters that define how the operation should be carried out,
            which are specific to <myTool> tool.
        :type configuration: dict
        """
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

        for k, v in self.configuration.items():
            if isinstance(v, list):
                self.configuration[k] = ' '.join(v)

        # Init variables
        self.current_dir = os.path.abspath(os.path.dirname(__file__))
        self.parent_dir = os.path.abspath(self.current_dir + "/../")
        self.execution_path = self.configuration.get('execution', '.')
        if not os.path.isabs(self.execution_path):
            self.execution_path = os.path.normpath(os.path.join(self.parent_dir, self.execution_path))

        self.arguments = dict(
            [(key, value) for key, value in self.configuration.items() if key not in self.DEFAULT_KEYS]
        )

    def run(self, input_files, input_metadata, output_files, output_metadata):
        """
        The main function to run the <myTool> tool.

        :param input_files: Dictionary of input files locations.
        :type input_files: dict
        :param input_metadata: Dictionary of input files metadata.
        :type input_metadata: dict
        :param output_files: Dictionary of output files locations expected to be generated.
        :type output_files: dict
        :param output_metadata: List of output files metadata expected to be generated.
        :type output_metadata: list
        :return: Generated output files and their metadata.
        :rtype: dict, dict
        """
        try:
            # Set and validate execution directory. If not exists the directory will be created
            os.makedirs(self.execution_path, exist_ok=True)

            # Set and validate execution parent directory. If not exists the directory will be created
            execution_parent_dir = os.path.dirname(self.execution_path)
            os.makedirs(execution_parent_dir, exist_ok=True)

            # Update working directory to execution path
            os.chdir(self.execution_path)

            # Tool Execution
            self.toolExecution(input_files)

            # Create and validate the output file from tool execution
            output_id = output_metadata[0]['name']
            #print ('Slef arugment', self.arguments)
            output_type = output_metadata[0]['file']['file_type'].lower()
            output_file_path = (self.execution_path + "/" + "output_" + output_id +"."+output_type)
            
            print ('Output specific', output_id, output_type, output_file_path)

            if os.path.isfile(output_file_path):
                output_files[output_id] = [(output_file_path, "file")]

                return output_files, output_metadata

            # TODO: add more output files to save, if it is necessary for you
            #   or create a method to manage more than one output file

            else:
                errstr = "Output file {} not created. See logs.".format(output_file_path)
                logger.fatal(errstr)
                raise Exception(errstr)

        except:
            errstr = "<myTool> tool execution failed. See logs."
            logger.fatal(errstr)
            raise Exception(errstr)

    def process_args(self):
        #conf_file = os.path.join(global_paths1["step0_remote"]["conf_path"],"BSC_" + prefix +".json")
        #check if it exists, and if should be modified or not
        steps = [key.split(',') for key in self.arguments.keys()]
        #print ('-----------keys', self.arguments.keys())
        #print ('-----------stesp', steps)
        steps_dict = defaultdict(list)  #New dictionary that would save the steps as key and the customized prop as a value;
        l_list = []   #New list for keeping track of the arguments 
        for slist in steps:
            #print ('---slist', slist)
            tmp = ' '.join(slist)
            l_list.append(tmp)
            step = [step.split(':') for step in slist]
            #print ('-----listsssss', l_list)
            #print ('------step', step)
            for mini_list in step:
                #print ('mini_list', mini_list)
                if len(mini_list) == 1:
                    #No specific properties is specified for this step, so no modification is needed;
                    #print ('Step does not require any further property, going forward...')
                    #print ('che cazzo vuol dire', mini_list)
                    continue
                    #Check for the other steps that present props; if necessary, split them by // or any other modifier;
                n_step = mini_list[0]
                key = mini_list[1]
                if '//' in n_step:
                    for n in n_step.split('//'):
                        if n in steps_dict.keys():
                            steps_dict[n].append(key)
                        else:
                            steps_dict[n] = [key]
                        #if n in steps_dict.keys():
                            #steps_dict[n].append(key)
                else:
                    #print ('cosa dovrebbe mette',steps_dict[n_step])
                    if n_step in steps_dict.keys():
                        steps_dict[n_step].append(key)
                    else:
                        steps_dict[n_step] = [key]


        return (l_list, steps_dict)

    def replace_conf(self, input_file, step, prop, sub):
    #Function which is taking in input the YAML config file, the precise step and property to substitue for it, with the given argument chosen by the user (taken from the JSON file)
        with open(input_file, 'r+') as yaml:
            yaml_text = yaml.readlines()
            new_line = ''
            original_line= ''
            for line in yaml_text:
                if str(prop) in line:
                    #Find the right step to modify
                    for line in yaml_text:
                        #if str(prop) in line:
                            #original_line = line
                            #original_line = line.strip()
                            #to_sub = original_line.split(': ')[1]
                            #new_line = original_line.replace(to_sub, sub)
                            #break
                        to_find = r"{0}".format(prop)
                        match = [to_sub.group() for to_sub in re.finditer( to_find, line)]
                        if match:
                 #           print (match)
                            original_line = line.strip()
                 #           print (original_line)
                            subbi = original_line.split(': ')[1]
                            #print (subbi, sub)
                            new_line = original_line.replace(subbi, sub)
                 #           print ('original', original_line, 'to substitute', subbi, sub , 'substituted', new_line)
                            #print ('original', original_line, 'to substitute', to_sub, 'substituted', new_line)
        with open(input_file, 'r') as yaml:
            yaml_text = yaml.read()
            yaml_text = yaml_text.replace(original_line, new_line)
        with open(input_file, 'w') as yaml:
            yaml.write(yaml_text)


    def toolExecution(self, input_files):
        """
        The main function to run the <myTool> tool.

        :param input_files: Dictionary of input files locations.
        :type input_files: dict
        """
        rc = None

        try:
            # Get input files
            # Nextflow config file
            input_file_1 = input_files.get('job_config')
            if not os.path.isabs(input_file_1):
                input_file_1 = os.path.normpath(os.path.join(self.parent_dir, input_file_1))
                print (input_file_1)
            # Params config for the profile

            # Get arguments
            argument_1 = self.arguments.get('profile')
            #if argument_1 is None:
            #    errstr = "Profile must be defined."
            #    logger.fatal(errstr)
            #    raise Exception(errstr)

            l_list, steps_dict = self.process_args()
            #print ('----------list_before_change',l_list)
            #print ('----------steps_before_change', steps_dict) 

            #Overwrite the YAML file given the JSON configuration for the customized MD setup
            if not l_list:
            #if not steps_dict:
               #Dictionary is empty, so no properties per step needs to be customized, and no need to modify the YAML file
               logger.info("Default configuration is maintained, starting Protein MD simulation.")
            else:
                for K,V in steps_dict.items():
              #      print ('---------k and V', K, V)
                    if type(V) == list:
                        for v in steps_dict[K]:
                            match1 = [l for l in l_list if K and v in l]
                            tmp1 = ' '.join(match1)
               #             print ('tmp1', tmp1)
                            value1 = self.arguments.get(tmp1)
               #             print ('value1', value1)
               #             print ('----------match1, tmp1, value1', match1, tmp1, value1)
                            self.replace_conf(input_file_1, K, v, value1)
                    else:
                        match = [l for l in l_list if K in l]
                        tmp = ' '.join(match)
                        value = self.arguments.get(tmp)
               #         print ('----------match, tmp, value', match,tmp, value)
                        self.replace_conf(input_file_1, K, V, value)

            # TODO: add more arguments to use, if it is necessary for you

            # <myApplication> execution
           # if os.path.isfile(input_file_1) and os.path.isfile(input_file_2):

            if os.path.isfile(input_file_1):
                cmd = [
                       'python', 
                       self.parent_dir + self.PYTHON_SCRIPT_PATH, 
                       '--config',
                       input_file_1,
                       '--output_netscore_path',
                       "output_netScore.txt"
                ]

                # TODO: change cmd command line to run <myApplication>

#                cmd = [
#                    'nextflow run',
#                    self.parent_dir + self.PYTHON_SCRIPT_PATH, # workflow.nf
#                    '-C',
#                    input_file_1,  # nextflow.config
#                    '-profile',
#                    argument_1,  # profile
#                ]
                    
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                # Sending the stdout to the log file
                for line in iter(process.stderr.readline, b''):
                    print(line.rstrip().decode("utf-8").replace("", " "))

                rc = process.poll()
                while rc is None:
                    rc = process.poll()
                    time.sleep(0.1)

                if rc is not None and rc != 0:
                    logger.progress("Something went wrong inside the <myApplication> execution. See logs", status="WARNING")
                else:
                    logger.progress("<myApplication> execution finished successfully", status="FINISHED")

            else:
                errstr = "Output file {} not created. See logs.".format(input_file_1)
                #logger.progress("{}".format(input_file_1))
                #errstr = "input_file_1 must be defined."
                logger.fatal(errstr)
                raise Exception(errstr)

        except:
            errstr = "<myApplication> execution failed. See logs."
            logger.error(errstr)
            if rc is not None:
                logger.error("RETVAL: {}".format(rc))
            raise Exception(errstr)
