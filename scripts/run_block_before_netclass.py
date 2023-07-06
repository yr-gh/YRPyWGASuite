#!/usr/bin/env python

import os
import json
from YRPyWGASuite.workflows import wga_workflow_1

if __name__ == '__main__':
    input_dict_file = "/gpfs/home/yangrui/mywd/align_mm10_mink/input/input.json"
    with open(input_dict_file, 'r') as jf:
        input_dict = json.load(jf)
    
    output_dict = wga_workflow_1(input_dict = input_dict,
                                 workflow_dict = {'block_before_netclass': True,
                                                  'return_before_netclass': True,
                                                  'skip_netclass': True,
                                                  'block_after_netclass': False})
    
    json_file = os.path.dirname(input_dict['root_dir'] + '/') + '/' + 'block_before_netclass_output.json'
    with open(json_file, 'w') as jf:
        json.dump(input_dict, jf)
