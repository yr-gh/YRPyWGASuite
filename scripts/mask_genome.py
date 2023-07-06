#!/usr/bin/env python

from YRPyWGASuite.components import mask_genome

if __name__ == '__main__':
    input_dict = {'fa': '/gpfs/home/yangrui/mywd/mask_mink_genome/input/mink.fa',
                  'prefix': '/gpfs/home/yangrui/mywd/mask_mink_genome/res/mink',
                  'cpus': 16}
    
    mask_genome(fa = input_dict['fa'],
                prefix = input_dict['prefix'],
                cpus = input_dict['cpus'])
