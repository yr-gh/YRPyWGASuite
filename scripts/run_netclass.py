#!/usr/bin/env python

from YRPyWGASuite.components import net_class

if __name__ == '__main__':
    input_dict = {'net_syntenic.net': '/home/yangrui/mywd/netclass/input/mm10.syned.net',
                  'net_class.ref_db': 'mm10',
                  'net_class.query_db': 'slbminkdb',
                  'net_class_params': '-noAr',
                  'out_dir': '/home/yangrui/mywd/netclass/res',
                  'log_dir': '/home/yangrui/mywd/netclass/log'}

    net_class(net = input_dict['net_syntenic.net'],
              ref_db = input_dict['net_class.ref_db'],
              query_db = input_dict['net_class.query_db'],
              net_class_params = input_dict['net_class_params'],
              out_dir = input_dict['out_dir'],
              log_dir = input_dict['log_dir'])
