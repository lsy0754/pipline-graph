import pickle
import urllib.error
import urllib.request

import time
import numpy as np
import os
import csv
import json
import math

from .. import droputils, utils
from ..drop import BarrierAppDROP, BranchAppDrop, ContainerDROP
from ..meta import dlg_float_param, dlg_string_param
from ..meta import dlg_bool_param, dlg_int_param
from ..meta import dlg_component, dlg_batch_input
from ..meta import dlg_batch_output, dlg_streaming_input

from dlg.apps.pyfunc import serialize_data, deserialize_data
from dlg.droputils import DROPFile

##
# @brief H3App
# @details A simple APP that implements the standard Hello World in DALiuGE.
# It allows to change 'World' with some other string and it also permits
# to connect the single output port to multiple sinks, which will all receive
# the same message. App does not require any input.
# @par EAGLE_START
# @param category PythonApp
# @param[in] param/path1 Path1/abc/String/readwrite/
#     \~English path for data
# @param[in] param/appclass Application Class/dlg.apps.makemodel1.H3App/String/readonly/
#     \~English Application class
# @param[in] port/Directions3/String/
#     \~English A CSV file 3
#     \~Chinese \n
# @param[out] port/hello Hello/String/
#     \~English The port carrying the message produced by the app.
# @par EAGLE_END
class H3App(BarrierAppDROP):
    """
    An App that writes 'Hello World!' or 'Hello <greet>!' to all of
    its outputs.
    Keywords:
    greet:   string, [World], whom to greet.
    """
    compontent_meta = dlg_component('H3App', 'H3 App.',
                                    [dlg_batch_input('binary/*', [])],
                                    [dlg_batch_output('binary/*', [])],
                                    [dlg_streaming_input('binary/*')])

    path1 = dlg_string_param('path1', 'World')

    def run(self):
        ins = self.inputs
        #os.system('bash /home/lu/graph/getdata.sh '+self.path1)
        data = self.readdata(self.inputs[0])    #data[0] /home/jywang02/multest/graph/corrected_data/scienceData_SB32039_NGC6744.beam00_averaged_cal.ms
        print('get data :',data)
        datapath=data[0]
        id1 = datapath.find("corrected_data")
        workpath = datapath[0:id1]  #beam00
        modelspath = workpath+'models'
        if not os.path.exists(modelspath):
            os.makedirs(modelspath)
            print('creat{}'.format(modelspath))
            
        id2 = datapath.find("beam")
        beamid = datapath[id2:id2+6]
        beamidnum = beamid[-2:]
        print('beamid',beamid)
        codepath = workpath+'vast-fastimager'
        #os.system('time srun --mpi=pmi2 casa  --logfile {}vast-fastimager_all/log/casa_pre_makemodel.py.{}.log --nologger --nogui -c {}vast-fastimager/pre_makemodel.py {}'.format(workpath,beamidnum,workpath,datapath))
        print('pre_makemodel done')
        
        cutout_model_path = workpath+'cutout_model'
        # if no inputs use the parameter else use the input
        if len(ins) == 0:
            self.greeting = 'H3 %s' % self.greet
        else: # the input is expected to be a vector. We'll use the first element
            # SB32039_beam00.image.tt0.fits  /cutout_model/SB32039_beam00.image.tt0.fits
            self.greeting1=  '{}/SB32039_{}.image.tt0.fits,{}/SB32039_{}.image.tt0.fits'.format(modelspath,beamid,cutout_model_path,beamid)
            self.greeting2=  '{}'.format(datapath)
            
            
        self.outputs[0].write(self.greeting1.encode())
        self.outputs[1].write(self.greeting2.encode())
        print('out2 done')
        
        

    def readdata(self, inDrop):
        a = []

        # NOTE: it appears csv.reader() can't use the DROPFile(inDrop) directly,
        #       since DROPFile is not a iterator. Instead, we read the whole
        #       inDrop to a string and pass that to csv.reader()
        with DROPFile(inDrop) as f:
            file_data = f.read()
            if type(file_data) is bytes:
                file_data=file_data.decode('utf-8')
            csvreader = csv.reader(file_data.split(','))
            for row in csvreader:
                # skip rows with incorrect number of values
                a.append(row[0])
                

        return a
