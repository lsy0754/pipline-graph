import pickle
import urllib.error
import urllib.request

import time
import numpy as np
import os
import csv
import json
import math

import numpy as np
import re
from casacore.tables import *
import sys

from .. import droputils, utils
from ..drop import BarrierAppDROP, BranchAppDrop, ContainerDROP
from ..meta import dlg_float_param, dlg_string_param
from ..meta import dlg_bool_param, dlg_int_param
from ..meta import dlg_component, dlg_batch_input
from ..meta import dlg_batch_output, dlg_streaming_input

from dlg.apps.pyfunc import serialize_data, deserialize_data
from dlg.droputils import DROPFile

RAD2DEG = 180.0/math.pi
DEG2RAD = math.pi/180.0


class Skypos:
    """Defines a class that works with spherical geometry, specifically points
    in a unit sphere, such as the sky.

    This is general spherical geometry, with little to tie it to astronomy. The
    exceptions are the naming of longitude and latitude as RA,Dec
    """
    def_precra = 3
    def_precde = 2

    def __init__(self, ra, dec, precra=def_precra, precdec=def_precde):
        """
        Initialise a Skypos object defining a point on a unit sphere with longitude ra and latitude dec
        :param ra: right ascension (radians or hh:mm:ss.ss)
        :type ra: float or str
        :param dec: declination (radians or dd:mm:ss.ss)
        :type dec: float or str
        :param precra:
        :param precdec:
        """
        if isinstance(ra, str):
            self.ra = ras_rad(ra)
            self.dec = decs_rad(dec)
        else:
            self.ra = ra
            self.dec = dec
        self.precra = precra
        self.precdec = precdec
        self.rn = 12+self.precra-Skypos.def_precra
        self.dn = 12+self.precdec-Skypos.def_precde
        self.ras = None
        self.decs = None
        ps = math.pi * 0.5 - self.dec
        sps = math.sin(ps)
        cps = math.cos(ps)
        sra = math.sin(self.ra)
        cra = math.cos(self.ra)
        self._dvecx = [cps * cra, cps * sra, -sps]
        self._dvecy = [-sra, cra, 0.0]
        self._dvecz = [sps * cra, sps * sra, cps]
        self._vec = [cra*sps, sra*sps, cps]

    def rotate_x(self, a):
        """return a skypos determined by rotating self about the X-axis by 
        angle a."""
        x, y, z = _rotate_v_x(self._vec, a)
        b2 = math.asin(z)
        b1 = (2 * math.pi + math.atan2(y, x)) % (2.0 * math.pi)
        return Skypos(b1, b2)

    def rotate_y(self, a):
        """return a skypos determined by rotating self about the X-axis by 
        angle a."""
        x, y, z = _rotatev_y(self._vec, a)
        b2 = math.asin(z)
        b1 = (2 * math.pi + math.atan2(y, x)) % (2.0 * math.pi)
        return Skypos(b1, b2)

    def rotate_z(self, a):
        """return a skypos determined by rotating self about the X-axis by 
        angle a."""
        x, y, z = _rotate_v_z(self._vec, a)
        b2 = math.asin(z)
        b1 = (2 * math.pi + math.atan2(y, x)) % (2.0 * math.pi)
        return Skypos(b1, b2)

    def shift(self, delta_lon, delta_lat):
        """
        Shift this direction (Skypos) in longitude and latitude.
        The longitude shift will be in radian units perpendicular to the direction to pole, along a great circle.
 
        :param float delta_lon: longitude (RA) offset in radians
        :param float delta_lat: latitude (DEC) offset in radians
        """
        lat = self.dec
        lon = self.ra
        # vector along X axis (first point of Aries)
        x0 = Skypos('0h0m0s', '0:0:0', 3, 3)
        shifted_direction = x0.rotate_z(delta_lon).rotate_y(lat + delta_lat).rotate_z(lon)
        return shifted_direction

    def get_ras(self):
        if self.ras is None:
            self.ras = ras(self.ra)
            self.decs = decs(self.dec)
        return self.ras[:self.rn]

    def get_decs(self):
        if self.ras is None:
            self.ras = ras(self.ra)
            self.decs = decs(self.dec)
        return self.decs[:self.dn]

    def __str__(self):
        return '{} {}'.format(self.get_ras(), self.get_decs())

def ras(ra):
    s = ra * (4.0 * 60.0 * RAD2DEG)
    hh = int(s / 3600.0)
    mm = int(s / 60.0) - hh * 60
    ss = s - 60 * (mm + 60 * hh)
    if "{:9.6f}".format(ss) == '60.000000':
        ss = 0.0
        mm += 1
        if mm == 60:
            mm = 0
            hh += 1
            if hh == 24:
                hh = 0
    return "%02d:%02d:%09.6f" % (hh, mm, ss)


def decs(dec):
    s = abs(dec) * (60.0 * 60.0 * RAD2DEG)
    dd = int(s / 3600.0)
    mm = int(s / 60.0) - dd * 60
    ss = s - 60 * (mm + 60 * dd)
    if "%8.5f" % ss == '60.00000':
        ss = 0.0
        mm += 1
        if mm == 60:
            mm = 0
            dd += 1
    sign = ' '
    if dec < 0.0:
        sign = '-'
    return "%s%02d:%02d:%08.6f" % (sign, dd, mm, ss)


def _rotate_v_x(vec, a):
    """Return a skypos determined by rotating vec about the X-axis by 
    angle a."""
    ca, sa = math.cos(a), math.sin(a)
    x = vec[0]
    y = vec[1] * ca - vec[2] * sa
    z = vec[1] * sa + vec[2] * ca
    return [x, y, z]


def _rotatev_y(vec, a):
    """Return a skypos determined by rotating vec about the Y-axis by 
    angle a."""
    ca, sa = math.cos(a), math.sin(a)
    x = vec[0] * ca - vec[2] * sa
    y = vec[1]
    z = vec[0] * sa + vec[2] * ca
    return [x, y, z]


def _rotate_v_z(vec, a):
    """Return a skypos determined by rotating vec about the Z-axis by 
    angle a."""
    ca, sa = math.cos(a), math.sin(a)
    x = vec[0] * ca - vec[1] * sa
    y = vec[0] * sa + vec[1] * ca
    z = vec[2]
    return [x, y, z]


def ras_rad(ra_string):
    """
    Convert right ascension string to radians
    :param ra_string: right ascension string (hh:mm:ss.ss)
    :type ra_string: str
    :return: right ascension in radians
    :rtype: float
    """
    if ra_string[0] == '-':
        raise(ValueError, 'Right ascension may not be negative: {}'.format(ra_string))
    (a, b, c) = re.findall("[0-9.]+", ra_string)
    hh, mm = map(int, [a, b])
    ss = float(c)
    return (ss + 60.0 * (mm + 60.0 * hh)) * 2.0 * math.pi / 86400.0


def decs_rad(dec_string):
    """
    Convert declination string to radians
    :param dec_string: declination string (dd:mm:ss.ss)
    :type dec_string: str
    :return: declination in radians
    :rtype: float
    """
    a, b, c = re.findall('[0-9.]+', dec_string)
    dd, mm = map(int, [a, b])
    ss = float(c)
    r = (ss + 60.0 * (mm + 60.0 * dd)) * 2.0 * math.pi / 1296000.0
    if dec_string[0] == '-':
        r = -r
    return r

def fix(ms):
  
    # Check that the observation wasn't in pol_fixed mode
    ta = table("%s/ANTENNA" %(ms), readonly=True, ack=False)
    ant_mount = ta.getcol("MOUNT", 0, 1)
    if ant_mount[0] != 'equatorial':
        sys.exit("%s doesn't support pol_fixed mode" %(sys.argv[0]))
    ta.close()

    # Work out which beam is in this MS
    t = table(ms, readonly=True, ack=False)
    vis_feed = t.getcol('FEED1', 0, 1)
    beam = vis_feed[0]

    if tableexists("%s/FIELD_OLD" %(ms)) == False:
        print("Making copy of original FIELD table")
        tablecopy(tablename="%s/FIELD" %(ms), newtablename="%s/FIELD_OLD" %(ms))
    else:
        print("Original copy of FIELD table is being used")

    if tableexists("%s/FEED_OLD" %(ms)) == False:
        print("Making copy of original FEED table")
        tablecopy(tablename="%s/FEED" %(ms), newtablename="%s/FEED_OLD" %(ms))
    else:
        print("Original copy of FEED table is being used")

    print("Reading phase directions")
    tp = table("%s/FIELD_OLD" %(ms), readonly=True, ack=False)
    ms_phase = tp.getcol('PHASE_DIR')
    tp.close()

    # Work out how many fields are in the MS.
    n_fields = ms_phase.shape[0]
    print("Found %d fields in FIELD table" %(n_fields))

    # Open up the MS FEED table so we can work out what the offset is for the beam.
    tf = table("%s/FEED" %(ms), readonly=False, ack=False)
    offset = tf.getcol("BEAM_OFFSET")
    offset = offset - offset
    offset = tf.putcol("BEAM_OFFSET", offset)
    tf.close()

    # Open up the MS FIELD table so it can be updated.
    tp = table("%s/FIELD" %(ms), readonly=False, ack=False)
    # Open up the MS FEED table so we can work out what the offset is for the beam.
    tf = table("%s/FEED_OLD" %(ms), readonly=True, ack=False)
    # The offsets are assumed to be the same for all antennas so get a list of all
    # the offsets for one antenna and for the current beam. This should return offsets
    # required for each field.
    t1 = taql("select from $tf where ANTENNA_ID==0 and FEED_ID==$beam")
    n_offsets = t1.getcol("BEAM_OFFSET").shape[0]
    offset_times = t1.getcol("TIME")
    offset_intervals = t1.getcol("INTERVAL")
    print("Found %d offsets in FEED table for beam %d" %(n_offsets, beam))
    for offset_index in range(n_offsets):
        offset = t1.getcol("BEAM_OFFSET")[offset_index]
        print("Offset %d : t=%f-%f : (%fd,%fd)" %(offset_index, offset_times[offset_index]-offset_intervals[offset_index]/2.0,offset_times[offset_index]+offset_intervals[offset_index]/2.0, -offset[0][0]*180.0/np.pi, offset[0][1]*180.0/np.pi))
    
    # Update the beam position for each field
    for field in range(n_fields):
        t = table(ms, readonly=True, ack=False)
        # Get times for the specified field
        tfdata = taql("select from $t where FIELD_ID==$field and FEED1==$beam and ANTENNA1==0 and ANTENNA2==0")
        time_data = tfdata.getcol("TIME")
        # Check if there is data for this field in the measurement set
        if len(time_data) == 0:
            continue        # If not, move on to the next field
        offset_index = -1
        for offset in range(n_offsets):
            if (time_data[0] > offset_times[offset]-offset_intervals[offset]/2.0) and (time_data[0] < offset_times[offset]+offset_intervals[offset]/2.0):
                offset_index = offset
                break

        # print("Field %d : t=%f : offset=%d" %(field, time_data[0], offset_index))
        # Obtain the offset for the current field.
        offset = t1.getcol("BEAM_OFFSET")[offset_index]

        # Get the pointing direction for the field
        p_phase = ms_phase[field]
    
        # Shift the pointing centre by the beam offset
        phase = Skypos(p_phase[0][0], p_phase[0][1], 9, 9)
        new_pos = phase.shift(-offset[0][0], offset[0][1])
        new_pos.rn = 15
        new_pos.dn = 15
        new_pos_str = "%s" %(new_pos)
        print("Setting position of beam %d, field %d to %s (t=%f-%f, offset=%d)" %(beam, field, new_pos_str, time_data[0], time_data[-1], offset_index))
        # Update the FIELD table with the beam position
        new_ra = new_pos.ra
        if new_ra > np.pi:
            new_ra -= 2.0 * np.pi
        ms_phase[field][0][0] = new_ra
        ms_phase[field][0][1] = new_pos.dec
    print('fix done')
    # Write the updated beam positions in to the MS.
    tp.putcol("DELAY_DIR", ms_phase)
    tp.putcol("PHASE_DIR", ms_phase)
    tp.putcol("REFERENCE_DIR", ms_phase)
    tp.close()
    tf.close()
    t.close()
    return phase


##
# @brief H1App
# @details A simple APP that implements the standard Hello World in DALiuGE.
# It allows to change 'World' with some other string and it also permits
# to connect the single output port to multiple sinks, which will all receive
# the same message. App does not require any input.
# @par EAGLE_START
# @param category PythonApp
# @param[in] param/path1 Path1/abc/String/readwrite/
#     \~English path for data
# @param[in] param/appclass Application Class/dlg.apps.downfix4.H1App/String/readonly/
#     \~English Application class
# @param[in] port/Directions3/String/
#     \~English A CSV file 3
#     \~Chinese \n
# @param[out] port/hello Hello/String/
#     \~English The port carrying the message produced by the app.
# @par EAGLE_END
class H1App(BarrierAppDROP):
    """
    An App that writes 'Hello World!' or 'Hello <greet>!' to all of
    its outputs.
    Keywords:
    greet:   string, [World], whom to greet.
    """
    compontent_meta = dlg_component('H1App', 'H1 App.',
                                    [dlg_batch_input('binary/*', [])],
                                    [dlg_batch_output('binary/*', [])],
                                    [dlg_streaming_input('binary/*')])

    path1 = dlg_string_param('path1', 'World')

    def run(self):
        ins = self.inputs
        #os.system('bash /home/lu/graph/getdata.sh '+self.path1)
        data = self.readdata(self.inputs[0])
        print('get data :',data)    
        workpath = self.path1
        path_data=data[2]
        if not os.path.exists(path_data):
            os.makedirs(path_data)
        #os.system('cd /home/lu/graph/{} && pwd && wget {}'.format('out',data[0]))
        os.system('cd {} && pwd && wget {} && tar xvf scienceData.NGC6744.SB32039.NGC6744.beam{}_averaged_cal.leakage.ms.tar && mv scienceData.NGC6744.SB32039.NGC6744.beam{}_averaged_cal.leakage.ms beam{}.ms'.format(data[2],data[0],data[1],data[1],data[1]))
        #os.system('cd {} && pwd && tar xvf beam{}.ms.tar -C {}'.format(data[0],data[1]),data[1],data[1],data[2])
        # if no inputs use the parameter else use the input
        if len(ins) == 0:
            self.greeting = 'H3 %s' % self.greet
        elif len(ins) != 1:
            raise Exception(
                'Only one input expected for %r' % self)
        else: # the input is expected to be a vector. We'll use the first element
            #data[2] /home/jywang02/multest/graph/data
            self.greeting=  '{}/beam{}.ms'.format(path_data,data[1])
          

        outs = self.outputs
        print('self.outputs',outs)
        if len(outs) < 1:
            raise Exception(
                'At least one output should have been added to %r' % self)
        for o in outs:
            o.len = len(self.greeting.encode())
            print('self.outputs length :',len(self.greeting.encode()))
            print('o:',o)
            o.write(self.greeting.encode())  # greet across all outputs

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

##
# @brief H2App
# @details A simple APP that implements the standard Hello World in DALiuGE.
# It allows to change 'World' with some other string and it also permits
# to connect the single output port to multiple sinks, which will all receive
# the same message. App does not require any input.
# @par EAGLE_START
# @param category PythonApp
# @param[in] param/path2 Path2/path/String/readwrite/
#     \~English path for data
# @param[in] param/appclass Application Class/dlg.apps.downfix4.H2App/String/readonly/
#     \~English Application class
# @param[in] port/Directions3/String/
#     \~English A CSV file 3
#     \~Chinese \n
# @param[out] port/hello Hello/String/
#     \~English The port carrying the message produced by the app.
# @par EAGLE_END
class H2App(BarrierAppDROP):
    """
    An App that writes 'Hello World!' or 'Hello <greet>!' to all of
    its outputs.
    Keywords:
    greet:   string, [World], whom to greet.
    """
    compontent_meta = dlg_component('H2App', 'H2 App.',
                                    [dlg_batch_input('binary/*', [])],
                                    [dlg_batch_output('binary/*', [])],
                                    [dlg_streaming_input('binary/*')])

    path2 = dlg_string_param('path2', 'World')

    def run(self):
        ins = self.inputs
        #os.system('bash /home/lu/graph/getdata.sh '+self.path1)
        data = self.readdata(self.inputs[0])   #/home/jywang02/multest/graph/data/beam00.ms
        print('***********************************get data :',data)
        #workpath = self.path2
        datapath=data[0]
        id1 = datapath.find("beam")
        beamid = datapath[id1:id1+6]  #beam00
        print('beamid',beamid)
        id2 = datapath.find("data")
        workpath = datapath[:id2]     #/home/jywang02/multest/graph/
        ms_file = '{}'.format(datapath)
        ms_file_out = '{}corrected_data/scienceData_SB32039_NGC6744.{}_averaged_cal.ms'.format(workpath,beamid)
        if not os.path.exists(workpath+'corrected_data'):
            os.makedirs(workpath+'corrected_data')
            print('creat{}corrected_data'.format(workpath))
        self.rescal(ms_file,ms_file_out)
        print('fix start')
        phase = fix(ms_file_out)
        print('do class done',phase)
        # if no inputs use the parameter else use the input
        if len(ins) == 0:
            self.greeting = 'H3 %s' % self.greet
        elif len(ins) != 1:
            raise Exception(
                'Only one input expected for %r' % self)
        else: # the input is expected to be a vector. We'll use the first element
            #self.greeting=  '{}/beam{}.ms'.format(data[2],data[1])
            self.greeting=  '{}'.format(ms_file_out)
          

        outs = self.outputs
        print('self.outputs',outs)
        if len(outs) < 1:
            raise Exception(
                'At least one output should have been added to %r' % self)
        for o in outs:
            o.len = len(self.greeting.encode())
            print('self.outputs length :',len(self.greeting.encode()))
            print('o:',o)
            o.write(self.greeting.encode())  # greet across all outputs

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
    def rescal(self,path1,path2):
        print("args_file    ="+str(type(path1))+" file="+path1)
        print("args_file_out="+str(type(path2))+" file="+path2) 
        
        os.system("cp -R %s %s" %(path1, path2))
        t = table(path2, readonly=False, ack=False)
        nrows = t.nrows()
        for row in range(nrows):
            if(row % 100000 == 0):
                print("%d/%d" %(row, nrows))
            cdata = t.getcol("DATA", startrow=row, nrow=1)
            cdata *= 2.0
            t.putcol("DATA", cdata, startrow=row, nrow=1)
        t.close()
            
        
        
        
        
        

