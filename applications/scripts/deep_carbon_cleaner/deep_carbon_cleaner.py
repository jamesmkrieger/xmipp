#!/usr/bin/env python2
"""/***************************************************************************
 *
 * Authors:    Ruben Sanchez Garcia rsanchez@cnb.csic.es
 *
 * CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
"""

import sys, os
import xmipp_base
from xmipp3 import Plugin
import pyworkflow.em.metadata as md

BAD_IMPORT_MSG='''
Error, tensorflow/keras is probably not installed. Install it with:\n  ./scipion installb deepLearnigToolkit
If gpu version of tensorflow desired, install cuda 8.0 or cuda 9.0
We will try to automatically install cudnn, if unsucesfully, install cudnn and add to LD_LIBRARY_PATH
add to SCIPION_DIR/config/scipion.conf
CUDA = True
CUDA_VERSION = 8.0  or 9.0
CUDA_HOME = /path/to/cuda-%(CUDA_VERSION)
CUDA_BIN = %(CUDA_HOME)s/bin
CUDA_LIB = %(CUDA_HOME)s/lib64
CUDNN_VERSION = 6 or 7
'''

class ScriptCarbonCleanerEm(xmipp_base.XmippScript):
    def __init__(self):

        xmipp_base.XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine('Compute goodness score for picked coordinates. Rule out bad coordinates')
        ## params
        self.addParamsLine(' -i <inputMicsMetada>               : metadata containing the micrograph(s) '+
                           'filenames where coordinates were picked (.mrc or .tif).\n')
                           
        self.addParamsLine(' [ -c <inputCoordsDir>  ] : input coordinates directory (.pos or tab separated x y). Filenames '+
                           'must agree with input micrographs except for file extension.')

        self.addParamsLine(' [ -o <outputCoordsDir> ] : output coordinates directory.')
        self.addParamsLine(' [ -d <deepLearningModel> ]  : (optional) deep learning model filename. If not provided, default model will be used')
        self.addParamsLine('-b <boxSize>     : particles box size in pixels')
        self.addParamsLine('-s <downFactor>   <F=1>   : (optional) micrograph downsampling factor to scale coordinates, Default no scaling')
        self.addParamsLine(' [ --deepThr <deepThr> ]: (optional) deep learning threshold to rule out a coordinate. The smaller the treshold '+
                           'the more coordiantes will be rule out. Ranges 0..1. Recommended 0.75')
        self.addParamsLine(' [--sizeThr <sizeThr> <F=0.8> ]: Failure threshold. Fraction of the micrograph predicted as contamination to ignore predictions. '+
                           '. Ranges 0..1. Default 0.8')
        self.addParamsLine('[ --predictedMaskDir <predictedMaskDir> ] : directory to store the predicted masks. If a given mask already existed, it will be used instead'+
                           ' of a new prediction')
        self.addParamsLine('[ -g <gpuIds>   <N=0> ] : GPU ids to employ. Comma separated list. E.g. "0,1". Default 0. use -1 for CPU-only computation')
        
        ## examples
        self.addExampleLine('xmipp_deep_carbon_cleaner -c path/to/inputCoords/ -o path/to/outputCoords -b $BOX_SIXE  -i  /path/to/micrographs/')
        
    def run(self):

        args={}
        gpusToUse="0"
        if self.checkParam('-g'):
          gpusToUse= self.getParam('-g')
          if "None" in gpusToUse or "-1" in gpusToUse:
            gpusToUse=None
        args["gpus"]=gpusToUse

        updateEnviron(gpusToUse); sys.stdout.flush()
        
        if not self.checkParam('-i'):
          raise Exception("Error, input micrographs fnames are requried as argument")
        else:
          mdObj= md.MetaData(os.path.expanduser( self.getParam('-i')))
          args["inputMicsPath"]= []
          for objId in mdObj:
            args["inputMicsPath"]+= [mdObj.getValue(md.MDL_IMAGE, objId)]

        if not self.checkParam('-b'):
          raise Exception("Error, box size in pixels is required as argument")
        else:
          args["boxSize"]=  self.getIntParam('-b')
          
        args["inputCoordsDir"], args["outputCoordsDir"], args["predictedMaskDir"]= None, None, None
        if self.checkParam('-c'):
          args["inputCoordsDir"]= os.path.expanduser( self.getParam('-c'))
        if self.checkParam('-o'):
          args["outputCoordsDir"]= os.path.expanduser( self.getParam('-o'))
        if self.checkParam('--predictedMaskDir'):
          args["predictedMaskDir"]= os.path.expanduser( self.getParam('--predictedMaskDir'))

        if self.checkParam('-s'):
          args["downFactor"]= self.getDoubleParam('-s')
        else:
          args["downFactor"]= 1

        if self.checkParam('--deepThr'):
          thr= self.getDoubleParam('--deepThr')
          if thr<=0 or thr>=1:
            thr=None
          args["deepThr"]=thr
        else:
          args["deepThr"]= None
                    
        if self.checkParam('--sizeThr'):
          thr=  self.getDoubleParam('--sizeThr')
          if thr<=0 or thr>=1:
            thr=None
          args["sizeThr"]=0.8
        else:
          args["sizeThr"]= 0.8
          
        if args["inputCoordsDir"] is None and args["predictedMaskDir"] is None:
          raise Exception("Either inputCoordsDir or predictedMaskDir (or both) must be provided")
          parser.print_help()
        if args["inputCoordsDir"] is None and args["outputCoordsDir"] is not None:
          raise Exception("Error, if inputCoordsDir provided, then outputCoordsDir must also be provided")
          parser.print_help()
          
        if args["outputCoordsDir"] is None and args["inputCoordsDir"] is not None:
          raise Exception("Error, if outputCoordsDir provided, then inputCoordsDir must also be provided")
    

        if self.checkParam('-d'):
          args["deepLearningModel"]= self.getParam('-d')
        else:
          args["deepLearningModel"]=Plugin.getModel('deepCarbonCleaner', 'defaultModel.keras')
          

        print(args)
        
        try:
          from xmippPyModules.carbon_cleaner_em.cleanMics import main
        except ImportError as e:
          print(e)
          raise ValueError(BAD_IMPORT_MSG)
        main(** args)
 

def updateEnviron(gpus=None):
  """ Create the needed environment for TensorFlow programs. """
  print("updating environ to select gpus: %s"%(gpus) )
  if gpus is not None or gpus is not "":
    os.environ['CUDA_VISIBLE_DEVICES']=str(gpus)
  else:
    os.environ['CUDA_VISIBLE_DEVICES']="-1"
    
if __name__ == '__main__':
    '''
scipion python `which xmipp_deep_carbon_cleaner` -g 0 
    '''
    exitCode=ScriptCarbonCleanerEm().tryRun()
    sys.exit(exitCode)
    