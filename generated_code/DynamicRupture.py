#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2016, SeisSol Group
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
#

from gemmgen import DB, Tools, Arch, Kernel
import numpy as np
import math

def addMatrices(db, matricesDir, order, dynamicRuptureMethod, numberOfElasticQuantities, numberOfQuantities):
  numberOfBasisFunctions = Tools.numberOfBasisFunctions(order)

  if dynamicRuptureMethod == 'quadrature':
    numberOfPoints = (order+1)**2
  elif dynamicRuptureMethod == 'cellaverage':
    numberOfPoints = int(4**math.ceil(math.log(order*(order+1)/2,4)))
  else:
    raise ValueError('Unknown dynamic rupture method.')

  clones = dict()

  # Load matrices
  db.update(Tools.parseMatrixFile('{}/dr_{}_matrices_{}.xml'.format(matricesDir, dynamicRuptureMethod, order), clones))

  # Determine matrices
  # Note: This does only work because the flux does not depend on the mechanisms in the case of viscoelastic attenuation
  godunovMatrixSpp = np.matlib.zeros((numberOfQuantities,numberOfQuantities))
  godunovMatrixSpp[0:numberOfElasticQuantities,0:numberOfElasticQuantities] = 1.
  db.insert(DB.MatrixInfo('godunovMatrix', numberOfQuantities, numberOfQuantities, matrix=godunovMatrixSpp))
  
  fluxSolverSpp = np.matlib.zeros((numberOfQuantities,numberOfQuantities))
  fluxSolverSpp[0:numberOfElasticQuantities,:] = 1.
  db.insert(DB.MatrixInfo('fluxSolver', numberOfQuantities, numberOfQuantities, matrix=fluxSolverSpp))
  
  godunovStateSpp = np.matlib.zeros((numberOfPoints,numberOfQuantities))
  godunovStateSpp[:,0:numberOfElasticQuantities] = 1.
  db.insert(DB.MatrixInfo('godunovState', numberOfPoints, numberOfQuantities, matrix=godunovStateSpp))

  stiffnessOrder = { 'Xi': 0, 'Eta': 1, 'Zeta': 2 }
  globalMatrixIdRules = [
    (r'^pP(\d{1})$', lambda x: (int(x[0])-1)*4),
    (r'^pM(\d{1})(\d{1})$', lambda x: (int(x[0])-1)*4 + int(x[1])),
    (r'^nP(\d{1})$', lambda x: 16 + (int(x[0])-1)*4),
    (r'^nM(\d{1})(\d{1})$', lambda x: 16 + (int(x[0])-1)*4 + int(x[1]))
  ]
  DB.determineGlobalMatrixIds(globalMatrixIdRules, db, 'dr')

def addKernels(db, kernels, dofMatrixName):
  # Kernels
  for i in range(0,4):
    godunovStatePlus = db['nP{}'.format(i+1)] * db[dofMatrixName] * db['godunovMatrix']
    kernels.append(Kernel.Prototype('godunovState[{}]'.format(i*4), godunovStatePlus, beta=0, prefetch=godunovStatePlus))
    
    flux = db['pP{}'.format(i+1)] * db['godunovState'] * db['fluxSolver']
    kernels.append(Kernel.Prototype('nodalFlux[{}]'.format(i*4), flux, prefetch=db['godunovState']))
    
    for h in range(1,4):
      godunovStateMinus = db['nM{}{}'.format(i+1,h)] * db[dofMatrixName] * db['godunovMatrix']
      kernels.append(Kernel.Prototype('godunovState[{}]'.format(i*4+h), godunovStateMinus, beta=1, prefetch=godunovStateMinus))

      flux = db['pM{}{}'.format(i+1,h)] * db['godunovState'] * db['fluxSolver']
      kernels.append(Kernel.Prototype('nodalFlux[{}]'.format(i*4+h), flux, prefetch=db['godunovState']))
