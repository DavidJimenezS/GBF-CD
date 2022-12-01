# -*- coding: utf-8 -*-
"""
Created on December 2022

graph based data fusion for change detection

FORMAT 
       model = generate_nystrom_graph.GBF_CD(t1, t2, n_samples)
       Change_map = model.gbf_cd()

INPUTS 
       n_samples   - Number of sample nodes (pixels)
       t1          - Single chanel image of time 1
       t2          - Single chanel image of time 2

OUTPUTS
       Change_map  - Binary image of the change zone detected
       labeled_map - Labeled matrix, where 1 is correct detected changes,
                     2 is Missed Alarms (MA), and 3 are False Alarms (FA)
__________________________________________________________________________
Copyright (C) 2022 Graph based data fusion
David Alejandro Jimenez Sierra

@author: david
"""

import numpy as np
import cv2
from numpy.matlib import repmat
from scipy.spatial.distance import cdist
from numpy.linalg import pinv, eig
from scipy.linalg import sqrtm
from skimage.util import img_as_ubyte
from matplotlib import colors

class GBF_CD():

     def __init__(self, t1, t2, n_samples):

          self.n_samples = int(n_samples)
          self.dim1, self.dim2 = t1.shape

          datos = [None] * 2

          ## Data normalization
          datos[0] = np.nan_to_num(np.divide(t1, t1.max()))
          datos[1] = np.nan_to_num(np.divide(t2, t2.max()))

          self.data = datos

     def get_nystrom_fuse_graph(self):
     
          n = len(self.samples)
          n_remainig = len(self.data[0].ravel()) - n
          W = [None] * 2

          for i in range(2):

               samples_d = (self.data[i].ravel())[self.samples]
               samples_d = samples_d[:, np.newaxis] 
               datos = (self.data[i].ravel())[self.complement] 
               datos = datos[:, np.newaxis] 

               distancias_AA = cdist(samples_d,samples_d)
               distancias_AB = np.power(cdist(datos,samples_d), 3)


               D1 = np.matmul(distancias_AA, np.ones((n,1))) + \
                    np.matmul(distancias_AB.T, np.ones((n_remainig,1)))
               distancias_AA = np.divide(distancias_AA, repmat((D1),1,n))

               D2 = np.matmul(distancias_AB, np.ones((n,1))) + \
                    np.matmul(np.matmul(distancias_AB, pinv(distancias_AA)), \
                    np.matmul(distancias_AB.T, np.ones((n_remainig,1))))
                         
               distancias_AB = np.divide(distancias_AB, repmat((D2),1,n))


               kernelstd = distancias_AB.mean()

               distancias_AA = np.exp(np.divide(-(np.power(distancias_AA, 2)) , kernelstd ** 2))
               distancias_AB = np.exp(np.divide(-(np.power(distancias_AB, 2)) , kernelstd ** 2))

               W[i] = np.concatenate((distancias_AA, distancias_AB), axis=0)          


          W_fusionado = np.dstack((W[0], W[1]))
          self.W_fusionado =  W_fusionado.min(axis = 2)

     def get_nystrom_eigs (self):
     
          # One shot method of Nystrom
          n = len(self.samples)
          W_AA = self.W_fusionado[0:n, :]
          W_BA = self.W_fusionado[n:, :]
          
          
          # Nystrom aproximation
          
          W_AA_sqrtinv = pinv(sqrtm(W_AA))
          
          S = W_AA + (W_AA_sqrtinv@(W_BA.T@W_BA)@W_AA_sqrtinv)
          
          [D_s, U_s] = eig(S)
          
          aux_val = np.concatenate((W_AA , (W_AA_sqrtinv @ (W_BA.T)).T), axis=0)
          
          Uhat_W = aux_val@(U_s@pinv(sqrtm(np.diag(D_s))))
          
          A = Uhat_W[0:n,:]
          AB = Uhat_W[n:, :]
          Uhat_W[self.samples, :] = A
          Uhat_W[self.complement, :] = AB
          
          return D_s, Uhat_W

     def mutual_information(self, hgram):

          """ Mutual information for joint histogram
          """
          # Convert bins counts to probability values
          pxy = hgram / float(np.sum(hgram))
          px = np.sum(pxy, axis=1) # marginal for x over y
          py = np.sum(pxy, axis=0) # marginal for y over x
          px_py = px[:, None] * py[None, :] # Broadcast to multiply marginals
          # Now we can do the calculation using the pxy, px_py 2D arrays
          nzs = pxy > 0 # Only non-zero pxy values contribute to the sum
          return np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))

     def get_samples(self):

          nSamplesCol = self.dim2 * np.sqrt(self.n_samples/(self.dim1 * self.dim2 ))
          nSamplesRow = np.round((self.dim1 / self.dim2) * nSamplesCol )
          nSamplesCol = np.round(nSamplesCol)

          dCol = ( self.dim2 / ( nSamplesCol + 1 ) )
          dRow = ( self.dim1 / ( nSamplesRow + 1 ) )

          n1 = np.round((self.dim2 - dCol) / dCol) + 1
          n2 = np.round((self.dim1 - dRow) / dRow) + 1
          samplePointsCol = np.linspace(dCol, self.dim2, n1.astype(np.int64)) - (dCol/2)
          samplePointsRow = np.linspace(dRow, self.dim1, n2.astype(np.int64)) - (dRow/2)

          nSamplesCol = np.size(samplePointsCol)
          nSamplesRow = np.size(samplePointsRow)
          samplePointsRow = samplePointsRow[:, np.newaxis]
          samplePointsCol = samplePointsCol[:, np.newaxis]

          SPRow = repmat(samplePointsRow, 1, nSamplesCol)
          SPCol = repmat(samplePointsCol.T, nSamplesRow, 1)

          samplePoints = np.zeros((2, np.size(SPCol)))

          samplePoints[0,:] = np.reshape(SPCol.T, (1, np.size(SPCol)))
          samplePoints[1,:] = np.reshape(SPRow.T, (1, np.size(SPRow)))
                         
          samplePoints[0,:] = np.minimum(np.maximum(1, samplePoints[0, : ]), self.dim2)
          samplePoints[1,:] = np.minimum(np.maximum(1, samplePoints[1, : ]), self.dim1)
          samples = np.ceil( samplePoints ).astype(np.int64)

          self.samples = np.ravel_multi_index(samples, (self.dim2, self.dim1), order = 'F')[::-1]
          self.complement = np.setdiff1d(np.arange(np.size(self.data[0])), self.samples)

     def data_normalize(self, data):

          normalized = np.divide((data - data.min()), (data.max() - data.min()))

          return normalized

     def gbf_cd(self):

          W_i = [None] * 2

          ## Prior bases in the difference image (DI) forward and backward

          prior1  = np.divide(self.data[0] - self.data[1], self.data[0] + self.data[1])
          prior1  = img_as_ubyte(prior1)

          prior2 = np.divide(self.data[1] - self.data[0], self.data[0] + self.data[1])
          prior2  = img_as_ubyte(prior2)

          _, b1  = cv2.threshold(prior1, 0, 255, cv2.THRESH_OTSU)
          _, b2  = cv2.threshold(prior2, 0, 255, cv2.THRESH_OTSU)
          prior = b1.astype(bool) + b2.astype(bool)

          del b1, b2, prior1, prior2

          #%%
          ## Graph fusion and sampling
          self.get_samples()
          self.get_nystrom_fuse_graph()

          #%%
          ## Change map detection

          Valores_propios, Vectores_propios =  self.get_nystrom_eigs()
          
          MI = np.zeros(self.n_samples)

          for i in range(1, self.n_samples):
          
               Iaux = Vectores_propios[:,i]*np.sqrt(Valores_propios[i])    
               Iaux = np.reshape(Iaux, (self.dim1, self.dim2))

               _, b1  = cv2.threshold(img_as_ubyte(self.data_normalize(np.abs(Iaux))) , 0, 255, cv2.THRESH_OTSU)


               b1 = b1.astype(bool)
               
               hist_2d, _, _ = np.histogram2d(
               prior.ravel(),
               b1.ravel(),
               bins=10)

               MI[i] =  self.mutual_information(hist_2d)

               del Iaux, b1, hist_2d


          idx = np.argmax(MI == np.amax(MI))

          Change_Map = Vectores_propios[:,idx]*np.sqrt(Valores_propios[idx])
          Change_Map = np.reshape(Change_Map, (self.dim1, self.dim2))

          _, b1  = cv2.threshold(img_as_ubyte(self.data_normalize(np.abs(Change_Map))) , 0, 255, cv2.THRESH_OTSU)


          self.Change_Map = b1.astype(bool)

          return self.Change_Map

     def get_rgb_label_map(self, y):

          Err = y.astype(int) - self.Change_Map.astype(int)
          Err[Err == -1] = 3 #FA
          Err[Err == 1] = 2 #MA
          aux_1 = y.astype(int) + self.Change_Map.astype(int)
          Err[aux_1 == 2]= 1 #Correct detections

          self.label_map = Err

          self.cmap = colors.ListedColormap(['white', '#5cff1c', 'blue', 'red'])

          return self.label_map, self.cmap