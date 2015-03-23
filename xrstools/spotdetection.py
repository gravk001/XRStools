import numpy as np
import scipy.misc
from scipy import signal

from scipy.ndimage import maximum_filter
import scipy.ndimage.morphology as morph
import scipy.ndimage.measurements as meas
LENMIN=10

def ReadFile(name):
  data=scipy.misc.imread(name)*1.0
  data=np.max(data, axis=-1)
  return data


def CercaAnelli(A):
  A=A.astype("d")
  # Nmarge=5
  # A[:Nmarge ,:]=1
  # A[:,:Nmarge ]=1
  # A[-Nmarge:,:]=1
  # A[:,-Nmarge:]=1


  Ax = np.concatenate( [A[1:], A[:1] ]  )- np.concatenate( [A[-1:], A[:-1] ]  )
  Ax = 2*Ax + np.concatenate( [Ax[:, 1:], Ax[:, :1] ], axis=1  ) + np.concatenate( [Ax[:,-1:], Ax[:,:-1] ] , axis=1  )


  Ay = np.concatenate( [A[:, 1:], A[:, :1] ] , axis=1  ) - np.concatenate( [A[:, -1:], A[:, :-1] ]  , axis=1   ) 
  Ay = 2*Ay + np.concatenate( [Ay[ 1:], Ay[ :1] ], axis=0  ) + np.concatenate( [Ay[-1:], Ay[:-1] ] , axis=0  )

  res=Canny(Ax,Ay)
  # print len(res)
  
  return res


def IsMaximum( i,j,slope):
  totry = [[i-1,j-1],[i,j-1],[i+1,j-1],[i-1,j],[i+1,j],[i-1,j+1],[i,j+1],[i+1,j+1], ]
  for punto in totry:
   if( slope[punto[0], punto[1]  ]> slope[i,j]):
     return 0
  return 1

def Canny(Ax,Ay):
  Angles = np.arctan2(    Ax, -Ay )

  Angles=Angles*180/np.pi

  Angles = np.floor( Angles/45  + 0.5)
   
  slope = Ax*Ax + Ay*Ay
  
  Nmarge=6
  MaximumRatio=100.0
  

  # maxfits = maximum_filter(slope, size=3)
  # indici = np.where(slope==maxfits)
  # ListLocalMaxima=[ [i,j] for i,j in zip(indici[0], indici[1])   ]
  
  Ni=Ax.shape[0]
  Nj=Ax.shape[1]
  ListLocalMaxima=[]
  for i in range( Nmarge, Ni-Nmarge):
   for j in range( Nmarge, Nj-Nmarge):
     if( IsMaximum( i,j,slope) ):
        ListLocalMaxima.append( [i,j] )

  # print len(ListLocalMaxima)
  # raise

  steps={-5:[-1,1], -4:[-1,0],-3:[-1,-1], -2:[0,-1],-1:[1,-1],0:[1,0],1:[1,1],2:[0,1],3:[-1,1],4:[-1,0], 5:[-1,-1]}
  
  EdgeList=[]

  endpoints=[]



  EdgePoints=np.zeros(Ax.shape  )
  EdgePoints[:Nmarge ,:]=-1
  EdgePoints[:,:Nmarge ]=-1
  EdgePoints[-Nmarge:,:]=-1
  EdgePoints[:,-Nmarge:]=-1
  StartPoints = ListLocalMaxima
  
  for starting in StartPoints:
    EdgePointsTmp = np.zeros(Ax.shape  )  
    s0=starting[0]
    s1=starting[1]
    pentevalue= slope[s0, s1 ]
    if(EdgePoints[ s0, s1  ]==0):
      Edge=[]
      Edge.append( [s0,s1] ) 
      prendi=0
      traccia=0
      while(1):        
        if( slope[s0, s1 ]/pentevalue < 1.0/MaximumRatio):
          prendi=0
          traccia=1
          break  
        EdgePointsTmp[ s0, s1  ]=1
        direction = Angles[s0,s1]
        ss =  [ steps[direction-1],  steps[direction],  steps[direction+1]]
        pttrs = [  [s0+ss[0][0],s1+ss[0][1]], [s0+ss[1][0],s1+ss[1][1]], [s0+ss[2][0],s1+ss[2][1]]  ]
        values =  [ [slope[ pttrs[i][0], pttrs[i][1]  ],i]  for i in range(3)  ]
        imax=    max(values)[1]           
        s0=pttrs[imax][0]
        s1=pttrs[imax][1]                
        if(EdgePoints[ s0, s1  ]):
          prendi=0
          traccia=1
          break        
        if( EdgePointsTmp[ s0, s1  ] ):
          prendi=1
          traccia=1
          Edge.reverse()
          newedge = [ ]
          for p in Edge :
            if tuple(p)== (s0,s1):
              break
            newedge.append(p)
          Edge = [[s0,s1]]+newedge
          break
        else:
          Edge.append( [s0,s1] ) 
      if(prendi):
        endpoints.append([s0,s1])
        EdgeList.append(Edge)
      if traccia:
        for p in Edge:
          EdgePoints[p[0],p[1]]=1
  

  # values =  [ [len(edge),edge]  for edge  in  EdgeList]
  # edge =    max(values)[1]           

  values =  [ edge   for edge  in  EdgeList   if   len(edge)>10   ]
  # edge =    max(values)[1]           

  # return [edge]  
  return  values




def divergenza( Ax,Ay, edge):
  if( edge[-1] != edge[0]):
      edge.append(edge[0])
  res=0
  for i in range(len(edge)-1):
     p1=edge[i]
     p2=edge[i+1]
     Ax1=Ax[p1[0],p1[1] ]
     Ax2=Ax[p2[0],p2[1] ]
     Ay1=Ay[p1[0],p1[1] ]
     Ay2=Ay[p2[0],p2[1] ]
     Dx = p2[0]-p1[0]
     Dy = p2[1]-p1[1]
    
     res=res+ ( Dx*(Ay1+Ay2) -Dy*(Ax1+Ax2) )*0.5
  return res 



def threshold_mask(mask, image,  value  ) :
  n=mask.max()
  npix=3
  
  for i in range(1,n+1):
  
    icount=0
    tmp_mask_old = 0
    while(icount<100):
      tmp_mask = np.equal( i,  mask  )
      
      if not np.any(tmp_mask):
        continue

      mask[:]=mask*(1-tmp_mask)
      
      massimo = (image*tmp_mask).max()
      print " MASSIMO " , massimo
      tmp_mask_2 =   morph.grey_dilation(tmp_mask,  footprint=np.ones([npix,npix]),
                                         structure=np.zeros([npix,npix]))

      print " Npunti , value " ,   tmp_mask.sum(), tmp_mask_2.sum() , value

      tmp_mask = np.less( value*massimo,  image  )*tmp_mask_2

      print " Npunti  " ,   tmp_mask.sum()
     
      
      mask[:]+=tmp_mask*i
     
      if ( tmp_mask-tmp_mask_old).sum()==0:
        break
      tmp_mask_old = tmp_mask


      print " ================ "
      print mask.sum()

      icount+=1

  return mask


def grow_mask(input, npix  ) :
  output =  morph.grey_dilation(input,  footprint=np.ones([npix,npix]),
                                structure=np.zeros([npix,npix]))
  return output


def get_spots_mask( A, median_size=None, nofroi=12,  give_borders=False) :


      A = signal.medfilt2d(A, kernel_size=median_size)
  

      mask=np.zeros(A.shape,"i")
      cerchi = CercaAnelli(A)
      for c in cerchi[:]:
        if len(c)<LENMIN: continue
        for y,x in c:
          mask[y,x] = 1

      filled=morph.binary_fill_holes(mask)

      newmask = relabelise(filled,A, nofroi)

      if give_borders:
        return newmask, mask
      else:
        return newmask


def relabelise(filled, A , nofroi):
      labels, nlabs = meas.label(filled, structure=np.ones([3,3]))
      aves = meas.mean(A ,labels=labels, index = range(1,nlabs+1) )
      avlabs = zip(aves, np.arange(1,nlabs+1))  
      avlabs.sort()
      intes = [ a for a,l in avlabs[-nofroi:]]
      amed = np.median(intes)
      labs = [ l for a,l in avlabs[-nofroi:] if a>amed/100]
      newmask = np.zeros(filled.shape,"i")
      i=1
      for l in labs:
        newmask[np.equal(labels,l)]=i
        i+=1
      return newmask



from PyMca5 import EdfFile
# from PyMca5 import EdfFile


if __name__=="__main__":

  thematrix=EdfFile.EdfFile("test.edf").GetData(0)

  # A=A[:256,:256]

  Nrows , Ncols =  thematrix.shape

  # thematrix=thematrix[Nrows-256:, Ncols-256:  ]
  # Nrows , Ncols =  thematrix.shape

  for i  in range(0,Nrows,256):
    for j  in range(0,Ncols,256):

      submatrix = thematrix[i:i+256,j:j+256]

      A = signal.medfilt2d(submatrix, kernel_size=5)
      submatrix[:]=A

      mask,borders  = get_spots_mask( submatrix, median_size=5, give_borders=True) 

      nspots = mask.max()

      if nspots!=12:
        print "WARNING: LESS SPOTS WERE FOUND :  " ,  nspots
      else:
        print "GOOD!:  FOUND :  " ,  nspots, "  SPOTS " 



      for l in range(1,nspots+1) :
        maskzone = np.equal(mask,l)*np.equal(borders,1)
        submatrix[maskzone] = submatrix.max()


  EdfFile.EdfFile("res.edf","w+").WriteImage({},thematrix)  


