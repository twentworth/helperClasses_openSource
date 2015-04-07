from time import time as getTime
import inspect
import threading
import numpy as np
import os
import math, copy, random, numbers
import operator
import subprocess, atexit
import multiprocessing
import bisect
try:
    import matplotlib.pyplot as pl
except:
    print 'MATPLOTLIB not LOADED'
import pickle#252gb, 207s #Can't do cPickle!!!
try:
    import cPickle
except:
    print 'cPickle NOT LOADED'



class pLambda(object): #a version of the lambda fn feature that is pickable
    def __init__(self, s):
        super(pLambda, self).__init__()
        self.fnStr = s
        self.fn = eval('lambda '+s)
    def __call__(self, *args):
        return self.fn(*args)
    def __getstate__(self):
        return self.fnStr
    def __setstate__(self,d):
        self.__init__(d)
    def __repr__(self):
        return self.fnStr
    def __str__(self):
        return self.fnStr

class modifiedObject(object):
    def __getstate__(self):
        try:
            self.__slots__
        except AttributeError:
            return self.__dict__
        #print 'SLOTS: '+str(self.__class__.__name__)#+' : '+repr(v)
        return tuple([getattr(self,k) for k in self.__slots__])#{k:getattr(self,k) for k in self.__slots__}


        #return {k:getattr(self,k) for k in self.__slots__}
    def __setstate__(self,d):
        if type(d)==tuple:
            i=-1
            #self.__emptyInit__()
            for k in self.__slots__:
                i+=1
                setattr(self,k,d[i])
        else:
            self.__dict__ = d

    def __emptyInit__(self):#This is to just get the class up and running so that the depickiler can do the rest
        a=inspect.getargspec(self.__init__)
        p=0
        if a.defaults is not None:
            p=len(a.defaults)
        ags=[None]*(len(a.args)-1-p)
        self.__init__(*ags)
def runCommand(cmdStr, grepStr=''):
    try:
        cmd = cmdStr.split(' ')
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=-1)
        if len(grepStr)>0:
            p2 = subprocess.Popen(['grep',grepString], stdout=subprocess.PIPE, bufsize=-1, stdin=p1.stdout)
            p=p2
        else:
            p=p1
        iterator = iter(p.stdout.readline,'')
        lines=[]
        for l in iterator:
            if l[-1]=='\n':
                lines.append(l[:-1])
            else:
                lines.append(l)
    except:
        raise
    finally:
        try:
            p1.terminate()
        except:
            pass
        try:
            p2.terminate()
        except:
            pass
    return lines

class fileStream():
    def __init__(self, verbose=True):
        self.flName = None
        #self.batchSize=20
        self.shortFlName = None

        self.sharedDict=None

        self.iterator = None
        self.header = None
        self.separator = None
        self.trailingCharRemove=None
        self.verbose = verbose


        self.lineCount=0
        self.nextDisplay=10
        self.firstDisplay=10
        self.startTime=0
        self.lastTime=0
        self.lastLine=0
    def openFile(self,flName,typ='Normal', grepString=''):
        self.flName=flName
        self.shortFlName = flName.split('/')[-1]

        if typ=='Normal':
            cmd = ['cat', flName]
            #self.process = subprocess.Popen(['cat', flName], stdout=subprocess.PIPE, bufsize=-1)
        elif typ=='bam':
            cmd = ['samtools', 'view', flName]
            #self.process = subprocess.Popen(['samtools', 'view', flName], stdout=subprocess.PIPE, bufsize=-1)
        else:
            raise "bad type"

        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=-1)
        if len(grepString)==0:
            self.processNoGrep=p1
            self.process = p1
        else:
            self.processNoGrep = p1
            self.process = subprocess.Popen(['grep',grepString], stdout=subprocess.PIPE, bufsize=-1, stdin=p1.stdout)
        atexit.register(self.process.terminate)
        self.iterator = iter(self.process.stdout.readline,'')
    def setHeader(self,header):
        self.header=header
        self.sharedDict=sharedDictHolder(keys=header, static=True)#want static=true so we don't store references to everything
    def setSeparator(self,sep):
        self.separator=sep
    def getHeader(self,skipChar=0):
        firstLine = self.iterator.next()[:-1]#remove newline
        if skipChar>0:
            firstLine=firstLine[skipChar:]
        self.setHeader(firstLine.split(self.separator))
    def removeTrailingChar(self,car):
        self.trailingCharRemove = car
    def _processFn(self,line):
        if line is None:
            return line

        line=line[:-1]#REMOVE NEW LINE
        if self.separator is not None:
            line=line.split(self.separator)
            if self.trailingCharRemove is not None:
                for i in xrange(len(line)):
                    if line[i][-len(self.trailingCharRemove):]==self.trailingCharRemove:
                        line[i]=line[i][:-len(self.trailingCharRemove)]
            if self.header is not None:
                line=sharedDict.getSharedDictWithValues(self.sharedDict,line[:len(self.header)])
                #line={self.header[i]:line[i] for i in xrange(len(self.header))}

        else:
            if self.trailingCharRemove is not None:
                if line[-len(self.trailingCharRemove):]==self.trailingCharRemove:
                    line = line[:-len(self.trailingCharRemove)]

        return line
    def __iter__(self):
        return self
    def next(self):
        try:
            self.lineDisp()
            return self._processFn(self.iterator.next())
        except StopIteration:
            if self.verbose:
                print '\r',
                print 'Read '+repr(self.lineCount-1)+' lines from file '+self.shortFlName+'                                                                    '
            raise StopIteration
    def lineDisp(self):
        self.lineCount+=1
        if self.lineCount>=self.nextDisplay and self.verbose:
            if self.lineCount==self.firstDisplay:
                print '\r',
                print 'Reading Line '+repr(self.lineCount)+' at ??? lines/min from file '+self.shortFlName+'                           ',
                self.startTime = getTime()
                self.lastTime = self.startTime
                self.lastLine = self.lineCount
                self.nextDisplay = 50
            else:
                curTime = getTime()
                avglinesPerMin = (self.lineCount-1-self.firstDisplay)/(curTime-self.startTime)*60
                linesPerSec = (self.lineCount-1-self.lastLine)/(curTime-self.lastTime)
                self.lastLine = self.lineCount
                self.lastTime = curTime
                self.nextDisplay+=math.floor(linesPerSec)
                print '\r',
                print 'Reading Line '+intToNiceStr(self.lineCount,2)+' at '+intToNiceStr(linesPerSec*60,2)+' lines/min (Avg: '+intToNiceStr(avglinesPerMin,2)+' lines/min ) from file '+self.shortFlName+'                                ',

class dictMap(object):#Mapps one dict to another
    def __init__(self, headerMap):
        self.otherDict=None
        self.headerMap=headerMap
    def setOtherDict(self,otherDict):
        self.otherDict=otherDict
    def __getitem__(self,key):
        return self.otherDict[self.headerMap[key]]
    def __setitem__(self, key, value):
        self.otherDict[self.headerMap[key]] = value

class sharedDict(modifiedObject):
    __slots__=['data','parent']

    def __init__(self,parent, defaultValue=None, **kwargs):
        #super(sharedDict, self).__init__(iterable, **kwargs)
        self.parent=parent
        try:
            self.data=[defaultValue]*len(self.parent.map)
        except:
            self.data=defaultValue#We need this for depickling
    @classmethod
    def getSharedDictWithValues(cls,parent,values):
        c=cls(parent)
        c.data=list(values)
        return c
    def values(self):
        return list(self.data)

    def __len__(self):
        return len(self.data)

    def keys(self):
        return self.parent.map.keys()

    def has_key(self, key):
        if key in self.parent.map:
            return True
        else:
            return False

    def __iter__(self):
        return self.parent.map.__iter__()

    def __eq__(self, y):
        try:
            if self.parent==y.parent and self.data==y.data:
                return True
            else:
                return False
        except:
            return False

    def items(self):
        return self.data

    def __ne__(self, y):
        return not self.__eq__(y)

    def __cmp__(self, y):
        try:
            if type(y)==type({}):
                return self.getDict()==y
            else:
                return self.getDict()==y.getDict()
        except:
            return False


    def getDict(self):
        map=self.parent.map
        data=self.data
        return {k:data[map[k]] for k in map}

    def setdefault(self, key, default=None):
        if key in self.parent.map:
            return self.data[self.parent.map[key]]
        else:
            return default
    def copy(self):
        return self.getSharedDictWithValues(self.parent,self.values())


    def itervalues(self):
        for key in self.__iter__():
            yield  self.get(key)
    def get(self, key, default=None):
        try:
            return self.__getitem__(key)
        except:
            return default
    def __getattribute__(self, name):
        return super(sharedDict, self).__getattribute__(name)
    def __getitem__(self, key):
        return self.data[self.parent.map[key]]
    def __setitem__(self, key, value):
        try:
            self.data[self.parent.map[key]]=value
        except KeyError:
            self.parent.addKey(key)
            self.data[self.parent.map[key]]=value
        except:
            raise
    def keys(self):
        return self.parent.map.keys()
    def __repr__(self):
        return '{'+', '.join([repr(k)+': '+repr(self[k]) for k in self.__iter__()])+'}'
    def __str__(self):
        return 'Shared Dictionary id:'+str(id(self))+' with map in: id:'+str(id(self.parent))
    def iteritems(self):
        for k in self.__iter__():
            yield (k, self[k])
    def iterkeys(self):
        return self.__iter__()
    def __delitem__(self, key):
        raise TypeError('Cannot delete items from a sharedDict.  You can only delete keys from the sharedDictHolder object. (id:'+str(id(self.parent))+')')

    def __delattr__(self, name):
        raise TypeError('Cannot delete attributes from a sharedDict.')

    def pop(self, key, default=None):
        raise TypeError('Cannot delete items from a sharedDict.  You can only delete keys from the sharedDictHolder object. (id:'+str(id(self.parent))+')')

    def __contains__(self, k):
        if k in self.parent.map:
            return True
        else:
            return False
    def __setattr__(self, name, value):
        return super(sharedDict, self).__setattr__(name, value)
    def update(self, other=None, **kwargs):
        raise TypeError('Cannot delete attributes from a sharedDict.')

    ##########################################################################################################

class sharedDictHolder(modifiedObject):
    __slots__=['sharedDicts','map','static']
    def __init__(self,keys=(),static=False):
        self.static=static
        if static:
            self.sharedDicts=None
        else:
            self.sharedDicts=[]
        self.map={}
        i=0
        if len(keys)>0:
            for k in keys:
                self.map[k]=i
                i+=1
    def renameKey(self,oldKey,newKey):
        if newKey in self.map:
            raise KeyError('New key already used')
        self.map[newKey]=self.map[oldKey]
        del self.map[oldKey]
    def getSharedDict(self, defaultValue=None):
        sd = sharedDict(self, defaultValue=defaultValue)
        if not self.static:
            self.sharedDicts.append(sd)
        return sd
    def getSharedDictWithValues(self,dIct):
        map_abc = self.map
        values=[None]*len(map_abc)
        for k in dIct:
            if k in map_abc:
                values[map_abc[k]]=dIct[k]
        sd = sharedDict.getSharedDictWithValues(self,values=values)
        if not self.static:
            self.sharedDicts.append(sd)
        return sd
    def addKey(self, key, default=None):
        try:
            assert not self.static
        except:
            raise TypeError
        if key not in self.map:
            self.map[key]=len(self.map)
            for sd in self.sharedDicts:
                sd.data.append(default)
    def removeKey(self,key):
        try:
            assert not self.static
        except:
            raise TypeError
        map_abc = self.map
        idx = map_abc[key]
        for k in map_abc:
            if map_abc[k]>idx:
                map_abc[k]-=1
        for sd in self.sharedDicts:
            del sd.data[idx]
        del map_abc[key]
    def keys(self):
        return self.map.keys()
    def __len__(self):
        return len(self.sharedDicts)
    def has_key(self,key):
        return self.map.has_key(key)
    def __iter__(self):
        for sd in self.sharedDicts:
            yield sd
    def getDicts(self):
        try:
            assert not self.static
        except:
            raise TypeError
        return tuple([sd.getDict() for sd in self.sharedDicts])
    def __str__(self):
        return 'sharedDictHolder with '+str(len(self.map))+' keys and '+str(len(self.sharedDicts))+' sharedDicts'
    def __repr__(self):
        return self.__str__()
    def iterkeys(self):
        return self.__iter__()
    def __contains__(self, item):
        if item in self.map:
            return True
        else:
            return False

def nullSpace(A, eps=1e-12):#1e-15 causes problems... I thought we were accurate to 16 digits...
    u, s, vh = np.linalg.svd(A)
    padding = max(0,np.shape(A)[1]-np.shape(s)[0])
    null_mask = np.concatenate(((s <= eps), np.ones((padding,),dtype=bool)),axis=0)
    null_space = np.compress(null_mask, vh, axis=0)
    return np.transpose(null_space)
def listChunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def allMotifIter(l, bases=('A','T','G','C')):

    idxV=[0]*l
    idxV[0]=-1
    try:
        while True:
            idxV[0]+=1
            for i in xrange(len(idxV)):
                if idxV[i]==len(bases):
                    idxV[i]=0
                    idxV[i+1]+=1
            yield ''.join(bases[i] for i in idxV)
    except:
        raise StopIteration
def motifCounter(s,k,dic=None):
    if dic is None:
        dic={}
    for i in xrange(len(s)-k):
        ss=s[i:i+k]
        try:
            dic[ss]+=1
        except KeyError:
            dic[ss]=1


def intToNiceStr(i,d):
    if i>=1000000000:
        return str(round(float(i)/1000000000,d))+' Billion'
    elif i>=1000000:
        return str(round(float(i)/1000000,d))+'M'
    elif i>=1000:
        return str(round(float(i)/1000,d))+'K'
    else:
        return str(round(i,d))

def clopper_pearson(k,n,alpha=0.05):
    """
    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    alpha confidence intervals for a binomial distribution of k expected successes on n trials
    Clopper Pearson intervals are a conservative estimate.
    """


    if k==0:
        lo=0
    else:
        lo = betaDist.ppf(alpha/2, k, n-k+1)
    if n-k==0:
        hi=1
    else:
        hi = betaDist.ppf(1 - alpha/2, k+1, n-k)

    return lo, hi
def clopper_pearson_0toP(k,n,alpha=0.05):
    """
    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    alpha confidence intervals for a binomial distribution of k expected successes on n trials
    Clopper Pearson intervals are a conservative estimate.
    """
    if n-k==0:
        hi=1
    else:
        hi = betaDist.ppf(1 - alpha, k+1, n-k)

    return hi
def clopper_pearson_Pto1(k,n,alpha=0.05):
    """
    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    alpha confidence intervals for a binomial distribution of k expected successes on n trials
    Clopper Pearson intervals are a conservative estimate.
    """
    if k==0:
        lo=0
    else:
        lo = betaDist.ppf(alpha, k, n-k+1)

    return lo
def saveData(d,fileNm, useCPickle=True):
    print 'Saving data to file: '+fileNm
    tmpFl = fileNm+'_tempFile_ASDHFAODSIGHASDKFAS'#+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
    if useCPickle:
        p=cPickle
    else:
        p=pickle
    f = open(tmpFl,'w')
    f.truncate()#delete whatever may have been there

    p.dump(d,f,pickle.HIGHEST_PROTOCOL)

    f.close()
    if os.path.isfile(fileNm):
        os.remove(fileNm)
        os.rename(tmpFl, fileNm)#Doing this shuffle prevents the original file from being deleted if the save program crashes part way through
    else:
        os.rename(tmpFl, fileNm)#Doing this shuffle prevents the original file from being deleted if the save program crashes part way through
def loadData(flName, useCPickle=True):
    print 'Loading data from file: '+flName
    f = open(flName,'r')
    if useCPickle:
        out = cPickle.load(f)
    else:
        out = pickle.load(f)
    f.close()
    return out

class SymArray(object):
    def __init__(self,n):
        self.data=[np.zeros((i+1),dtype=float) for i in xrange(n)]
    def __setitem__(self, (i, j), value):
        if j>i:
            self.data[j][i]=value
        else:
            self.data[i][j]=value
    def __getitem__(self, (i,j)):
        if j>i:
            return self.data[j][i]
        else:
            return self.data[i][j]
    def __repr__(self):
        return 'A '+str(len(self.data))+'x'+str(len(self.data))+' symmetric array.'
    def __str__(self):
        s=[]
        for i in xrange(len(self.data)):
            v=[self[i,j] for j in xrange(len(self.data))]
            s.append(str(v))
        return '['+'\n '.join(s)+']'


def barHist(x,yv, showMean=True, showMeanConf=True, detail=10, tics=None, intervalRange='min-max', nProc=4):
    width=.45
    midPoints = [float(i)-.5 for i in xrange(len(x))]
    try:
        pool = multiprocessing.Pool(nProc)
        rvec=[]
        for i in xrange(len(x)):
            rvec.append( pool.apply_async(gaussSmoothFn, (yv[i],), {'detail':detail} ) )
            #rvec.append(gaussSmoothFn(yv[i], detail=detail))

        for i in xrange(len(x)):
            print float(i)/len(x)
            h,yplts = rvec[i].get() #gaussSmoothFn(yv[i], detail=detail)
            hPlotPos = [midPoints[i]+width*h[k] for k in xrange(len(h))]
            hPlotNeg = [midPoints[i]-width*h[k] for k in xrange(len(h))]
            pl.fill_betweenx(yplts,hPlotNeg,hPlotPos)
            if showMean:
                m=np.mean(yv[i])
                pl.plot([midPoints[i]-width-.02, midPoints[i]+width+.02],[m,m], color='red')
                #bootstrap.ci(data=yBins[i], statfunction=numpy.mean, n_samples=5000)
        if tics is not None:
            pl.xticks(midPoints, tics )
    except KeyboardInterrupt:
        pool.close()
        print 'EXITING!!!'
    finally:
        pool.close()
def hist2d(x,y,nbins=10, logHeight=False, doPlShow=True):
    assert len(x)==len(y)>0
    mx=min(x)
    my=min(y)
    xEdges=np.linspace(mx,max(x),nbins+1).tolist()
    yEdges=np.linspace(my,max(y),nbins+1).tolist()
    cnts = matrix(zerosShape=(nbins,nbins))
    dx=xEdges[1]-xEdges[0]
    dy=yEdges[1]-yEdges[0]
    nbinsm1=nbins-1
    for i in xrange(len(x)):
        xidx=min(int((x[i]-mx)/dx),nbinsm1)
        yidx=min(int((y[i]-my)/dy),nbinsm1)
        if x<xEdges[xidx] and xidx<nbins:
            xidx+=1#fix rounding errors!
        if y<yEdges[yidx] and yidx<nbins:
            yidx+=1#fix rounding errors!

        cnts[xidx,yidx]+=1

    if logHeight:
        mxValLog=math.log(max(cnts))
        mnVal=min(x for x in cnts if x>0)
        mnValLog=math.log(mnVal)
        d=mxValLog-mnValLog
        def fn(x_f):
            return (x_f-mnValLog)/d+1e-16
        for i in xrange(nbins**2):
            if cnts[i]>=mnVal:
                cnts[i]=fn(math.log(cnts[i]))

    else:
        cnts/=max(cnts)


    #NOW DO PLOTTING!

    cmap = pl.get_cmap('jet')
    for i in xrange(nbins):
        for j in xrange(nbins):
            if cnts[i,j]>0:
                pl.fill_between(xEdges[i:i+2],[yEdges[j],yEdges[j]],[yEdges[j+1],yEdges[j+1]],color=cmap(cnts[i,j]))
    if doPlShow:
        pl.show()
def hist(v, nbins=10, edges=None, logHeight=False, normalize=False, outputData=False, makePlot=True, doPlShow=True, color=None, alpha=1, bumpy=False, label=None):
    assert not (logHeight and normalize)
    xmin=min(v)
    xmax=max(v)
    if edges is None:
        edges=np.linspace(xmin,xmax,nbins+1)[1:]
    bins=[0]*(len(edges)+1)
    for x in v:
        i = bisect.bisect_left(edges,x)
        bins[i]+=1
    if logHeight:
        bins = [math.log(x) for x in bins]
    if normalize:
        sumb=sum(bins)
        bins = [float(x)/sumb for x in bins]
    if makePlot:
        if bumpy:
            x=[xmin, xmin, edges[0]]
            y=[0, bins[0], bins[0]]
            for i in xrange(1,len(edges)):
                x+=[edges[i-1],edges[i]]
                y+=[bins[i],bins[i]]
            x+=[ xmax, xmax]
            y+=[ bins[-1],0]
        else:
            x=[xmin, edges[0]]
            y=[0, bins[0]]
            for i in xrange(1,len(edges)):
                x.append(edges[i])
                y.append(bins[i])
            x+=[ xmax, xmax]
            y+=[ bins[-1],0]

        if color is not None:
            pl.plot(x,y, color=color, alpha=alpha, label=label)
        else:
            pl.plot(x,y, alpha=alpha, label=label)
        if doPlShow:
            pl.show()
    if outputData:
        return edges, bins
def movingConfidenceRange(x,y,minWidth=.05, minWidthPts=200,topAlpha=2.5,botAlpha=2.5, doPlShow=True):
    topAlpha=(100.-topAlpha)/100.
    botAlpha/=100.
    xy=zip(x,y)
    xy.sort(key=lambda x_:x_[0])
    xp=[]
    yph=[]
    ypl=[]
    curs=0
    lxy=len(xy)

    cure=min(minWidthPts,lxy)-1

    while cure<=lxy and xy[cure-1][1]-xy[0][1]<minWidth:
        cure+=1
    cursum=sum(x_[0] for x_ in xy[curs:cure])
    sortedList = [x_[1] for x_ in xy[curs:cure]]
    sortedList.sort()
    while cure<lxy:

        cure+=1
        cursum+=xy[cure-1][0]
        bisect.insort(sortedList, xy[cure-1][1])
        while cure-curs>minWidthPts and xy[cure-1][0]-xy[curs][0]>minWidth:
            cursum-=xy[curs][0]
            del sortedList[bisect.bisect_left(sortedList,xy[curs][1])]
            curs+=1
        lenList=cure-curs
        #print sortedList[int(math.ceil(lenList*topAlpha))]
        #ypl.append(xy[curs][0])
        #yph.append( xy[cure-1][0])
        yph.append(sortedList[max(0,int(math.ceil(lenList*topAlpha))-1)])
        ypl.append(sortedList[int(lenList*botAlpha)])
        #xp.append(xy[cure-1][0])
        xp.append(cursum/(cure-curs))

    pl.fill_between(xp,ypl,yph,color='g',alpha=.5)
    if doPlShow:
        pl.show()


def grep(itr,*args):
    for i in itr:
        if type(i)==str and all(strng in i for strng in args):
            yield i
def grepv(itr,*args):
    return [x for x in grep(itr,*args)]
def gaussHistNormalized(y, showMean=False, detail=10):
    h,yplts = gaussSmoothFn(y, detail=detail)
    #Normalize
    sumh=float(sum(h))
    h=[x/sumh for x in h]
    pl.plot(yplts,h)
    if showMean:
        m=np.mean(y)
        pl.plot([m,m],[0,max(h)])


def gaussSmoothFn(v, detail=10, nProc=1):
    minyi=min(v)
    maxyi=max(v)
    d=np.std(v)/detail# centers[1]-centers[0]
    centers = np.linspace(minyi,maxyi,8*math.ceil((maxyi-minyi)/d))
    h=[0]*len(centers)


    for i in xrange(len(centers)):
        center_f=centers[i]

        #def kFn(v_f):
        h[i]=sum(max(0, (d - abs(center_f - x_f))/d) for x_f in v)#kFn(v_f)
        #wtSum=0

        #for x in v:
        #    dif=abs(centers[i]-x)
        #    if dif>d:
        #        continue
        #    wt = (d-dif)/d
        #    wtSum+=wt
        #h[i]=wtSum
    maxH=max(h)
    if maxH>0:
        h=[x/maxH for x in h]
    return h, centers
def triangleKernelFn(x, center, d):
    return max(0, (d-abs(center-x))/d)

def getGaussLine(xx,yy,n,d=None,zeroSeparate=True):
    minxx=min(xx)
    maxxx=max(xx)
    if d is None:
        d=((float(maxxx)-float(minxx))/100)**2
    xxx=np.arange(minxx,maxxx+2*float(maxxx-minxx)/n,float(maxxx-minxx)/n)
    yyy=np.zeros(np.shape(xxx))
    yyyVar=np.zeros(np.shape(xxx))
    for i in xrange(len(xxx)):

        s=0.
        wt=0.
        for j in xrange(len(xx)):
            if zeroSeparate and ((i>0 and xx[j]==0) or (i==0 and xx[j]>0)):
                continue
            curWt=math.exp(-(xxx[i]-xx[j])**2/d)
            wt+=curWt
            s+=curWt*yy[j]
        yyy[i]=s/wt
    for i in xrange(len(xxx)):
        s=0.
        wt=0.
        for j in xrange(len(xx)):
            if zeroSeparate and ((i>0 and xx[j]==0) or (i==0 and xx[j]>0)):
                continue
            curWt=math.exp(-(xxx[i]-xx[j])**2/d)
            #if curWt<.05:
            #    continue
            wt+=curWt
            s+=curWt*(yy[j]-yyy[i])**2
        yyyVar[i]=math.sqrt(s/wt)
    return xxx,yyy,yyyVar

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


class Error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)



class matrix(list):
    def __init__(self, data=None, zerosShape=None, zeroType=0.):
        self.horizontal = True
        self.digitsToShow=8

        if data is not None and zerosShape is not None:
            raise Error('You cannot give both data and zerosShape')
        elif data is None and zerosShape is None:
            raise Error('You Must Definecannot give both data and zerosShape')
        if data is not None:
            self.data = data
            if len(data)>0:
                try:
                    firstLen=len(data[0])

                except:
                    self.data = data = [data]
                    firstLen=len(data[0])
                try:
                    x=(i for i,x in enumerate(data) if not len(x)==firstLen).next()
                    raise Error('Bad Shape!  Row '+str(x)+'is not the correct length!')
                except StopIteration:
                    pass
        elif zerosShape is not None:
            if isinstance(zeroType, numbers.Number):
                self.data = [[zeroType]*zerosShape[1] for _ in xrange(zerosShape[0])]
            elif str.lower(zeroType)=='rand':
                self.data = [[random.random() for __ in xrange(zerosShape[1])] for _ in xrange(zerosShape[0])]
            elif str.lower(zeroType)=='randn':
                self.data = np.random.normal(size=tuple(zerosShape))
            else:
                raise Error('Bad zeroType, must be a number or the string \'rand\' or \'randn\'')
    @property
    def shape(self):
        if self.horizontal:
            return len(self.data),len(self.data[0])
        else:
            return len(self.data[0]),len(self.data)
    @property
    def T(self):
        return self.transpose()
    def transpose(self):
        if self.horizontal:
            m=matrix(data=self.data)
            m.horizontal=False
            return m
        else:
            m=matrix(data=self.data)
            m.horizontal=True
            return m


    def __lt__(self, y):
        raise Error('Unsupported Method \'<\'')


    def __ge__(self, y):
        raise Error('Unsupported Method \'>=\'')


    def __imul__(self, y):### multiplication with assignment
        if isinstance(y, numbers.Number):
            selfShape=self.shape
            for i in xrange(selfShape[0]):
                for j in xrange(selfShape[1]):
                    self[i,j]=self[i,j]*y
            return self
        else:
            raise Error('In place multiplication is not implemented for other matrix or list.')


    def __iadd__(self, y): ### addition with assignment
        if type(self)==type(y):

            try:
                assert self.shape == y.shape
            except:
                raise Error('Incompatible Shapes, '+str(self.shape)+' and '+str(y.shape))
            if self.shape==y.shape:
                for i in xrange(len(y)):
                    self[i]+=y[i]
        elif isinstance(y, numbers.Number):
            selfData=self.data
            for i in xrange(len(selfData)):
                for j in xrange(len(selfData[i])):
                    selfData[i][j]+=y

            #for i in xrange(self.shape[0]*self.shape[1]):
            #    print 'setting'
            #    self[i]=self[i]+y
        else:
            raise Error('Bad input type')
        return self

    def __gt__(self, y):
        return self.l2Norm >= y.l2Norm
    def extend(self, t):
        raise Error('Unsupported Method \'extend,\' use extendDown or extendLeft instead')

    def insert(self, i, x):
        raise Error('Unsupported Method \'insert\'')

    def extendRight(self, t):
        try:
            assert self.shape[0]==t.shape[0]
        except:
            raise Error('Incompatible Shapes, '+str(self.shape)+' and '+str(t.shape))
        selfData=self.data
        if self.horizontal:
            if t.horizontal:

                for i in xrange(len(selfData)):
                    selfData[i].extend(copy.deepcopy(t.data[i]))
            else:
                tData=t.data
                for i in xrange(len(selfData)):
                    selfData[i].extend(x[i] for x in tData)
        else:
            tT=t.T
            self.horizontal=True
            self.extendDown(tT)
            self.horizontal=False

    def extendDown(self, t):
        try:
            assert self.shape[1]==t.shape[1]
        except:
            raise Error('Incompatible Shapes, '+str(self.shape)+' and '+str(t.shape))
        if self.horizontal:
            if t.horizontal:
                tData=t.data #non-deep copy
                self.data.extend(copy.deepcopy(tData[i]) for i in xrange(len(tData)))
            else:
                selfData=self.data
                for i in xrange(t.shape[0]):
                    selfData.append([x[i] for x in tData])
        else:
            tT=t.T
            self.horizontal = not self.horizontal
            self.extendRight(tT)
            self.horizontal = not self.horizontal
    def index(self, x, i=None, j=None):
        raise Error('Unsupported Method \'index\'')


    def reverse(self):
        raise Error('Unsupported Method \'reverse\'')

    def __add__(self, y):
        try:
            assert self.shape == y.shape
        except:
            raise Error('Incompatible Shapes, '+str(self.shape)+' and '+str(y.shape))
        selfCopy = matrix(zerosShape=self.shape)
        for i in xrange(self.shape[0]):
            for j in xrange(self.shape[1]):
                selfCopy[i,j]=self[i,j]+y[i,j]
        return selfCopy
    def __neg__(self):
        selfCopy = matrix(zerosShape=self.shape)
        for i in xrange(self.shape[0]):
            for j in xrange(self.shape[1]):
                selfCopy[i,j]=-self[i,j]
        return selfCopy
    def __sub__(self,y):
        try:
            assert self.shape == y.shape
        except:
            raise Error('Incompatible Shapes, '+str(self.shape)+' and '+str(y.shape))
        selfCopy = matrix(zerosShape=self.shape)
        for i in xrange(self.shape[0]):
            for j in xrange(self.shape[1]):
                selfCopy[i,j]=self[i,j]-y[i,j]
        return selfCopy
    def __isub__(self, y):

        if type(self)==type(y):
            try:
                assert self.shape == y.shape
            except:
                raise Error('Incompatible Shapes, '+str(self.shape)+' and '+str(y.shape))
            for i in xrange(self.shape[0]*self.shape[1]):
                self[i]-=y[i]
        elif isinstance(y, numbers.Number):
            selfData=self.data
            for i in xrange(len(selfData)):
                for j in xrange(len(selfData[i])):
                    selfData[i][j]-=y
        return self

    def append(self, x):
        raise Error('Unsupported Method \'append,\' use extendDown or extendLeft instead')




    def __eq__(self, y):
        if not self.shape==y.shape:
            return False
        shp=self.shape
        for i in xrange(shp[0]):
            for j in xrange(shp[1]):
                if not self[i,j] == y[i,j]:
                    return False
        return True

    def __getslice__(self, i, j):
        return self.__getitem__(i,j)

    def __hash__(self):
        return super(matrix, self).__hash__()

    def __getitem__(self, cord):
        if type(cord)==int:

            selfShape=self.shape
            i=cord/selfShape[1]
            j=cord-i*selfShape[1]
            if not (0<=i<=selfShape[0] and 0<=j<=selfShape[1]):
                raise Error('Bad Coordinates!')
        elif len(cord)==2:
            (i,j)=cord
        else:
            raise Error('Bad Coordinates!')

        if not self.horizontal:
            t=i
            i=j
            j=t
        if not type(i)==slice and not type(j)==slice:
            return self.data[i][j]

        selfData=self.data
        if not type(i)==slice:
            i=slice(i,i+1)
        if not type(j)==slice:
            j=slice(j,j+1)
        if i.stop is None:
            i = slice(i.start,len(selfData), i.step)
        if j.stop is None:
            j = slice(j.start,len(selfData[0]), j.step)
        if i.step is None:
            i=slice(i.start,i.stop,1)
        if j.step is None:
            j=slice(j.start,j.stop,1)
        m = matrix(data=[v[j] for v in selfData[i]])
        m.horizontal=self.horizontal
        return m


    def toLists(self):
        if self.horizontal:
            return copy.deepcopy(self.data)
        else:
            selfData=self.data
            selfShape = self.shape
            return [[selfData[i][j] for i in xrange(selfShape[1])] for j in xrange(selfShape[0])]#selfshape indices reversed!!!
        #d=[[0]*self.shape[0] for _ in xrange(self.shape[1])]
        #for i in xrange(self.shape[0]):
        #    for j in xrange(self.shape[1]):
        #        d[i][j] = self[i,j]
        #return d
    def copy(self):
        return matrix(data=self.toLists())
    def iterindices(self):
        for i,j in ((i,j) for i in self.shape[0] for j in self.shape[1]):
            yield i,j
    def __div__(self, y):
        if isinstance(y, numbers.Number):
            try:
                assert y>0
            except:
                raise ZeroDivisionError
            selfDataCopy=copy.deepcopy(self.data)
            for i in xrange(len(selfDataCopy)):
                for j in xrange(len(selfDataCopy[0])):
                    selfDataCopy[i][j]/=y
            m=matrix(data=selfDataCopy)
            m.horizontal=self.horizontal
            return m
        else:
            raise Error('In place division is not implemented for other matrix or list.')
    def __idiv__(self,y):
        if isinstance(y, numbers.Number):
            try:
                assert y>0
            except:
                raise ZeroDivisionError

            selfShape=self.shape
            selfData = self.data
            for i in xrange(selfShape[0]):
                for j in xrange(selfShape[1]):
                    selfData[i][j]/=y
            return self
        else:
            raise Error('In place division is not implemented for other matrix or list.')
    def __mul__(self, t):
        if isinstance(t, numbers.Number):
            selfData=copy.deepcopy(self.toLists())
            for v in selfData:
                for j in xrange(len(v)):
                    v[j]*=t
            return matrix(data=selfData)

        try:
            inSz = self.shape[1]
            assert self.shape[1]==t.shape[0]
        except:
            raise Error('Inner dimensions must be equal: '+str(self.shape)+', '+str(y.shape))
        m = matrix(zerosShape=(self.shape[0], t.shape[1]))
        for i in xrange(self.shape[0]):
            for j in xrange(t.shape[1]):
                m[i,j]=sum(self[i,k]*t[k,j] for k in xrange(inSz))
        return m

    def __rmul__(self, t):
        if isinstance(t, numbers.Number):
            selfData=copy.deepcopy(self.toLists())
            for v in selfData:
                for j in xrange(len(v)):
                    v[j]*=t
            return matrix(data=selfData)
        raise Error('Possibly bad order of operations, or undefined operation.')
    def __radd__(self, other):
        if type(other)==type(self):
            return self+other
        else:
            raise Error('Possibly bad order of operations, or undefined operation.')
    def __format__(self, *args, **kwargs):
        raise Error('Unsupported Method \'format\'')

    def __ne__(self, y):
        return not self.__eq__(y)

    def __setitem__(self, cord, y):
        if type(cord)==int:
            selfShape=self.shape
            i=cord/selfShape[1]
            j=cord-i*selfShape[1]
            if not (0<=i<=selfShape[0] and 0<=j<=selfShape[1]):
                raise Error('Bad Coordinates!')
        elif len(cord)==2:
            (i,j)=cord
        else:
            raise Error('Bad Coordinates!')
        if not self.horizontal:
            t=i
            i=j
            j=t
        selfData = self.data
        i,j=self.__indicesToSliceForHorizontal(i,j)

        ilen = math.ceil(float(i.stop-i.start)/i.step)
        jlen = math.ceil(float(j.stop-j.start)/j.step)

        if isinstance(y, numbers.Number):
            for ii in xrange(i.start,i.stop,i.step):
                for jj in xrange(j.start,j.stop,j.step):
                    selfData[ii][jj]=y
            return
        if type(y)==list:
            yMtx = matrix(data=y)
        else:
            yMtx=y
        if ilen==yMtx.shape[0] and jlen==yMtx.shape[1]:
            thisIterator = ((x,y) for x in xrange(i.start,i.stop,i.step) for y in xrange(j.start,j.stop,j.step))
            thatIterator = yMtx.__iter__()
            for (ii,jj),val in zip(thisIterator, thatIterator):
                selfData[ii][jj]=val
        elif (ilen==1 or i.step==1) and jlen==yMtx.shape[1]:
            ist=i.start
            for ii in reversed(xrange(ist,i.end,i.step)):
                del selfData[ii]
            for ii in reversed(xrange(yMtx.shape[0])):
                selfData.insert(ist, yMtx[ii,:].toLists()[0])
        elif (jlen==1 or j.step==1) and ilen==yMtx.shape[0]:
            jst=j.start
            for ii in xrange(len(selfData)):
                for jj in reversed(xrange(jst,j.end,j.step)):
                    del selfData[ii][jj]
            for ii in xrange(len(selfData)):
                for val in (yMtx[ii][jj] for jj in reversed(xrange(jst,j.end,j.step))):
                    selfData[ii].insert(jst, val)
        else:

            raise Error('Slicing error, possibly unsupported operation.')



    def __setslice__(self, i, j, y):
        self.__setitem__((i,j),y)

    def __delslice__(self, i, j):
        self.__delitem__((i,j),y)

    def __indicesToSliceForHorizontal(self,i,j):
        if not type(i)==slice:
            i=slice(i,i+1,1)
        if not type(j)==slice:
            j=slice(j,j+1,1)
        if i.start is None:
            i = slice (0,i.stop,i.step)
        if j.start is None:
            j = slice (0,j.stop,j.step)
        if i.stop is None:
            i = slice(i.start,len(self.data), i.step)
        if j.stop is None:
            j = slice(j.start,len(self.data[0]), j.step)
        if i.step is None:
            i=slice(i.start,i.stop,1)
        if j.step is None:
            j=slice(j.start,j.stop,1)
        return i,j
    def __delitem__(self, cord):
        if type(cord)==int:
            if self.shape[0]==1:
                i=0
                j=cord
            elif self.shape[1]==1:
                i=cord
                j=0
            else:
                raise Error('Bad Coordinates!')
        elif len(cord)==2:
            (i,j)=cord
        else:
            raise Error('Bad Coordinates!')
        if not self.horizontal:
            t=i
            i=j
            j=t
        i,j=self.__indicesToSliceForHorizontal(i,j)
        selfData = self.data
        ilen = math.ceil(float(i.stop-i.start)/i.step)
        jlen = math.ceil(float(j.stop-j.start)/j.step)
        shp0=len(selfData)
        shp1=len(selfData[0])
        if not ((ilen==shp0 and i.step==1) or (jlen==shp1 and j.step==1)):
            raise Error('Can only delete whole columns or rows.')
        if ilen==shp0 and jlen==shp1:
            if not i.step==j.step==1:
                raise Error('Can only delete whole columns or rows.')
            self.data=[[]]
        elif ilen==shp0:#del whole column
            for ii in xrange(shp0):
                selfDataii = selfData[ii]
                for jj in reversed(xrange(j.start,j.stop,j.step)):
                    del selfDataii[jj]
        elif jlen==shp1:
            for ii in reversed(xrange(i.start,i.stop,i.step)):
                del selfData[ii]


    def __le__(self, y):
        raise Error('Unsupported Method \'<=\'')

    def __delattr__(self, name):
        raise Error('Unsupported Method \'__delattr__\'')


    def pop(self, i=-1):
        raise Error('Unsupported Method \'pop\'')

    def __reversed__(self):
        raise Error('Unsupported Method \'reversed\'')


    def count(self, x):
        cnt=0
        for y in self:
            if y==x:
                cnt+=1
        return cnt
    def reshape(self, newShape, direction='down'):
        try:
            assert self.shape[0]*self.shape[1]==newShape[0]*newShape[1] and type(newShape[1])==type(newShape[0])==int
        except:
            raise Error('Number of entries in new shape must be '+str(self.shape[0]*self.shape[1])+'.')

        mData=[[0]*newShape[1] for _ in xrange(newShape[0])]
        if str.lower(direction)=='down':
            entries=self.iterDown()
            for j in xrange(newShape[1]):
                for i in xrange(newShape[0]):
                    mData[i][j]=entries.next()
            return matrix(data=mData)
        elif str.lower(direction)=='right':
            entries=self.iterRight()
            for i in xrange(newShape[0]):
                for j in xrange(newShape[1]):
                    mData[i][j]=entries.next()
            return matrix(data=mData)
        else:
            raise Error('Direction must be eitehr \'down\' or \'right\'.')

    def sort(self, cmp=None, key=None, reverse=False):
        l=[x for x in self]
        l.sort(cmp=cmp,key=key,reverse=reverse)
        liter=l.__iter__()
        for i in xrange(newShape[0]):
            for j in xrange(newShape[1]):
                mData[i][j]=liter.next()
        self.data=mData
        self.horizontal=True

    def __str__(self, forCopy=False):
        if max(self.shape)==0:
            return '[[ ]]'
        digits = self.digitsToShow
        fmtS = '.'+str(digits)+'g'
        spaceForNumber = self.__getSpaceForNumberNeeded(fmtS)
        spaceSep=2
        fn = self.__getNumberString

        shp1=self.shape[1]
        fn2 = self.__getRowString
        if forCopy:
            nl=',\n '
        else:
            nl='\n '
        return '['+nl.join(fn2(i, fmtS, fn,shp1, spaceForNumber, spaceSep, forCopy) for i in xrange(self.shape[0]))+']'


    def __getRowString(self, rowN, fmtS, fn,shp1, spaceForNumber, spaceSep, forCopy):
        if forCopy:
            nl=','
        else:
            nl=''
        return '['+(nl+' '*spaceSep).join(fn(rowN, j, fmtS, spaceForNumber) for j in xrange(shp1))+']'
    def __getNumberString(self,rowN, colN, fmtS, spaceForNumber):

        num = self[rowN,colN]
        numStr = format(num,fmtS)
        return numStr+' '*(spaceForNumber-len(numStr))
    def __getSpaceForNumberNeeded(self, fmtS):
        mx = max(max(abs(y) for y in v) for v in self.data )

        try:
            mn = min(min(abs(y) for y in v if abs(y)>0) for v in self.data )
        except:
            mn=mx
        return max(len(format(mx,fmtS)),len(format(mn,fmtS)))+1#add one for minus sign

    def __repr__(self):
        return self.__str__()

    def __contains__(self, y):
        return any(y==x for x in self)

    def __setattr__(self, name, value):
        return super(matrix, self).__setattr__(name, value)
    def __len__(self):
        return max(self.shape)
    def __iter__(self):
        return self.iterRight()
    def iterDown(self):
        selfshape=self.shape
        for j in xrange(selfshape[1]):
            for i in xrange(selfshape[0]):
                yield self[i,j]
    def iterRight(self):
        selfshape=self.shape
        for i in xrange(selfshape[0]):
            for j in xrange(selfshape[1]):
                yield self[i,j]
    def toNumpyMatrix(self):
        return np.matrix(self.toLists())
    def solve(self,b):
        if not self.shape[0]==self.shape[1]:
            raise Error('Matrix is not square, use lstsq function')
        Anp = np.matrix(self.toLists())
        if type(b)==matrix:
            bnp = np.matrix(b.toLists())
        else:
            bnp = np.matrix(b)
        return matrix(data=np.linalg.solve(Anp,bnp).tolist())
    def lstsq(self,b):
        Anp = np.matrix(self.toLists())
        if type(b)==matrix:
            bnp = np.matrix(b.toLists())
        else:
            bnp = np.matrix(b)
        return matrix(data=np.linalg.lstsq(Anp,bnp)[0].tolist())

        #solves Ax=b for x using numpy!
    def nullspace(self, eps=1e-12):
        A = np.matrix(self.toLists())
        return matrix(data=nullSpace(A, eps=1e-12).tolist())
    def svd(self):
        #gives SVD using numpy!
        assert False
    def l1norm(self):
        #gives l1norm
        sz=self.shape
        return max( sum(abs(self[i,j]) for i in xrange(sz[0])) for j in xrange(sz[1]))

    def l2norm(self):
        #gives l1norm
        sz=self.shape
        if sz[0]==1:
            return math.sqrt(sum(self[0,j]**2 for j in xrange(sz[1]) ))
        elif sz[1]==1:
            return math.sqrt(sum(self[i,0]**2 for i in xrange(sz[0]) ))
        else:
            Anp = np.matrix(self.toLists())
            return numpy.linalg.norm(np.matrix(Anp))
    def fronorm(self):
        return math.sqrt(sum(x**2 for x in self))
    def linfnorm(self):
        #gives linfnorm
        sz=self.shape
        return max( sum(abs(self[i,j]) for j in xrange(sz[1])) for i in xrange(sz[0]))
    def applyScalarFun(self,fun):
        newSelf = self.copy()
        for i in xrange(self.shape[0]*self.shape[1]):
            newSelf[i] = fun(newSelf[i])
        return newSelf

def trueEquality(obj1,obj2, objId=(), stack=(), printStack=False, depth=0):
    if id(obj1)==id(obj2):
        return True
    if depth==0:#only do the first time!
        objId={id(obj1):set([id(obj2)])}
        stack=[]
    if depth>0:
        try:
            if id(obj2) in objId[id(obj1)]:
                return True #we already compared these two!
        except:
            pass

    hasDic=True
    hasSlots=True
    if not type(obj1)==type(obj2):
        if printStack:
            stack.append('  <- type difference')
            if depth==0:
                print 'obj'+''.joing(reversed(stack))
        return False
    try:
        obj1.__dict__
    except:
        hasDic=False
    try:
        obj1.__slots__
    except:
        hasSlots=False
    if hasSlots:
        if not obj1.__dict__=={}:
            raise Error('Obj1 has slots and non-empty __dict__, obj1 = '+str(id(obj1)) + ', obj2 = '+str(id(obj2)) )
        if not obj2.__dict__ == {}:
            if printStack:
                stack.append('.__dict__  <-- Obj2 has a non-empty __dict__, but obj1 has __slots__ , obj1 = '+str(id(obj1)) + ', obj2 = '+str(id(obj2)))
                if depth==0:
                    print 'Obj'+''.join(reversed(stack))
            return False
        if not obj1.__slots__==obj2.__slots__:
            if printStack:
                stack.append('.__slots__')
                if depth==0:
                    print 'Obj'+''.join(reversed(stack))
            return False
        try:
            objId[id(obj1)].add(id(obj2))
        except:
            objId[id(obj1)] = set([id(obj2)])
        for i in obj1.__slots__:
            if not trueEquality(getattr(obj1,i), getattr(obj2,i), stack=stack, objId=objId, depth=depth+1, printStack=printStack):
                if printStack:
                    stack.append('.'+i)
                    if depth==0:
                        print 'Obj'+''.join(reversed(stack))
                return False

    else:
        if hasDic:
            if not len(obj1.__dict__)==len(obj2.__dict__):
                if printStack:
                    stack.append('.__dict__  <-- lengths not equal')
                    if depth==0:
                        print 'Obj'+''.join(reversed(stack))
                return False
            j=obj2.__dict__
            for i in obj1.__dict__.iterkeys():
                if not i in j:
                    if printStack:
                        stack.append('.__dict__  <-- missing key '+str(i))
                        if depth==0:
                            print 'Obj'+''.join(reversed(stack))
                    return False

            try:
                objId[id(obj1)].add(id(obj2))
            except:
                objId[id(obj1)] = set([id(obj2)])

            for ik,i in obj1.__dict__.iteritems():
                jj=j[ik]
                if not trueEquality(i, jj, stack=stack, objId=objId, depth=depth+1, printStack=printStack):
                    if printStack:
                        if type(ik)==str:
                            ik='\''+ik+'\''
                        stack.append('.'+str(ik))
                        if depth==0:
                            print 'Obj'+''.join(reversed(stack))
                    return False
        else:
            if type(obj1)==dict:
                if not len(obj1)==len(obj2):
                    if printStack:
                        stack.append('  <-- different lengths')
                        if depth==0:
                            print 'Obj'+''.join(reversed(stack))
                    return False
                for i in obj1:
                    if not i in obj2:
                        if printStack:
                            stack.append('  <-- keys are not the same')
                            if depth==0:
                                print 'Obj'+''.join(reversed(stack))
                        return False
                try:
                    objId[id(obj1)].add(id(obj2))
                except:
                    objId[id(obj1)] = set([id(obj2)])
                for ik,i in obj1.iteritems():
                    if not trueEquality(i,obj2[ik], stack=stack, objId=objId, depth=depth+1, printStack=printStack):
                        if printStack:
                            if type(ik)==str:
                                ik='\''+ik+'\''
                            stack.append('['+str(ik)+']')
                            if depth==0:
                                print 'Obj'+''.join(reversed(stack))
                        return False
            elif type(obj1)==tuple or type(obj1)==list:
                if not len(obj1)==len(obj2):
                    if printStack:
                        stack.append('  <-- different lengths')
                        if depth==0:
                            print 'Obj'+''.join(reversed(stack))
                    return False
                try:
                    objId[id(obj1)].add(id(obj2))
                except:
                    objId[id(obj1)] = set([id(obj2)])
                for i in xrange(len(obj1)):
                    if not trueEquality(obj1[i],obj2[i], stack=stack, objId=objId, depth=depth+1, printStack=printStack):
                        if printStack:
                            stack.append('['+str(i)+']')
                            if depth==0:
                                print 'Obj'+''.join(reversed(stack))
                        return False
            else:
                if not obj1==obj2:
                    if printStack:
                        stack.append('  <-- not equal according to __eq__ fn')
                        if depth==0:
                            print 'Obj'+''.join(reversed(stack))
                    return False
    return True
