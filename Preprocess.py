import matplotlib.pyplot as plt
import math
import copy
import numpy as np
import csv
def main(f1,f2):
    index=0
    newind6=0
    newtempind=0
    largest=-10000
    smallest=10000
    smallest2=10000
    largest2=-10000
    maxcorrelation=0
    maxcorrelation2=0
    maxcorrelation3=0
    maxcorrelation4=0
    arraytotal=0
    totalpowerdifference=0
    tempnewind10=0
    itemindex=0
    finalcount=0
    tempoutlierindex=0
    finalindex=0
    tempcount=0
    anotempcount=0
    start=-1
    start2=-1
    a=-1
    mylist=[]
    peptide=[]
    anoindex=[]
    controlindex=[]
    controlitem=[]
    newcontrolitem=[]
    averagecontrollist=[]
    newaveragecontrollist=[]
    totalaveragecontrollist=[]
    samplename=[]
    typelist=[]
    array=[]
    temparray=[]
    finalarray=[]
    finalresult=[]
    fullwidthlist=[]
    column=[]
    column2=[]
    arraymedian=[]
    arraymedian2=[]
    arraymedian3=[]
    arraymedian4=[]
    combinedmedian=[]
    numcountzero=[]
    numcountmax=[]
    negativecontrol=[]
    normalizednegativecontrol=[]
    serumcontrol=[]
    normalizedserumcontrol=[]
    serumsample=[]
    normalizedserumsample=[]
    nonambiguous=[]
    ambiguous=[]
    averagelinkeronlysignal=[]
    arrayaveragelist=[]
    anoarrayaveragelist=[]
    variationcoefficientlist=[]
    anovariationcoefficientlist=[]
    anovariationcoefficientlist2=[]
    anoreference=[]
    azcvlist=[]
    cacvlist=[]
    bothcvlist=[]
    correlationlist=[]
    correlationlist2=[]
    correlationlist3=[]
    newcorrelationlist=[]
    newcorrelationlist2=[]
    newcorrelationlist3=[]
    newcorrelationlist4=[]
    item=""
    data=f1.read()
    data2=f2.read()
    totalindex=findindex(data2)
    getindex=totalindex[0]
    tempgetindex=copy.copy(getindex)
    azindex=totalindex[1]
    caindex=totalindex[2]
    bothindex=totalindex[3]
    reference=totalindex[4]
    tempoutlierlist=[]
    tempoutlierlist2=[]
    tempoutlierlistsumlist=[]
    totaltempoutlierlist=[]
    totaltempoutlierlist2=[]
    outliermedianlist=[]
    newcombinedmedian=[]
    normalizedoutlierlist=[]
    totalnormalizedoutlierlist=[]
    newoutliermedianlist=[]
    newncmedianlist=[]
    newscmedianlist=[]
    newssmedianlist=[]
    newcombinedmedianlist=[]
    meanlist=[]
    newncmeanlist=[]
    newscmeanlist=[]
    newssmeanlist=[]
    newoutliermeanlist=[]
    combinedmeanlist=[]
    newnonoutliermedian=[]
    scvariationcoefficientlist=[]
    ncvariationcoefficientlist=[]
    nonoutliersscvlist=[]
    outliersscvlist=[]
    newcombinedcvlist=[]
    anofinalresult=[]
    threshold=math.log(100,10)
    threshold2=math.log(65635,10)
    
    for ind2 in azindex:
        anoindex.append(ind2)

    for ind3 in caindex:
        anoindex.append(ind3)

    for ind4 in bothindex:
        anoindex.append(ind4)

    anoindex.sort()
    print("length of anoindex "+str(len(anoindex)))

    while index<len(data):
        while data[index]!='\t'and data[index]!='\n' and index<len(data)-1:
            item=item+data[index]
            index=index+1
        if index==len(data)-1:
            item=item+data[len(data)-1]
        if "XY" in item:
            subnewitem=item[1:len(item)-1]
            peptide.append(subnewitem)
        elif "13-1000" in item:
            samplename.append(item)
        elif "control" in item or "library" in item:
            typelist.append(item)
        elif "control" not in item and "library" not in item and "ID" not in item and "_" not in item and "type" not in item:
            tempitem=float(item)
            logtransform=math.log(tempitem+100,10)
            if logtransform<threshold:
                logtransform=threshold
            elif logtransform>threshold2:
                logtransform=threshold2
            peptidelastindex=peptide[len(peptide)-1]
            if "C" in peptidelastindex:
                mylist.append(str(logtransform)+"C")
            else:
                mylist.append(str(logtransform))
        index=index+1
        item=""
    length=len(peptide)
    end=length-1
    end2=length-1
    b=length-1
    samplenumber=len(samplename)
    print("samplename[0] "+str(samplename[0]))
   
    for newindex in range(0,len(peptide)):
        if "C" in peptide[newindex]:
            controlindex.append(newindex)
    
    for newind2 in range(0,samplenumber):
        newind3=newind2
        while newind3<len(mylist):
            newitem2=mylist[newind3]
            array.append(newitem2)
            newind3=newind3+samplenumber
    print("array[len(array)-1] "+str(array[len(array)-1]))
    print("length of array "+str(len(array)))

    temparray=copy.copy(array)
    mylist=[]
    anomylist=[]

    for tempind in range(0,samplenumber):
        while newtempind>a and newtempind<b+1:
            if "C" in array[newtempind]:
                array[newtempind]=float(array[newtempind][0:len(array[newtempind])-1])
            else:
                array[newtempind]=float(array[newtempind])
                column2.append(array[newtempind])
            originaltempitem=math.pow(10,array[newtempind])-100
            column.append(originaltempitem)
            newtempind=newtempind+1
        anostandarddeviation=np.std(column)
        anoarrayaverage=np.mean(column)
        anoarrayaveragelist.append(anoarrayaverage)
        anotempvariationcoefficient=anostandarddeviation/anoarrayaverage
        anovariationcoefficientlist.append(anotempvariationcoefficient)
        if tempind in anoindex:
            anoreference.append(reference[tempind])
            serumcontrol.append(column)
        elif tempind in getindex:
            negativecontrol.append(column)
   
        countzero=column2.count(threshold)
        countmax=column2.count(4.81713548958)
        column2.sort()
        percentileindex=int(round(len(column2)/float(4)*3,0))
        fullwidth=column2[len(column2)-1]-column2[percentileindex]
        fullwidthlist.append(fullwidth)
        median=np.median(column)
        median2=np.median(column2)
        arraymedian.append(median)
        arraymedian2.append(median2)
        numcountzero.append(countzero/float(len(column2)))
        numcountmax.append(countmax/float(len(column2)))
        
        a=a+length
        b=b+length
        column=[]
        column2=[]
        controlitem=[]

    outlier=findoutlier(getindex,azindex,caindex,bothindex,anoarrayaveragelist,anovariationcoefficientlist,arraymedian,samplenumber)
    outlier.sort()

    for outlierlist in range(0,len(outlier)):
        startoutlierindex=outlier[outlierlist]*length
        tempoutlierlist.append(outlierlist)
        for outlierlist2 in range(startoutlierindex,startoutlierindex+length):
            tempoutlierlist.append(array[outlierlist2])
            tempoutlieritem=int(round(math.pow(10,array[outlierlist2])-100))
            tempoutlierlist2.append(tempoutlieritem)
        tempoutlierlistsum=sum(tempoutlierlist[1:len(tempoutlierlist)-1])
        tempoutlierlistsumlist.append(tempoutlierlistsum)
        totaltempoutlierlist.append(tempoutlierlist)
        totaltempoutlierlist2.append(tempoutlierlist2)
        tempoutlierlist=[]
        tempoutlierlist2=[]

    f3=open('E:/health_tell/Expt_864/outliers.csv', 'wb')
    spamwriter = csv.writer(f3)
    for y in range(0,len(totaltempoutlierlist[0])):
        spamwriter.writerow([str(x[y]) for x in totaltempoutlierlist])
    f3.close()
    totaltempoutlierlist=[]

    normfactors=[1.0036555,0.9995349,0.9989067,0.9999503,1.0020235,1.0038480,0.9994047,0.9989573,1.0002816,1.0032265,0.9994666,0.9932539,0.9990770,0.9988909,0.9995695]
    
    for outlierlist3 in range(0,len(outlier)):
        startoutlierindex2=outlier[outlierlist3]*length
        for outlierlist4 in range(startoutlierindex2,startoutlierindex2+length):
            array[outlierlist4]=array[outlierlist4]/normfactors[outlierlist3]
            normalizedoutlierlist.append(array[outlierlist4])
        totalnormalizedoutlierlist.append(normalizedoutlierlist)
        newoutliermedian=np.median(normalizedoutlierlist)
        newoutliermean=np.mean(normalizedoutlierlist)
        newoutliermedianlist.append(newoutliermedian)
        newoutliermeanlist.append(newoutliermean)
        normalizedoutlierlist=[]

    for ncindex in range(0,len(tempgetindex)):
        for outlierindex in range(0,len(outlier)):
            if outlier[outlierindex]<tempgetindex[ncindex]:
                tempcount=tempcount+1
        tempgetindex[ncindex]=tempgetindex[ncindex]-tempcount
        tempcount=0

    for outlierindex in outlier:
        fullwidthlist[outlierindex]=-1

    while -1 in fullwidthlist:
        fullwidthlist.remove(-1)

    for tempindex in range(0,len(outlier)):
        anotempitem=outlier[tempindex]
        outliermedianlist.append(arraymedian[anotempitem])
        arraymedian[anotempitem]=-1
    
    for tempindex2 in arraymedian:
        if tempindex2!=-1:
            arraymedian3.append(tempindex2)
    
    combinedmedian.append(arraymedian3)
    newcombinedmedian.append(outliermedianlist)
    newcombinedmedian.append(newoutliermedianlist)

    for anotempind in range(0,len(anovariationcoefficientlist)):
        if anotempind in anoindex:
            anovariationcoefficientlist2.append(anovariationcoefficientlist[anotempind])
    
    for tempind2 in range(0,len(serumcontrol)-1):
        for tempind3 in range(tempind2+1,len(serumcontrol)):
            correlation=np.corrcoef(serumcontrol[tempind2],serumcontrol[tempind3])[0,1]
            correlationlist.append(correlation)
    serumcontrol=[]

    for tempind4 in range(0,len(negativecontrol)-1):
        for tempind5 in range(tempind4+1,len(negativecontrol)):
            correlation2=np.corrcoef(negativecontrol[tempind4],negativecontrol[tempind5])[0,1]
            correlationlist2.append(correlation2)
    negativecontrol=[]

    for tempind6 in range(0,len(totaltempoutlierlist2)-1):
        for tempind7 in range(tempind6+1,len(totaltempoutlierlist2)):
            correlation3=np.corrcoef(totaltempoutlierlist2[tempind6],totaltempoutlierlist2[tempind7])[0,1]
            correlationlist3.append(correlation3)
    totaltempoutlierlist2=[]

    fig=plt.boxplot(numcountmax)
    plt.xticks([1], [' '])
    plt.title("fraction of library peptides at scanner maximum of 65,535 before normalization")
    plt.show()
    plt.close()
    plt.clf()
    del fig

    fig=plt.boxplot(numcountzero)
    plt.xticks([1], [' '])
    plt.title("fraction of library peptides with zero signal before normalization")
    plt.show()
    plt.close()
    plt.clf()
    del fig

    fig2=plt.boxplot(fullwidthlist)
    plt.xticks([1], [' '])
    plt.title("full width at one quarter of the maximum of the transformed intensity distribution")
    plt.show()
    plt.close()
    plt.clf()
    del fig2
    
    numcountzero=[]
    numcountmax=[]
    fullwidthlist=[]

    for newind4 in range(0,samplenumber):
        getlibrarymedian=arraymedian2[newind4]
        while newind6>start2 and newind6<end2+1:
            if newind4 not in outlier:
                modifieditem=array[newind6]-getlibrarymedian
                array[newind6]=modifieditem
            newind6=newind6+1
        start2=start2+length
        end2=end2+length
    
    for newind8 in range(0,samplenumber):
        start3=newind8*length
        if newind8 in getindex:
            for newind9 in range(start3, start3+length):
                if "C" in temparray[newind9]:
                    tempitem2=temparray[newind9]
                    newtempitem2=float(tempitem2[0:len(tempitem2)-1])
                else:
                    newtempitem2=float(temparray[newind9])
                ambiguous.append(newtempitem2)
                temparray[newind9]=-10000
    ambiguousmedian=np.median(ambiguous)
    ambiguous=[]
        
    temparray.sort()
    for newind9 in range(0,len(temparray)):
        if temparray[newind9]==-10000:
            tempstart=newind9
    newtemparray=temparray[tempstart+1:len(temparray)]
    
    if len(newtemparray)%2==0:
        item1=newtemparray[len(newtemparray)/2-1]
        item2=newtemparray[len(newtemparray)/2]
        total=item1+item2
        nonambiguousmedian=float(total)/2
    else:
        nonambiguousmedian=newtemparray[len(newtemparray)/2]
    newtemparray=[]
    
    for newind10 in range(0,samplenumber):
        subnewarraystart=newind10*length
        if newind10 in getindex:
            for newind10j in range(subnewarraystart, subnewarraystart+length): 
                array[newind10j]=float(array[newind10j])+float(ambiguousmedian)
                if array[newind10j]<threshold:
                    array[newind10j]=threshold
                elif array[newind10j]>threshold2:
                    array[newind10j]=threshold2
        elif newind10 not in outlier:  
            for newind10j in range(subnewarraystart, subnewarraystart+length):
                array[newind10j]=float(array[newind10j])+float(nonambiguousmedian)
                if array[newind10j]<threshold:
                    array[newind10j]=threshold
                elif array[newind10j]>threshold2:
                    array[newind10j]=threshold2
        column3=array[subnewarraystart:subnewarraystart+length-1]
        arrayaverage=np.mean(column3)
        arrayaveragelist.append(arrayaverage)
        standarddeviation=np.std(column3)
        tempvariationcoefficient=standarddeviation/arrayaverage
        variationcoefficientlist.append(tempvariationcoefficient)
        median3=np.median(column3)
        arraymedian4.append(median3)
        mean=np.mean(column3)
        arraymedian.append(median)
        meanlist.append(mean)

    for newind9j in range(0,len(arraymedian4)):
        if newind9j in getindex:
            newncmedianlist.append(arraymedian4[newind9j])
            newncmeanlist.append(meanlist[newind9j])
            
        elif newind9j not in outlier:
            if newind9j in anoindex:
                newscmedianlist.append(arraymedian4[newind9j])
                newscmeanlist.append(meanlist[newind9j])
            else:
                newssmedianlist.append(arraymedian4[newind9j])
                newssmeanlist.append(meanlist[newind9j])

    newcombinedmedianlist.append(newncmedianlist)
    newcombinedmedianlist.append(newscmedianlist)
    newcombinedmedianlist.append(newssmedianlist)
    newcombinedmedianlist.append(newoutliermedianlist)

    tempfig2=plt.boxplot(newcombinedmedianlist)
    plt.xticks([1,2,3,4], ['nc', 'sc', 'non-outliers of ss', 'outliers of ss'])
    plt.title('distribution of median of different groups after normalization')
    plt.savefig('tempfig2.png')
    plt.show()
    plt.close()
    plt.clf()
    del tempfig2

    newncmedianlist=[]
    newscmedianlist=[]
    newssmedianlist=[]
    newcombinedmedianlist=[]

    for tempindex3 in range(0,len(outlier)):
        anotempitem2=outlier[tempindex3]
        arraymedian4[anotempitem2]=-1

    for tempindex4 in range(0,len(arraymedian4)):
        if arraymedian4[tempindex4]!=-1:
            newnonoutliermedian.append(arraymedian4[tempindex4])
        
    temparray=[]
    combinedmedian.append(newnonoutliermedian)
    fig3=plt.boxplot(combinedmedian)
    plt.xticks([1, 2], ['before normalization', 'after normalization'])
    plt.title('median intensity of non-outliers')
    plt.show()
    plt.close()
    plt.clf()
    del fig3

    arraymedian3=[]
    arraymedian4=[]
    newnonoutliermedian=[]
    combinedmedian=[]

    combinedmeanlist.append(newncmeanlist)
    combinedmeanlist.append(newscmeanlist)
    combinedmeanlist.append(newssmeanlist)
    combinedmeanlist.append(newoutliermeanlist)

    tempfig3=plt.boxplot(combinedmeanlist)
    plt.xticks([1,2,3,4], ['nc', 'sc', 'non-outliers of ss', 'outliers of ss'])
    plt.title('distribution of average of different groups after normalization')
    plt.savefig('tempfig3.png')
    plt.show()
    plt.close()
    plt.clf()
    del tempfig3
    newncmeanlist=[]
    newscmeanlist=[]
    newssmeanlist=[]
    newoutliermeanlist=[]
    combinedmeanlist=[]

    for newind7 in range(0,len(variationcoefficientlist)):
        if newind7 in anoindex:
            scvariationcoefficientlist.append(variationcoefficientlist[newind7])
            if "AZ" in reference[newind7]:
                azcvlist.append(variationcoefficientlist[newind7])
            elif "CA" in reference[newind7]:
                cacvlist.append(variationcoefficientlist[newind7])
            else:
                bothcvlist.append(variationcoefficientlist[newind7])
        elif newind7 in getindex:
            ncvariationcoefficientlist.append(variationcoefficientlist[newind7])
        elif newind7 in outlier:
            nonoutliersscvlist.append(variationcoefficientlist[newind7])
        else:
            outliersscvlist.append(variationcoefficientlist[newind7])

    newcombinedcvlist.append(ncvariationcoefficientlist)
    newcombinedcvlist.append(scvariationcoefficientlist)
    newcombinedcvlist.append(nonoutliersscvlist)
    newcombinedcvlist.append(outliersscvlist)

    tempfig4=plt.boxplot(newcombinedcvlist)
    plt.xticks([1,2,3,4], ['nc', 'sc', 'non-outliers of ss', 'outliers of ss'])
    plt.title('distribution of cv of different groups after normalization')
    plt.savefig('tempfig4.png')
    plt.show()
    plt.close()
    plt.clf()
    del tempfig4

    ncvariationcoefficientlist=[]
    nonoutliersscvlist=[]
    outliersscvlist=[]
    newcombinedcvlist=[]
    
    plt.hist(azcvlist,label='PooledFlu2016_AZ')
    plt.hist(cacvlist,label='PooledFlu2016_CA')
    plt.hist(bothcvlist,label='PooledFlu2016_Both')
    plt.legend(loc='upper right')
    plt.ylabel('frequency')
    plt.title("distribution of coefficient of variation of serum controls after normalization")
    plt.show()
    plt.close()
    plt.clf()
    azcvlist=[]
    cacvlist=[]
    bothcvlist=[]

    n=len(scvariationcoefficientlist)
    fig, ax=plt.subplots()
    index=np.arange(n)
    bar_width=0.35
    opacity=0.8
    rects1 = plt.bar(index, anovariationcoefficientlist2, bar_width,alpha=opacity,color='b',label='before normalization')
    rects2 = plt.bar(index+bar_width, scvariationcoefficientlist, bar_width,alpha=opacity,color='g',label='after normalization')
    plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14], ['PooledFlu2016_AZ', 'PooledFlu2016_CA', 'PooledFlu2016_AZ', 'PooledFlu2016_CA', 'PooledFlu2016_Both', 'PooledFlu2016_AZ', 'PooledFlu2016_CA', 'PooledFlu2016_AZ', 'PooledFlu2016_CA', 'PooledFlu2016_Both', 'PooledFlu2016_AZ', 'PooledFlu2016_CA', 'PooledFlu2016_AZ', 'PooledFlu2016_CA', 'PooledFlu2016_Both'], rotation ='vertical')
    plt.xlabel('sample')
    plt.ylabel('intra-array coefficient of variation of serum control peptides')
    for rect in rects1:
        h=rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2.,1.01*h,"%.3f" % round(h,3),ha='center',va='bottom',rotation ='vertical')
    for rect in rects2:
        h=rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2.,1.01*h,"%.3f" % round(h,3),ha='center',va='bottom',rotation ='vertical')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.close()
    plt.clf()

    arraymedian2=[]
    scvariationcoefficientlist=[]
    anovariationcoefficientlist2=[]
    
    for tempnewind10 in range(0,samplenumber):
        if tempnewind10 not in outlier:
            for anotempnewind10 in range(0,len(controlindex)):
                anotempnewind10j=tempnewind10*length
                controlitem.append(array[controlindex[anotempnewind10]+anotempnewind10j])
            averagecontrol=np.mean(controlitem)
            averagecontrollist.append(averagecontrol)
            controlitem=[]
        else:
            for anotempnewind10 in range(0,len(controlindex)):
                anotempnewind10j=tempnewind10*length
                newcontrolitem.append(array[controlindex[anotempnewind10]+anotempnewind10j])
            newaveragecontrol=np.mean(newcontrolitem)
            newaveragecontrollist.append(newaveragecontrol)
            newcontrolitem=[]
    totalaveragecontrollist.append(averagecontrollist)
    totalaveragecontrollist.append(newaveragecontrollist)

    tempfig=plt.boxplot(totalaveragecontrollist)
    plt.xticks([1,2], ['non-outlers', 'outliers'])
    plt.title('average linker-only signals after normalization')
    plt.show()
    plt.close()
    plt.clf()
    del tempfig
    averagecontrollist=[]
    newaveragecontrollist=[]
    totalaveragecontrollist=[]
    

    for newind11 in range(0,samplenumber):
        tempstart2=newind11*length
        if newind11 in anoindex:
            normalizedserumcontrol.append(array[tempstart2:tempstart2+length])
        elif newind11 in getindex:
            normalizednegativecontrol.append(array[tempstart2:tempstart2+length])
        elif newind11 not in anoindex and newind11 not in getindex and newind11 not in outlier:
            normalizedserumsample.append(array[tempstart2:tempstart2+length])

    for newind12 in range(0,len(normalizedserumcontrol)-1):
        for newind13 in range(newind12+1,len(normalizedserumcontrol)):
            newcorrelation=np.corrcoef(normalizedserumcontrol[newind12],normalizedserumcontrol[newind13])[0,1]
            newcorrelationlist.append(newcorrelation)

    for newind14 in range(0,len(normalizednegativecontrol)-1):
        for newind15 in range(newind14+1,len(normalizednegativecontrol)):
            newcorrelation2=np.corrcoef(normalizednegativecontrol[newind14],normalizednegativecontrol[newind15])[0,1]
            newcorrelationlist2.append(newcorrelation2)

    for newind16 in range(0,len(normalizedserumsample)-1):
        for newind17 in range(newind16+1,len(normalizedserumsample)):
            newcorrelation3=np.corrcoef(normalizedserumsample[newind16],normalizedserumsample[newind17])[0,1]
            newcorrelationlist3.append(newcorrelation3)

    for newind18 in range(0,len(totalnormalizedoutlierlist)-1):
        for newind19 in range(newind18+1,len(totalnormalizedoutlierlist)):
            newcorrelation4=np.corrcoef(totalnormalizedoutlierlist[newind18],totalnormalizedoutlierlist[newind19])[0,1]
            newcorrelationlist4.append(newcorrelation4)

    
    n=len(correlationlist)
    fig, ax=plt.subplots()
    index=np.arange(n)
    bar_width=0.35
    opacity=0.8
    rects1 = plt.bar(index, correlationlist, bar_width,alpha=opacity,color='b',label='before normalization')
    rects2 = plt.bar(index+bar_width, newcorrelationlist, bar_width,alpha=opacity,color='g',label='after normalization')
    plt.title('Pearson correlations of serum controls')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.close()
    plt.clf()

    n=len(correlationlist2)
    fig, ax=plt.subplots()
    index=np.arange(n)
    bar_width=0.35
    opacity=0.8
    newrects1 = plt.bar(index, correlationlist2, bar_width,alpha=opacity,color='b',label='before normalization')
    newrects2 = plt.bar(index+bar_width, newcorrelationlist2, bar_width,alpha=opacity,color='g',label='after normalization')
    plt.title('Pearson correlations of negative controls')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.close()
    plt.clf()

    fig6=plt.boxplot(newcorrelationlist3)
    plt.xticks([1], [' '])
    plt.title('Pearson correlations of non-outliers of serum samples after normalization')
    plt.show()
    plt.close()
    plt.clf()
    del fig6

    normalizednegativecontrol=[]
    normalizedserumcontrol=[]
    normalizedserumsample=[]
    newcorrelationlist=[]
    newcorrelationlist2=[]
    newcorrelationlist3=[]

    fig7=plt.boxplot(newcombinedmedian)
    plt.xticks([1, 2], ['before normalization', 'after normalization'])
    plt.title('median intensity of outliers of serum samples')
    plt.show()
    plt.close()
    plt.clf()
    del fig7

    n=len(correlationlist3)
    fig, ax=plt.subplots()
    index=np.arange(n)
    bar_width=0.35
    opacity=0.8
    rects1 = plt.bar(index, correlationlist3, bar_width,alpha=opacity,color='b',label='before normalization')
    rects2 = plt.bar(index+bar_width, newcorrelationlist4, bar_width,alpha=opacity,color='g',label='after normalization')
    plt.title('Pearson correlations of outliers of serum samples')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()
    plt.close()
    plt.clf()
        
    f4=open('E:/health_tell/Expt_864/result.csv', 'wb')
    spamwriter=csv.writer(f4)
    peptide.insert(0,"peptideID")
    typelist.insert(0,"type")
    anofinalresult.append(peptide)
    anofinalresult.append(typelist)
    while finalindex<len(array)-length+1:
        tempanofinalresult=array[finalindex:finalindex+length]
        tempanofinalresult.insert(0,samplename[anotempcount])
        anofinalresult.append(tempanofinalresult)
        finalindex=finalindex+length
        anotempcount=anotempcount+1
    array=[]
    for y in range(0,len(anofinalresult[0])):
        spamwriter.writerow([str(x[y]) for x in anofinalresult])
    
def findindex(data2):
    index2=0
    mylist2=[]
    reference=[]
    reference2=[]
    result=[]
    result2=[]
    result3=[]
    result4=[]
    result5=[]
    totalresult=[]
    item2=""
    while index2<len(data2):
        while data2[index2]!='\t'and data2[index2]!='\n' and index2<len(data2)-1:
            item2=item2+data2[index2]
            index2=index2+1
        if index2==len(data2)-1:
            item2=item2+data2[len(data2)-1]
        subitem2=item2[1:len(item2)-1]
        mylist2.append(subitem2)
        index2=index2+1
        item2=""

    for index3 in range(0,len(mylist2)):
        if "13-1000" in mylist2[index3]:
            sample=mylist2[index3]
            sampletype=mylist2[index3+1]
            reference.append(sample)
            reference2.append(sampletype)

    for index4 in range(0,len(reference2)):
        if reference2[index4]=="Secondary Only (IgG)":
            result.append(index4)
        elif reference2[index4]=="PooledFlu2016_AZ":
            result2.append(index4)
        elif reference2[index4]=="PooledFlu2016_CA":
            result3.append(index4)
        elif reference2[index4]=="PooledFlu2016_Both":
            result4.append(index4)
        result5.append(reference2[index4])
    totalresult.append(result)
    totalresult.append(result2)
    totalresult.append(result3)
    totalresult.append(result4)
    totalresult.append(result5)
    return totalresult

def findoutlier(getindex,azindex,caindex,bothindex,averagelist,cvlist,medianlist,samplenumber):
    ssindex=[]
    nccvlist=[]
    azcvlist=[]
    cacvlist=[]
    bothcvlist=[]
    sscvlist=[]
    ncmedianlist=[]
    azmedianlist=[]
    camedianlist=[]
    bothmedianlist=[]
    ssmedianlist=[]
    ncaveragelist=[]
    azaveragelist=[]
    caaveragelist=[]
    bothaveragelist=[]
    ssaveragelist=[]
    outlier=[]
    totaloutlier=[]
    
      
    for ind in range(0,samplenumber):
        if ind not in getindex and ind not in azindex and ind not in caindex and ind not in bothindex:
            ssindex.append(ind)
    
    for newind4 in getindex:
        ncmedianlist.append(medianlist[newind4])
        nccvlist.append(cvlist[newind4])
        ncaveragelist.append(averagelist[newind4])

    for newind5 in azindex:
        azmedianlist.append(medianlist[newind5])
        azcvlist.append(cvlist[newind5])
        azaveragelist.append(averagelist[newind5])

    for newind6 in caindex:
        camedianlist.append(medianlist[newind6])
        cacvlist.append(cvlist[newind6])
        caaveragelist.append(averagelist[newind6])

    for newind7 in bothindex:
        bothmedianlist.append(medianlist[newind7])
        bothcvlist.append(cvlist[newind7])
        bothaveragelist.append(averagelist[newind7])
        
    for newind8 in range(0,len(medianlist)):
        if newind8 in ssindex:
            ssmedianlist.append(medianlist[newind8])
            sscvlist.append(cvlist[newind8])
            ssaveragelist.append(averagelist[newind8])

    fig2=plt.hist(ncmedianlist,width=0.5)
    plt.ylabel('frequency')
    plt.title("distribution of median intensity of negative controls before normalization")
    plt.show()
    plt.close()
    plt.clf()
    del fig2
    ncmedianlist=[]
    
    plt.hist(azmedianlist,label='PooledFlu2016_AZ')
    plt.hist(camedianlist,label='PooledFlu2016_CA')
    plt.hist(bothmedianlist,label='PooledFlu2016_Both')
    plt.legend(loc='upper right')
    plt.ylabel('frequency')
    plt.title("distribution of median intensity of serum controls before normalization")
    plt.show()
    plt.close()
    plt.clf()
    azmedianlist=[]
    camedianlist=[]
    bothmedianlist=[]

    findpercentile=np.percentile(ssmedianlist, 95)
    for newind9 in range(0,len(ssmedianlist)):
        if ssmedianlist[newind9]>findpercentile:
            totaloutlier.append(ssindex[newind9])
            
    fig4=plt.hist(ssmedianlist)
    plt.ylabel('frequency')
    plt.title("distribution of median intensity of serum samples before normalization")
    plt.show()
    plt.close()
    plt.clf()
    del fig4
    ssmedianlist=[]

    fig5=plt.hist(nccvlist)
    plt.ylabel('frequency')
    plt.title("distribution of coefficient of variation of negative controls before normalization")
    plt.show()
    plt.close()
    plt.clf()
    del fig5
    nccvlist=[]

    plt.hist(azcvlist,label='PooledFlu2016_AZ')
    plt.hist(cacvlist,label='PooledFlu2016_CA')
    plt.hist(bothcvlist,label='PooledFlu2016_Both')
    plt.legend(loc='upper right')
    plt.ylabel('frequency')
    plt.title("distribution of coefficient of variation of serum controls before normalization")
    plt.show()
    plt.close()
    plt.clf()
    azcvlist=[]
    cacvlist=[]
    bothcvlist=[]

    findpercentile2=np.percentile(sscvlist, 95)
    for newind10 in range(0,len(sscvlist)):
        if sscvlist[newind10]>findpercentile2:
            totaloutlier.append(ssindex[newind10])
    fig6=plt.hist(sscvlist)
    plt.ylabel('frequency')
    plt.title("distribution of coefficient of variation of serum samples before normalization")
    plt.show()
    plt.close()
    plt.clf()
    del fig6
    sscvlist=[]

    fig7=plt.hist(ncaveragelist,width=0.5)
    plt.ylabel('frequency')
    plt.title("distribution of average of negative controls before normalization")
    plt.show()
    plt.close()
    plt.clf()
    del fig7
    ncaveragelist=[]

    plt.hist(azaveragelist,label='PooledFlu2016_AZ')
    plt.hist(caaveragelist,label='PooledFlu2016_CA')
    plt.hist(bothaveragelist,label='PooledFlu2016_Both')
    plt.legend(loc='upper right')
    plt.ylabel('frequency')
    plt.title("distribution of average of serum controls before normalization")
    plt.show()
    plt.close()
    plt.clf()
    azaveragelist=[]
    caaveragelist=[]
    bothaveragelist=[]
    
    findpercentile3=np.percentile(ssaveragelist, 95)
    for newind11 in range(0,len(ssaveragelist)):
        if ssaveragelist[newind11]>findpercentile3 and ssindex[newind11] not in totaloutlier:
            totaloutlier.append(ssindex[newind11])
            
    fig8=plt.hist(ssaveragelist)
    plt.ylabel('frequency')
    plt.title("distribution of average of serum samples before normalization")
    plt.show()
    plt.close()
    plt.clf()
    del fig8
    ssaveragelist=[]

    return totaloutlier

   
    
