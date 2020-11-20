import math
def main(f1,f2,f3,f4,f5):
    index=0
    index2=0
    index3=0
    index7=0
    index8=0
    index12=0
    ind3=0
    ind5=0
    ind6=0
    ind10=0
    countzero=0
    countmax=0
    tempcount=0
    countserumcontrol=0
    rowtotal=0
    totalsquareitem=0
    largest=-10000
    largest2=-10000
    start=-1
    start2=-1
    start6=-1
    samplenumber=38
    mylist=[]
    anomylist=[]
    mylist2=[]
    reference=[]
    totallist=[]
    newtotallist=[]
    peptide=[]
    array=[]
    subarray=[]
    arrayitem=[]
    arraymedian=[]
    arraynumzero=[]
    arraynummax=[]
    control=[]
    normalizedcontrol=[]
    library=[]
    controlindex=[]
    totalcontrol=[]
    totalnormalizedcontrol=[]
    rowtotallist=[]
    item1list=[]
    totallibrary=[]
    variationcoefficientlist=[]
    controlpearsoncorrelation=[]
    normalizedcontrolpearsoncorrelation=[]
    controlpearsoncorrelationindex=[]
    normalizedcontrolpearsoncorrelationindex=[]
    normalizedambiguous=[]
    normalizednonambiguous=[]
    normalizedambiguousnonambiguous=[]
    newnormalizedambiguousnonambiguous=[]
    item=""
    item2=""
    item3=""
    item4=""
    item5=""
    data=f1.read()
    data2=f2.read()
    data3=f3.read()
    data4=f4.read()
    data5=f5.read()
    
    while index<len(data):
        while data[index]!='\t'and data[index]!='\n' and index<len(data)-1:
            item=item+data[index]
            index=index+1
        if index==len(data)-1:
            item=item+data[len(data)-1]
        mylist.append(item)
        index=index+1
        item=""

    while index2<len(data2):
        while data2[index2]!=' 'and index2<len(data2)-1:
            item2=item2+data2[index2]
            index2=index2+1
        if index2==len(data2)-1:
            item2=item2+data2[len(data2)-1]
        item2=float(item2)
        anomylist.append(item2)
        index2=index2+1
        item2=""    
    
    for ind in range(0,len(mylist)):
        newitem=mylist[ind]
        if "XY" in newitem:
            subnewitem=newitem[1:len(newitem)-1]
            peptide.append(subnewitem)
            if "C" in newitem:
                countserumcontrol=countserumcontrol+1
                newitemindex=peptide.index(subnewitem)
                controlindex.append(newitemindex)
        if "control" not in newitem and "library" not in newitem and "ID" not in newitem and "_" not in newitem and "type" not in newitem: 
            if '.' not in newitem:
                tempitem=int(newitem)
            else:
                tempitem=float(newitem)
            totallist.append(tempitem)
    length=len(peptide)
    end=length-1
    end2=length-1
    print("controlindex "+str(controlindex))
    print("peptide[0] "+str(peptide[0]))
    print("peptide[1] "+str(peptide[1]))
    print("peptide[2] "+str(peptide[2]))
    print("countserumcontrol "+str(countserumcontrol))
    
    for newind in range(0,len(totallist)):
        peptidetype=peptide[newind/samplenumber]
        subpeptidetype2=peptidetype[len(peptidetype)-1]
        newstring=str(totallist[newind])+' '+subpeptidetype2
        if newstring=="65535 L":
            tempcount=tempcount+1
        newtotallist.append(newstring)

    for ind2 in range(0,samplenumber):
        ind3=ind2
        while ind3<len(totallist):
            newitem2=newtotallist[ind3]
            array.append(newitem2)
            ind3=ind3+samplenumber

    for ind4 in range(0,samplenumber):
        while ind5>start and ind5<end+1:
            newitem3=array[ind5][0:len(array[ind5])-2]
            arrayitem.append(newitem3)
            ind5=ind5+1
        arrayitem.sort()
        median=findmedian(arrayitem)
        arraymedian.append(int(median))
        start=start+length
        end=end+length
        arrayitem=[]
        print("array median "+str(arraymedian))
        
        while ind6>start2 and ind6<end2+1:
            newitem4=array[ind6][0:len(array[ind6])-2]
            if array[ind6]=="0 L":
                countzero=countzero+1
            elif array[ind6]=="65535 L":
                countmax=countmax+1
            ind6=ind6+1
        fractionzero=float(countzero/126051)
        fractionmax=str(countmax)+"/"+"126051"
        arraynumzero.append(fractionzero)
        arraynummax.append(fractionmax)
        start2=start2+length
        end2=end2+length
        countzero=0
        countmax=0

    while index3<len(data3):
        while data3[index3]!='\t'and data3[index3]!='\n' and index3<len(data3)-1:
            item3=item3+data3[index3]
            index3=index3+1
        if index3==len(data3)-1:
            item3=item3+data3[len(data3)-1]
        subitem3=item3[1:len(item3)-1]
        mylist2.append(subitem3)
        index3=index3+1
        item3=""
    
    for index4 in range(0,len(mylist2)):
        if "13-1024" in mylist2[index4]:
            sampletype=mylist2[index4+5]
            reference.append(sampletype)

    for index5 in range(0,samplenumber):
        subarraystart=index5*length
        if reference[index5]!="ambiguous":
            for index6 in range(subarraystart, subarraystart+length):
                subarray.append(array[index6][0:len(array[index6])-2])

    variationcoefficient=findvariationcoefficient(peptide,subarray,length)
    variationcoefficientnormalization=findvariationcoefficient(peptide,anomylist,length)
    print("variation coefficient for original serum control peptides "+str(variationcoefficient))
    print("variation coefficient for normalized serum control peptides "+str(variationcoefficientnormalization))

    print("controlindex "+str(controlindex))       
    for ind9 in range(0,126051):
        if ind9 in controlindex:
            start5=ind9*38-1
            end5=start5+38
            ind10=start5+1
            while ind10>start5 and ind10<end5+1:
                control.append(totallist[ind10])
                ind10=ind10+1
            totalcontrol.append(control)
            control=[]
            
    for ind11 in range(0,len(totalcontrol)-1):
        for ind12 in range(ind11+1,len(totalcontrol)):
            if ind11!=ind12:
                correlation=findcorrelation(totalcontrol[ind11],totalcontrol[ind12])
                controlpearsoncorrelation.append(correlation)
                controlpearsoncorrelationindex.append((ind11,ind12))
                if correlation>largest:
                    largest=correlation
                    largestindex1=controlindex[ind11]
                    largestindex2=controlindex[ind12]

    for ind13 in range(0,len(controlpearsoncorrelation)):
        if  controlpearsoncorrelation[ind13]==largest:
            print("control1 "+peptide[largestindex1])
            print("control2 "+peptide[largestindex2])
            print("correlation "+str(largest))
   
    while index7<len(data4):
        while data4[index7]!=' 'and index7<len(data4)-1:
            item4=item4+data4[index7]
            index7=index7+1
        if index7==len(data4)-1:
            item4=item4+data4[len(data4)-1]
        item4=float(item4)
        normalizedambiguous.append(item4)
        index7=index7+1
        item4=""

    while index8<len(data5):
        while data5[index8]!=' 'and index8<len(data5)-1:
            item5=item5+data5[index8]
            index8=index8+1
        if index8==len(data5)-1:
            item5=item5+data5[len(data5)-1]
        item5=float(item5)
        normalizednonambiguous.append(item5)
        index8=index8+1
        item5=""

    for index9 in range(0,len(normalizedambiguous)):
        normalizedambiguousnonambiguous.append(normalizedambiguous[index9])

    for index10 in range(0,len(normalizednonambiguous)):
        normalizedambiguousnonambiguous.append(normalizednonambiguous[index10])

    for index11 in range(0,126051):
        index12=index11
        while index12<len(normalizedambiguousnonambiguous):
            item6=normalizedambiguousnonambiguous[index12]
            newnormalizedambiguousnonambiguous.append(item6)
            index12=index12+126051

    for ind14 in range(0,126051):
        if ind14 in controlindex:
            start6=ind14*38-1
            end6=start6+38
            ind15=start6+1
            while ind15>start6 and ind15<end6+1:
                normalizedcontrol.append(newnormalizedambiguousnonambiguous[ind15])
                ind15=ind15+1
            totalnormalizedcontrol.append(normalizedcontrol)
            normalizedcontrol=[]

    for ind16 in range(0,len(totalcontrol)-1):
        for ind17 in range(ind16+1,len(totalnormalizedcontrol)):
            if ind16!=ind17:
                correlation2=findcorrelation(totalnormalizedcontrol[ind16],totalnormalizedcontrol[ind17])
                normalizedcontrolpearsoncorrelation.append(correlation2)
                normalizedcontrolpearsoncorrelationindex.append((ind16,ind17))
                if correlation2>largest2:
                    largest2=correlation2
                    newlargestindex1=controlindex[ind16]
                    newlargestindex2=controlindex[ind17]

    print("largest2 "+str(largest2))
    for ind18 in range(0,len(normalizedcontrolpearsoncorrelation)):
        if  normalizedcontrolpearsoncorrelation[ind18]==largest2:
            print("normalizedcontrol1 "+peptide[newlargestindex1])
            print("normazliedcontrol2 "+peptide[newlargestindex2])
            print("correlation2 "+str(largest2))
    
def findmedian(listl):
    if len(listl)%2==0:
        item1=listl[len(listl)/2-1]
        item2=listl[len(listl)/2]
        total=item1+item2
        median=float(total)/2
    else:
        median=listl[len(listl)/2]
    return median

def findcorrelation(list1,list2):
    list1total=0
    list2total=0
    totalsquarelist1item=0
    totalsquarelist2item=0
    totalproduct=0
    for index in range(0,len(list1)):
        list1total=list1total+list1[index]
        squarelist1item=math.pow(list1[index],2)
        totalsquarelist1item=totalsquarelist1item+squarelist1item
    
    for index2 in range(0,len(list2)):
        list2total=list2total+list2[index2]
        squarelist2item=math.pow(list2[index2],2)
        totalsquarelist2item=totalsquarelist2item+squarelist2item

    for index3 in range(0,len(list1)):
        product=list1[index3]*list2[index3]
        totalproduct=totalproduct+product
    top=len(list1)*totalproduct-list1total*list2total
    squarelist1total=math.pow(list1total,2)
    squarelist2total=math.pow(list2total,2)
    item1=math.sqrt(len(list1)*totalsquarelist1item-squarelist1total)
    item2=math.sqrt(len(list1)*totalsquarelist2item-squarelist2total)
    bottom=item1*item2
    correlation=top/bottom
    return correlation

def findvariationcoefficient(peptide,anomylist,length):
    ind7=0
    ind8=0
    arraytotal=0
    totalpowerdifference=0
    start3=-1
    end3=length-1
    countserumcontrol=42
    variationcoefficientlist=[]
    
    for newind2 in range(0,34):
        while ind7>start3 and ind7<end3+1:
            tempind=ind7%126051
            if "C" in peptide[tempind]:
                newitem5=anomylist[ind7]
                arraytotal=arraytotal+float(newitem5)
            ind7=ind7+1
        arrayaverage=arraytotal/countserumcontrol
        arraytotal=0
        
        while ind8>start3 and ind8<end3+1:
            tempind=ind8%126051
            if "C" in peptide[tempind]:
                difference=float(anomylist[ind8])-arrayaverage
                powerdifference=math.pow(difference,2)
                totalpowerdifference=totalpowerdifference+powerdifference
            ind8=ind8+1
        ratiototalpowerdifference=totalpowerdifference/countserumcontrol
        standarddeviation=math.sqrt(ratiototalpowerdifference)
        variationcoefficient=standarddeviation/arrayaverage
        variationcoefficientlist.append(variationcoefficient)
        start3=start3+length
        end3=end3+length
        totalpowerdifference=0
    return variationcoefficientlist

def findproduct(list1,list2):
    total=0
    for index in range(0,len(list1)):
        product=list1[index]*list2[index]
        total=total+product
    return total
        
    
    
