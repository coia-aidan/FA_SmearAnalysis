# Fragment Analysis - Circ RNA Smear Analysis Script - 2
#
#   Operator needs to set threshold so that there are 2 peaks-> peak 1 = uncirc, peak 2 = circ
#
# 2024
#
# Aidan Coia - Strand Therapeutics
# 
#
#
#

capd = {'"A1"': 1, '"A2"': 2, '"A3"': 3, '"A4"': 4, '"A5"': 5, '"A6"': 6, '"A7"': 7, '"A8"': 8, '"A9"': 9, '"A10"': 10, '"A11"': 11, '"A12"': 12, 
            '"B1"': 13, '"B2"': 14, '"B3"': 15, '"B4"': 16, '"B5"': 17, '"B6"': 18, '"B7"': 19, '"B8"': 20, '"B9"': 21, '"B10"': 22, '"B11"': 23, '"B12"': 24, 
            '"C1"': 25, '"C2"': 26, '"C3"': 27, '"C4"': 28, '"C5"': 29, '"C6"': 30, '"C7"': 31, '"C8"': 32, '"C9"': 33, '"C10"': 34, '"C11"': 35, '"C12"': 36, 
            '"D1"': 37, '"D2"': 38, '"D3"': 39, '"D4"': 40, '"D5"': 41, '"D6"': 42, '"D7"': 43, '"D8"': 44, '"D9"': 45, '"D10"': 46, '"D11"': 47, '"D12"': 48, 
            '"E1"': 49, '"E2"': 50, '"E3"': 51, '"E4"': 52, '"E5"': 53, '"E6"': 54, '"E7"': 55, '"E8"': 56, '"E9"': 57, '"E10"': 58, '"E11"': 59, '"E12"': 60, 
            '"F1"': 61, '"F2"': 62, '"F3"': 63, '"F4"': 64, '"F5"': 65, '"F6"': 66, '"F7"': 67, '"F8"': 68, '"F9"': 69, '"F10"': 70, '"F11"': 71, '"F12"': 72, 
            '"G1"': 73, '"G2"': 74, '"G3"': 75, '"G4"': 76, '"G5"': 77, '"G6"': 78, '"G7"': 79, '"G8"': 80, '"G9"': 81, '"G10"': 82, '"G11"': 83, '"G12"': 84, 
            '"H1"': 85, '"H2"': 86, '"H3"': 87, '"H4"': 88, '"H5"': 89, '"H6"': 90, '"H7"': 91, '"H8"': 92, '"H9"': 93, '"H10"': 94, '"H11"': 95, '"H12"': 96}
    
# function: anai
#   takes anai file and peak values array
#   reads smear analysis lines, inserting +/-10% range values according to peak size for each sample
#   writes edited lines over original anai file, returns True
def anai(anai,vals):
    try:
        f = open(anai, mode='r')
    except FileNotFoundError:
        return False
    lines = f.readlines()

    capnum = 0
    for x in range(len(lines)-30):
        if lines[x][:10] == '[Capillary':
            capnum += 1
            if capnum in vals.keys():
                uLow = str(round(vals[capnum][0],7))
                uHigh = str(round(vals[capnum][1],7))
                cLow = str(round(vals[capnum][2],7))
                cHigh = str(round(vals[capnum][3],7))

                lines[x+27] = 'Smear Analysis = "'+ uLow + ',' + uHigh + ',' + cLow + ','+ cHigh + ',0,0,0,0,0,0,0,0,0,0,0,0"\n'             

    f = open(anai, mode='w')
    f.writelines(lines)
    f.close()  
    return True
    


# for use by other programs - without calling main func
def inputs(csv,a):
    try: 
        c = open(csv, mode='r')
    except FileNotFoundError:
        return False
    data = []
    for row in c:
        vals = row.split(',')
        data += [vals]

    # create array of well data
    peaks = []
    for x in data:
        peaks.append([x[0],x[3]])
    print(peaks)
    wellDict = {}
    cWell = '"A1"'
    for y in range(1,len(peaks)):
        if peaks[y][0] == cWell:
            if cWell not in wellDict.keys():
                wellDict[cWell] = [peaks[y][1]]
                print(peaks[y][1])
            else:
                vals = list(wellDict.pop(cWell))
                print(peaks[y][1])
                if peaks[y][1] != '':
                    vals.append(peaks[y][1])
                wellDict[cWell] = vals
                print(vals)
        else:
            cWell = peaks[y][0]

    print(wellDict.keys())
    tPass = {}
    for x in wellDict.keys():
        cap = capd[x]
        pks = []

       
        i = 0
        j = 0
        while i<2:
            if len(wellDict[x]) <= j or wellDict[x][j] == '':
                i += 1
                break
            tmp = int(str(wellDict[x][j]))
            print(tmp)
            j += 1
            if tmp >= 1000:
                pks.append(tmp)
                i += 1

    
        if len(pks) == 2:
            mid = (pks[0]+pks[1])/2
            uLow = pks[0]*0.9
            cHigh = pks[1]*1.1
            tPass[cap] = [uLow, mid, mid, cHigh]
        elif len(pks) == 1:
            low = pks[0]*0.9
            high = pks[0]*1.1
            tPass[cap] = [low, high, 0, 0]
        else:
            tPass[cap] = [0,0,0,0]

    if anai(a,tPass) == True:
        return True
 
    
    
    
    
    
# Main
def main():
    file =  input('Paste FA Peak Table file for analysis: ')   
    csv = open(file, mode="r",encoding='utf8')
    data = []
    for row in csv:
        vals = row.split(',')
        data += [vals]

    # create array of well data [[well str, len, %], [well str, len, %], [well str, len, %],etc.] 
    peaks = []
    for x in data:
        peaks.append([x[0],x[2]])
    print(peaks)
    print('\n')
    return True

    

# calculate smear %'s and append to peaks dict
    mps = {}
    for b in range(len(lengths)):
        i1=float(lengths[b][1])
        i2=float(lengths[b][2])
        mnp=float(lengths[b][3])
        cnc=float(lengths[b][4])
        inL = 0.9 * i1
        inH = 1.1 * i2
        ML = 0.9 * mnp
        MH = 1.1 * mnp
        CL = 0.9 * cnc
        CH = 1.1 * cnc
       # if (x[0] == '"E1"'):
        Temp = [j for j in peaks[lengths[b][0]][0]]
        if len(Temp)>=2:
            if len(Temp) == 3:
                TempS = str('"'+Temp[0]+Temp[1]+Temp[2]+'"')
            else:
                TempS = str('"'+Temp[0]+Temp[1]+'"')
           #     cap = capd[TempS]-48 ### -48 for second half of plate if two halves run on 48cap FA
           # else:
            cap = capd[TempS]
            mps[cap] = [inL, inH, ML, MH, CL, CH]

        
    # Call anai function for .anai file
    a = input('Paste .anai file corresponding to peak table: ')
    if anai(a,mps) == True:
        print('\nSuccessful Smear Analysis.')
    

    
if __name__ == "__main__":
    main()
    
