# Fragment Analysis - Circ RNA Smear Analysis Script
#
# Version 2: Calculates specific range for introns, main, concatemers
#   introns = -10% of low intron to +10% of high intron
#   main = +/- 10% of main peak
#   concatemers = +/- 10% of concatemer peak
#
#   Operator needs to set threshold so that there are 2 intron peaks, 1 main peak, 1 concatemer peak ONLY
#
# 2023
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

    print(vals)
    capnum = 0
    for x in range(len(lines)-30):
        if lines[x][:10] == '[Capillary':
            capnum += 1
            if capnum in vals.keys():
                inLow = str(round(vals[capnum][0],7))
                inHigh = str(round(vals[capnum][1],7))
                MnLow = str(round(vals[capnum][2],7))
                MnHigh = str(round(vals[capnum][3],7))
                ConctLow = str(round(vals[capnum][4],7))
                ConctHigh = str(round(vals[capnum][5],7))
                
                # rough approximation for concatemers
                high2 = vals[capnum][2]*1.5
                h2 = str(round(high2,7))

                lines[x+27] = 'Smear Analysis = "'+ inLow + ',' + inHigh + ',' + MnLow + ','+ MnHigh + ',' + ConctLow + ',' + ConctHigh + ',0,0,0,0,0,0,0,0,0,0"\n'             

    f = open(anai, mode='w')
    f.writelines(lines)
    f.close()  
    return True
    


# for use by other programs - without calling main func
# NEED TO FIX FOR CIRC GUI APP
def inputs(csv,a):
    try: 
        c = open(csv, mode='r')
    except FileNotFoundError:
        return False
    data = []
    for row in csv:
        vals = row.split(',')
        data += [vals]

    # create array of well data [[well str, len, %], [well str, len, %], [well str, len, %],etc.] 
    peaks = []
    for x in data:
        peaks.append([x[0],x[2],x[3]])
    print(peaks)
    print('\n')

    # get boundaries between each well's peaks
    blks = []
    for y in peaks:
        if y[1] == '':
            blks.append(peaks.index(y))
      
    # create array of main peak lengths
    lengths = []
    for i in range(len(blks)):
        mp = 40.0
        if i > 0:
            for j in range(blks[i-1],blks[i]):
                if (peaks[j][2] != '') and (float(peaks[j][2]) >= float(mp)):
                    mp = peaks[j][1]
                    in1 = peaks[j-2][1]
                    in2 = peaks[j-1][1]
                    cnct = peaks[j+1][1]
                    lengths.append([j,in1,in2,mp,cnct])
        else:
            for j in range(1,blks[i]):
                if (peaks[j][2] != '') and (float(peaks[j][2]) >= float(mp)):
                    mp = peaks[j][1]
                    in1 = peaks[j-2][1]
                    in2 = peaks[j-1][1]
                    cnct = peaks[j+1][1]
                    lengths.append([j,in1,in2,mp,cnct])
    print(lengths,'\n')

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
        peaks.append([x[0],x[2],x[3]])
    print(peaks)
    print('\n')

    # get boundaries between each well's peaks
    blks = []
    for y in peaks:
        if y[1] == '':
            blks.append(peaks.index(y))
      
    # create array of main peak lengths
    lengths = []
    for i in range(len(blks)):
        mp = 40.0
        if i > 0:
            for j in range(blks[i-1],blks[i]):
                if (peaks[j][2] != '') and (float(peaks[j][2]) >= float(mp)):
                    mp = peaks[j][1]
                    in1 = peaks[j-2][1]
                    in2 = peaks[j-1][1]
                    cnct = peaks[j+1][1]
                    lengths.append([j,in1,in2,mp,cnct])
        else:
            for j in range(1,blks[i]):
                if (peaks[j][2] != '') and (float(peaks[j][2]) >= float(mp)):
                    mp = peaks[j][1]
                    in1 = peaks[j-2][1]
                    in2 = peaks[j-1][1]
                    cnct = peaks[j+1][1]
                    lengths.append([j,in1,in2,mp,cnct])
    print(lengths,'\n')

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
    
