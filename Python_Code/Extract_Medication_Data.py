f_m = open("/Volumes/ATUL_6TB/Data/BioFinder/Medication/Medication_Data_AD_PD.txt", 'w', 1)
with open ("/Volumes/ATUL_6TB/Data/BioFinder/Medication/Medication_Data.txt", 'r') as myFile:
    line = myFile.readline ()
    f_m.write (line.rstrip() + "\t" + "AD_PD_Medication" + "\n")
    for line in myFile:
        if (("Memantin" in line) or ("memantin" in line)):
            f_m.write (line.rstrip() + "\t" + "1" + "\n")
        else:
            f_m.write (line.rstrip() + "\t" + "0" + "\n")
